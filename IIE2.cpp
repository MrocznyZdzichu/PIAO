#include <ogx/Plugins/EasyPlugin.h>
#include <ogx/Data/Clouds/CloudHelpers.h>
#include <ogx/Data/Clouds/KNNSearchKernel.h>
#include <ogx/Color.h>
#include <ogx/Common.h>
#include <cmath>

using namespace ogx;
using namespace ogx::Data;



struct IIE2 : public ogx::Plugin::EasyMethod
{
	// parameters
	Data::ResourceID	m_node_id;
	Integer				neighborCount;	

	//methods

#pragma region methods
			//inputs:	normalsNeighboorhood - normal vectors of given neighborhood,
//			neighborsCount - size of neighborhood
//output:	vector which is mean of all normals in given neighboorhood
	Clouds::Point3D meanVector(const std::vector<Clouds::Point3D>& normalsNeighborhood,
		const Integer& neighborCount)
	{
		float sumX = 0;
		float sumY = 0;
		float sumZ = 0;

		for (auto& normal : normalsNeighborhood)
		{
			sumX += normal(0, 0);
			sumY += normal(1, 0);
			sumZ += normal(2, 0);
		}

		return Clouds::Point3D(sumX / neighborCount, sumY / neighborCount, sumZ / neighborCount);
	}


	//computes angular deviation basing on equation below:
	//u dot v = |u||v|*cos(x) -> x = arccos((u dot v)/(|u||v|))
	//inputs:	u, v: 3D vectors
	//output:	angularDeviation in rad
	inline float angularDeviation(const Clouds::Point3D& u, const Clouds::Point3D& v)
	{
		return acos((u.dot(v)) / (u.norm() * v.norm()));
	}

	//debug method which verboses 'count' first deviations
	inline void logDeviations(const std::vector<float>& deviations, const int& count)
	{
		int i = 0;
		for (auto& dev : deviations)
		{
			OGX_LINE.Format(ogx::Debug, L"%f", dev);
			i++;
			if (i == count)
				return;
		}
	}

	//debug method which verboses normals, dot, norms and deviation of given vector pair
	inline void debugLog(Clouds::Point3D u, Clouds::Point3D v)
	{
		OGX_LINE.Format(ogx::Debug, L" u - x:%f, y:%f, z:%f", u(0, 0), u(1, 0), u(2, 0));
		OGX_LINE.Format(ogx::Debug, L" v - x:%f, y:%f, z:%f", v(0, 0), v(1, 0), v(2, 0));
		OGX_LINE.Format(ogx::Debug, L"u dot v - %f", u.dot(v));
		OGX_LINE.Format(ogx::Debug, L"|u|: %f", u.norm());
		OGX_LINE.Format(ogx::Debug, L"|V|: %f", v.norm());
		OGX_LINE.Format(ogx::Debug, L"phi: %f", angularDeviation(u, v));
	}
#pragma endregion

#pragma region params,constructor,init
	// constructor
	IIE2() : EasyMethod(L"Maksym Zawrotny", L"Computes angular deviation from mean normal vector")
	{
	}

	// add input/output parameters
	virtual void DefineParameters(ParameterBank& bank)
	{
		bank.Add(L"node_id", m_node_id).AsNode();
		bank.Add(L"Neighbors count", neighborCount = 3).Min(3);
	}

	bool Init(Execution::Context& context)
	{
		OGX_SCOPE(log);
		OGX_LINE.Msg(User, L"Initialization succeeded");
		return EasyMethod::Init(context);
	}
#pragma endregion

	virtual void Run(Context& context)
	{
		auto subtree = context.Project().TransTreeFindNode(m_node_id);
		// report error if give node was not found, this will stop execution of algorithm
		if (!subtree) ReportError(L"Node not found");

		// run with number of threads available on current machine, optional
		auto const thread_count = std::thread::hardware_concurrency();

		//prepare kernel
		Clouds::KNNSearchKernel knnKernel(Math::Point3D::Zero(), neighborCount);

		// perform calculations for each cloud in given subtree
		Clouds::ForEachCloud(*subtree, [&](Data::Clouds::ICloud & cloud, Data::Nodes::ITransTreeNode & node)
		{
			// access points in the cloud
			Clouds::PointsRange points_all;
			cloud.GetAccess().GetAllPoints(points_all);

			//get normal vectors into prereserved container
			std::vector<Clouds::Point3D> normalVectors;
			normalVectors.reserve(points_all.size());

#pragma region computingNormalVectors
			//somehow I couldnt pass cloud as argument so it is not separate method :(
			//define neighborhood seach kernel
			Clouds::KNNSearchKernel knnKernelNormals(Math::Point3D::Zero(), 3);

			for (auto& xyz : Data::Clouds::RangeLocalXYZConst(points_all))
			{
				//find neighborhood
				Clouds::PointsRange neighborhood;
				knnKernel.GetPoint() = xyz.cast<Math::Point3D::Scalar>();
				cloud.GetAccess().FindPoints(knnKernel, neighborhood);

				//get xyz coords of neighbors
				std::vector<Clouds::Point3D> neighboor_xyz;
				neighboor_xyz.reserve(3);
				neighborhood.GetXYZ(neighboor_xyz);

				//fit best 3D plane and get its normal vector
				Math::Point3D fittedPlaneNormal = Math::CalcBestPlane3D(neighboor_xyz.begin(), neighboor_xyz.end()).normal();
				normalVectors.push_back(fittedPlaneNormal.cast<StoredPoint3D::Scalar>());
			}
			points_all.SetNormals(normalVectors);

#pragma endregion

#pragma region computingDeviations
			//somehow I couldnt pass cloud as argument so it is not separate method :(
			//reserve vector for meanVectors
			std::vector<Clouds::Point3D> meanVectors;
			meanVectors.reserve(points_all.size());

			//reserve vector for output deviations
			std::vector<float> deviations;
			deviations.reserve(points_all.size());

			//compute deviations for all points' neighborhoods
			int pointNumber = 0;
			for (auto& xyz : Data::Clouds::RangeLocalXYZConst(points_all))
			{
				//prepare containers for neighborhood and their normals
				Clouds::PointsRange neighborhood;
				std::vector<Clouds::Point3D> normalsNeighbors;
				normalsNeighbors.reserve(neighborCount);

				//search for neighbors and their normal vectors
				knnKernel.GetPoint() = xyz.cast<Math::Point3D::Scalar>();
				cloud.GetAccess().FindPoints(knnKernel, neighborhood);
				neighborhood.GetNormals(normalsNeighbors);

				//compute mean normal vector and deviations				
				meanVectors.push_back(meanVector(normalsNeighbors, neighborCount));
				deviations.push_back(angularDeviation(meanVectors[pointNumber], normalVectors[pointNumber]));

				//debug verboses
				if (!pointNumber)
					debugLog(meanVectors[0], normalVectors[0]);
				pointNumber++;
			}
			logDeviations(deviations, 100);
#pragma endregion

#pragma region addDeviationsToLayer
			// declare layer name and interface to layer
			String layerName = L"Deviations layer";
			Data::Layers::ILayer *deviationsLayer;

			// check if layer with that name exit already. If - points it, else create it
			auto layers = cloud.FindLayers(layerName);
			!layers.empty() ? deviationsLayer = layers[0] :
				deviationsLayer = cloud.CreateLayer(layerName, 0);

			points_all.SetLayerVals(deviations, *deviationsLayer);
#pragma endregion

		}		
		, thread_count); // run with given number of threads, optional parameter, if not given will run in current thread
	}
};

OGX_EXPORT_METHOD(IIE2)
