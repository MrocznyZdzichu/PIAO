#include <ogx/Plugins/EasyPlugin.h>
#include <ogx/Data/Clouds/CloudHelpers.h>
#include <ogx/Data/Clouds/KNNSearchKernel.h>
#include <ogx/Color.h>
#include <ogx/Common.h>
#include <cmath>
#include <vector>

using namespace ogx;
using namespace ogx::Data;

struct IIIE2 : public ogx::Plugin::EasyMethod
{
	// parameters
	Data::ResourceID	m_node_id;
	float neighSize = 20;
	float maxDist = 0.2;
#pragma region params,constructor,init
	// constructor
	IIIE2() : EasyMethod(L"Maksym Zawrotny", L"Show planes and compute thier area")
	{
	}

	// method computing point's distance to plane: normal*xyz + offset
	float PlanePointDist(Clouds::Point3D normal, Data::Clouds::Point3D point, float offset)
	{
		return fabs(normal(0, 0) * point(0, 0) +
			normal(1, 0) * point(1, 0) +
			normal(2, 0) * point(2, 0) +
			offset);
	}

	Eigen::Matrix<float, 3, 3> concVector(Clouds::Point3D u,
										  Clouds::Point3D v,
										  Clouds::Point3D w)
	{
		Eigen::Matrix<float, 3, 3> originBase;
		originBase(0, 0) = u(0, 0);
		originBase(0, 1) = v(0, 0);
		originBase(0, 2) = w(0, 0);
		originBase(1, 0) = u(1, 0);
		originBase(1, 1) = v(1, 0);
		originBase(1, 2) = w(1, 0);
		originBase(2, 0) = u(2, 0);
		originBase(2, 1) = v(2, 0);
		originBase(2, 2) = w(2, 0);

		return originBase;
	}

	Clouds::Point3D transformPoint(Clouds::Point3D point,
								   Clouds::Point3D u,
								   Clouds::Point3D v,
								   Clouds::Point3D w)
	{
		return  concVector(u, v, w).inverse() * point;
	}
	// method computing plane approx area
	float planeArea(std::vector<Clouds::Point3D> plane)
	{
		// compute new coords' vuw (w is normal and uv lies in plane)
		auto fittedPlane = Math::CalcBestPlane3D(plane.begin(), plane.end());
		Clouds::Point3D normal = fittedPlane.normal().cast<StoredPoint3D::Scalar>();
		Clouds::Point3D u = Clouds::Point3D(0.33, 0.66, 0.66);
		Clouds::Point3D v = u.cross(normal);
	

		// transform points xyz into new line-space base
		for (int i = 0; i < plane.size(); i++)
		{
			Clouds::Point3D newPoint = transformPoint(plane[i], u, v, normal);
			plane[i] = newPoint;
		}

		// compute range along u and v axis
		float minU = plane[0](0, 0);
		float maxU = plane[0](0, 0);
		float minV = plane[0](1, 0);
		float maxV = plane[0](1, 0);

		for (const auto& point : plane)
		{
			if (point(0, 0) < minU)
				minU = point(0, 0);

			if (point(1, 0) < minV)
				minV = point(1, 0);

			if (point(0, 0) > maxU)
				maxU = point(0, 0);

			if (point(1, 0) > maxV)
				maxV = point(1, 0);
		}

		float DU = maxU - minU;
		float DV = maxV - minV;

		// compute plane area as a mean of areas of cicrumscribed retangle area and inscribed elipse
		// approx method but in keeps performance quite good
		return (PI * 0.25 * DU * DV + DU * DV) / 2;
	}
	// add input/output parameters
	virtual void DefineParameters(ParameterBank& bank)
	{
		bank.Add(L"node_id", m_node_id).AsNode();
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

		// perform calculations for each cloud in given subtree
		Clouds::ForEachCloud(*subtree, [&](Data::Clouds::ICloud & cloud, Data::Nodes::ITransTreeNode & node)
		{
			// access points in the cloud
			Clouds::PointsRange points_all;
			cloud.GetAccess().GetAllPoints(points_all);;

			// create and 0-initialise layer holding segment's number
#pragma region Creating layer for segment numbers
			String layerName1 = L"Segment number layer";
			Data::Layers::ILayer *segmentsLayer;

			// check if layer with that name exit already. If - points it, else create it
			auto layers = cloud.FindLayers(layerName1);
			!layers.empty() ? segmentsLayer = layers[0] :
				segmentsLayer = cloud.CreateLayer(layerName1, 0);

			// initialize segment numbers with zeros
			std::vector<StoredReal> segmentNumbers;
			for (int i = 0; i < points_all.size(); i++)
				segmentNumbers.push_back(0);

			points_all.SetLayerVals(segmentNumbers, *segmentsLayer);
#pragma endregion

#pragma region SegmentPointCloudIntoPlanes
			ogx::StoredReal		groupNumber = 1;
			//prepare iterators and ranges to sweep through pc
			auto xyzRange = Clouds::RangeLocalXYZConst(points_all);
			auto xyz = xyzRange.begin();				
			auto segNumRange = Clouds::RangeLayer(points_all, *segmentsLayer);
			auto segNum = segNumRange.begin();

			for (; xyz != xyzRange.end(); xyz++, segNum++)
			{

				// I tried segmentation by growth but updating segNumRange killed performance :(			
				// define KNN kernel and find neighbors
				Clouds::KNNSearchKernel knnKernel(Math::Point3D::Zero(), neighSize);
				Clouds::PointsRange neighborhood;
				knnKernel.GetPoint() = (*xyz).cast<Math::Point3D::Scalar>();
				cloud.GetAccess().FindPoints(knnKernel, neighborhood);

				//get xyz coords of neighbors
				Data::Clouds::RangeLocalXYZConst neighbor_xyz(neighborhood);

				//fit plane 
				auto fittedPlane = Math::CalcBestPlane3D(neighbor_xyz.begin(), neighbor_xyz.end());
				Clouds::Point3D normal = fittedPlane.normal().cast<StoredPoint3D::Scalar>();
				float offset = fittedPlane.offset();

				//compute max point distance
				ogx::Real maxPointDistance = 0;
				ogx::Real currentDistance;

				//plane is fit so its distance to each point is equal :( so I break after one iteration
				for (auto neighXYZ : neighbor_xyz)
				{
					currentDistance = PlanePointDist(normal, *xyz, offset);
					currentDistance > maxPointDistance
						? maxPointDistance = currentDistance
						: maxPointDistance = maxPointDistance;
					break;
				}

				//if neighborhood is valid and big enough, assign new semgent number
				if (maxPointDistance <= maxDist / 100)
				{
					std::vector<ogx::StoredReal> segNumbers;

					for (int i = 0; i < neighSize; i++)
						segNumbers.push_back(groupNumber);

					neighborhood.SetLayerVals(segNumbers, *segmentsLayer);
					groupNumber++;
				}
			}

#pragma endregion

#pragma region SelectingPlanesOnGUI
			//go back with segment number iterator and prepare state iterator
			auto segNumRange2 = Clouds::RangeLayer(points_all, *segmentsLayer);
			auto segNum2 = segNumRange2.begin();
			auto statesRange = Clouds::RangeState(points_all);
			auto state = statesRange.begin();

			// select all points segmented as planes
			for (; segNum2 != segNumRange2.end(); segNum2++, state++)
				if (fabs(*segNum2) > 0.5)
					state->set(Clouds::PS_SELECTED);
#pragma endregion

#pragma region ComputingFlatSurfacesArea
			// create vectors for each plane
			std::vector<std::vector<Clouds::Point3D> > planePoints;
			planePoints.resize(groupNumber - 1);

			segNum2 = segNumRange2.begin();
			xyz = xyzRange.begin();

			// get points assigned to each plane
			for (; xyz != xyzRange.end(); segNum2++, xyz++)
			{
				if (fabs(*segNum2 - 0) > 0.5)
				{
					planePoints[int(std::round(*segNum2-1))].push_back(*xyz);
				}
			}

			float totalArea = 0;
			// count each plane area and add it to totalArea
			for (int i = 0; i < groupNumber - 1; i++)
			{
				if (planePoints[i].size() >= 3)
					totalArea += planeArea(planePoints[i]);
			}
			OGX_LINE.Format(ogx::Info, L"pole obszarów p³askich wynosi: %f", totalArea);
#pragma endregion

		}
		, thread_count); // run with given number of threads, optional parameter, if not given will run in current thread
	}
};

OGX_EXPORT_METHOD(IIIE2)
