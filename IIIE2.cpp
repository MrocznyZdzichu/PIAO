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
	float				maxDist;
	float				neighSize;
#pragma region params,constructor,init
	// constructor
	IIIE2() : EasyMethod(L"Maksym Zawrotny", L"Show planes and compute thier area")
	{
	}

	float PlanePointDist(Clouds::Point3D normal, Data::Clouds::Point3D point, float offset)
	{
		return fabs(normal(0, 0) * point(0, 0) +
			normal(1, 0) * point(1, 0) +
			normal(2, 0) * point(2, 0) +
			offset);
	}
	// add input/output parameters
	virtual void DefineParameters(ParameterBank& bank)
	{
		bank.Add(L"node_id", m_node_id).AsNode();
		bank.Add(L"neighbors count", neighSize);
		bank.Add(L"max point distance to fir plane", maxDist);
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

			ogx::StoredReal		groupNumber = 1;
			int					iterations = 0;
			//prepare iterators and ranges to sweep through pc
			auto xyzRange = Clouds::RangeLocalXYZConst(points_all);
			auto segNumRange = Clouds::RangeLayer(points_all, *segmentsLayer);
			auto xyz = xyzRange.begin();
			auto segNum = segNumRange.begin();
		
#pragma region SegmentPointCloudIntoPlanes
			for (; xyz != xyzRange.end(); xyz++, segNum++)
			{
				// go through all points without plane-segment yet
				// segNum is reference to float so to avoid float == i am doing it this way
				
				if (fabs(*segNum) < 0.5)
				{
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

					for (auto neighXYZ : neighbor_xyz)
					{
						currentDistance = PlanePointDist(normal, *xyz, offset);
						currentDistance > maxPointDistance
							? maxPointDistance = currentDistance
							: maxPointDistance = maxPointDistance;
					}

					//if neighborhood is valid and big enough, assign new semgent number
					if (maxPointDistance <= maxDist/100)
					{
						std::vector<ogx::StoredReal> segNumbers;

						for (int i = 0; i < neighSize; i++)
							segNumbers.push_back(groupNumber);

						neighborhood.SetLayerVals(segNumbers, *segmentsLayer);
						groupNumber++;
					}
				}
			}
#pragma endregion

			//go back with segment number iterator and prepare state iterator
			auto segNumRange2 = Clouds::RangeLayer(points_all, *segmentsLayer);
			auto segNum2 = segNumRange2.begin();
			auto statesRange = Clouds::RangeState(points_all);
			auto state = statesRange.begin();

			// select all points segmented as planes
			for (; segNum2 != segNumRange2.end(); segNum2++, state++)
				if (fabs(*segNum2) > 0.5)
					state->set(Clouds::PS_SELECTED);
		}
		, thread_count); // run with given number of threads, optional parameter, if not given will run in current thread
	}
};

OGX_EXPORT_METHOD(IIIE2)
