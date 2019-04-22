#include <ogx/Plugins/EasyPlugin.h>
#include <ogx/Data/Clouds/CloudHelpers.h>
#include <ogx/Data/Clouds/KNNSearchKernel.h>
#include <ogx/Color.h>
#include <ogx/Common.h>

using namespace ogx;
using namespace ogx::Data;

struct IE4 : public ogx::Plugin::EasyMethod
{
	// parameters
	Data::ResourceID	m_node_id;
	Integer				m_neighbors_count;
	String				m_error_layer_name;
	float				Rmax, Rmin,
						Gmax, Gmin,
						Bmax, Bmin;
	// constructor
	IE4() : EasyMethod(L"Maksym Zawrotny", L"Filters PC that only matching colour remains.")
	{
	}

	// add input/output parameters
	virtual void DefineParameters(ParameterBank& bank)
	{
		bank.Add(L"node_id", m_node_id).AsNode();
		bank.Add(L"Max red value", Rmax = 255).Max(255);
		bank.Add(L"Min red value", Rmin = 0).Min(0);
		bank.Add(L"Max green value", Gmax = 255).Max(255);
		bank.Add(L"Min green value", Gmin = 0).Min(0);
		bank.Add(L"Max blue value", Bmax = 255).Max(255);
		bank.Add(L"Min blue value", Bmin = 0).Min(0);
	}

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
			Data::Clouds::PointsRange points_all;
			cloud.GetAccess().GetAllPoints(points_all);
		
			// initialize buffer for the output states
			std::vector<ogx::Data::Clouds::State> outputStateBuffer;
			ogx::Data::Clouds::State currentState;
			// reserve points in advance to speed things up
			outputStateBuffer.reserve(points_all.size());

			auto const ignored_mask = context.Feedback().GetIgnoredChannels();
			auto const selected_mask = context.Feedback().GetActiveChannels();
			// setup variables for color computations
			ogx::Byte b, g, r;

			// get access to range of colors
			auto pointColors = Data::Clouds::RangeColor(points_all);
			// iterate through each Color and check filtering conditions
			for (ogx::Color color : pointColors)
			{	
				r = color(0, 0);
				g = color(1, 0);
				b = color(2, 0);

				// checking filtering conditions and change state depending on color values
				if (((r < Rmax && r > Rmin) && 
					(g < Gmax && g > Gmin) && 
					(b < Bmax && b > Bmin)))
					currentState |= selected_mask;

				// add current state to output state buffer
				outputStateBuffer.push_back(currentState);
				// reset 
				currentState.reset();
			}
			// set respective states
			points_all.SetStates(outputStateBuffer);
		}, thread_count); // run with given number of threads, optional parameter, if not given will run in current thread
	}
};

OGX_EXPORT_METHOD(IE4)
