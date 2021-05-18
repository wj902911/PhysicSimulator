#pragma once

#include "../PhysicsHook.h"
#include <core/fluid/SimParameters.h>
#include <deque>

namespace fluid {
	class FluidCore;

	struct MouseClick
	{
		double x;
		double y;
		//SimParameters::ClickMode mode;
	};

	class FluidHook : public PhysicsHook {
	public:
		FluidHook(FluidCore* core);
		~FluidHook();

		virtual void drawGUI(igl::opengl::glfw::imgui::ImGuiMenu& menu) override;

		virtual void mouseClicked(double x, double y, int button) override;

		virtual void tick();

	private:
		FluidCore* core_;
		SimParameters& params_;
		double time_;

		std::mutex message_mutex;
		std::deque<MouseClick> mouseClicks_;
	};
}