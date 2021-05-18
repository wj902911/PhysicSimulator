#include "FluidHook.h"
#include <core/fluid/FluidCore.h>

namespace fluid {
	FluidHook::FluidHook(FluidCore* core)
		: core_(core), PhysicsHook(core), params_(*core->getPointerToSimParameters())
	{
	}

	FluidHook::~FluidHook()
	{
	}

	void FluidHook::drawGUI(igl::opengl::glfw::imgui::ImGuiMenu& menu)
	{
	}

	void FluidHook::mouseClicked(double x, double y, int button)
	{
	}

	void FluidHook::tick()
	{
	}
}