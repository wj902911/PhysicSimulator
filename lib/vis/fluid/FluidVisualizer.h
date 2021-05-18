#pragma once

#include "../IglVisualizer.h"
#include <thread>
#include <memory>

static PhysicsHook* hook = NULL;

namespace fluid {

	class FluidCore;
	class FluidHook;

	class FluidVisualizer :public IglVisualizer {
		std::unique_ptr<FluidCore> core_;
		std::unique_ptr<FluidHook> hook_;
	public:
		FluidVisualizer();
		~FluidVisualizer();
	};
}