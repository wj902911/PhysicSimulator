#include "FluidHook.h"
#include "FluidVisualizer.h"
#include <core/fluid/FluidCore.h>

namespace fluid {

	FluidVisualizer::FluidVisualizer()
	{
		core_.reset(new FluidCore);
		hook_.reset(new FluidHook(core_.get()));
		hook_->reset();
		init(core_.get(), hook_.get());
	}

	FluidVisualizer::~FluidVisualizer()
	{
	}

}