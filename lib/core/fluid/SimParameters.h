#pragma once

namespace fluid {

	struct SimParameters {
		float airPressure = 0;
		float density = 1000;
		int cellNumber = 20;
		float gravity = -9.8;
		float timeStep = 0.001;
	};
}