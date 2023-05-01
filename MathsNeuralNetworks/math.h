#ifndef MATH_H
#define MATH_H
#include <chrono>
//#define ONLY_POSITIVE
#define HALF_RAND_MAX RAND_MAX / 2

using frac = float;

namespace mcore { //math core
	struct Timer {
		std::chrono::time_point<std::chrono::high_resolution_clock> point;

		Timer() {
			reset();
		}

		void reset() {
			point = std::chrono::high_resolution_clock::now();
		}

		double get_millis() {
			return std::chrono::duration<double, std::milli>(std::chrono::high_resolution_clock::now() - point).count();
		}
	};

	frac get_random() {
#ifdef POSITIVE
		return frac(rand()) / RAND_MAX;
#else
		return frac(rand() - HALF_RAND_MAX) / HALF_RAND_MAX;
#endif
	}

	template<size_t AZ>
	void fill_with_random(frac (&arr)[AZ]) {
		for (size_t i = 0; i < AZ; ++i)
			arr = get_random();
	}
}
#endif