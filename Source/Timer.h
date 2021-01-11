#pragma once
#include <chrono>
#include <iostream>
class Timer
{
private:
		std::chrono::time_point<std::chrono::high_resolution_clock> start;
public:
		Timer() {
				start = std::chrono::high_resolution_clock::now();
		}
		auto operator()() {
				return start;
		}
		~Timer() {
				std::chrono::time_point<std::chrono::high_resolution_clock> end
						= std::chrono::high_resolution_clock::now();
				std::chrono::duration<double>timeElapsed = end - start;
				std::cout << "The timeElapsed: " << timeElapsed.count() << std::endl;
		}
};

