#ifndef TIMER_HPP
#define TIMER_HPP

#include <chrono>

class Timer
{
public:
    Timer()
    {}

    void Start() {
	start = std::chrono::system_clock::now();
    }
    double Stop() {
	auto end = std::chrono::system_clock::now();

	double elap = static_cast<double>(std::chrono::duration_cast<std::chrono::microseconds>(end-start).count());
	return elap / 1000.0;
    }

    double getSec() {
	double ms = Stop();
	return (ms / 1000.0);
    }

private:
    std::chrono::system_clock::time_point start;
};


#endif
