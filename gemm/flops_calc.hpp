#ifndef FLOPS_CALC_HPP
#define FLOPS_CALC_HPP

#include <vector>

class FlopsCalc
{
public:
    FlopsCalc(int m, int n, int k)
	: m_m(m), m_n(n), m_k(k), m_recorded_flops(0)
    {}

    void setFlops(double sec) {
	double flops = (2.0 * m_m * m_n * m_k) / (1000.0 * 1000.0 * 1000.0 * sec);
	//std::cout << sec << " sec " << flops << " GFlops. " << std::endl;
	m_recorded_flops.push_back(flops);
    }

    const std::vector<double>& getFlops(void) {
	return m_recorded_flops;
    }

    void getAverage(double& average, double& standard) const {
	double sum = 0.0;
	double sum2 = 0.0;
	for(const auto& f : m_recorded_flops) {
	    sum += f;
	    sum2 += f * f;
	}

	average = sum / static_cast<double>(m_recorded_flops.size());
	double r = sum2 / static_cast<double>(m_recorded_flops.size()) - (average * average);
	if(r < 0.0) r = 0;
	standard = sqrt(r);
    }

private:
    int m_m;
    int m_n;
    int m_k;

    std::vector<double> m_recorded_flops;
};


#endif
