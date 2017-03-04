#include "dgemm.hpp"
#include <iostream>
#include <random>
#include "timer.hpp"
#include "flops_calc.hpp"
#include <random>
#include <ctime>
#include <cblas.h>

static void* allocate(size_t size)
{
    return malloc(size);
}

static void destruct(void* ptr)
{
    free(ptr);
}

static void dgemm(const int m,
		  const int n,
		  const int k,
		  const double alpha,
		  const double* A,
		  const int lda,
		  const double* B,
		  const int ldb,
		  const double beta,
		  double* C,
		  const int ldc)
{
    cblas_dgemm(CblasRowMajor,
		CblasNoTrans,
		CblasNoTrans,
		m, n, k,
		alpha,
		A, lda,
		B, ldb,
		beta,
		C, ldc);
}


void exec_dgemm(int n, int m, int k, double alpha, double beta, double time_limit)
{
    std::cout << "Dgemm start" << std::endl;

    int lda = k;
    int ldb = n;
    int ldc = n;

    double* pA = reinterpret_cast<double*>(allocate(sizeof(double) * lda * m));
    double* pB = reinterpret_cast<double*>(allocate(sizeof(double) * ldb * k));
    double* pC = reinterpret_cast<double*>(allocate(sizeof(double) * ldc * m));

    size_t memory_use = sizeof(double) * lda * m;
    memory_use += sizeof(double) * ldb * k;
    memory_use += sizeof(double) * ldc * m;

    std::cout << "memory use " << memory_use / 1024.0 / 1024.0 / 1024.0 << " GB" << std::endl;

    {
	std::mt19937 mt(1);
	std::uniform_real_distribution<double> rand01(0, 1);
	// initialize
	for(size_t i = 0; i < lda * m; ++i) {
	    pA[i] = rand01(mt);
	}

	for(size_t i = 0; i < ldb * k; ++i) {
	    pB[i] = rand01(mt);
	}

	for(size_t i = 0; i < ldc * m; ++i) {
	    pC[i] = 0;
	}
    }

    {
	time_t now = time(nullptr);
	struct tm* pNow = localtime(&now);
	std::cout << pNow->tm_hour << ":" << pNow->tm_min << ":" << pNow->tm_sec << " initialize done." << std::endl;
    }

    FlopsCalc fp(m, n, k);
    double elap = 0.0;
    Timer timer;
    while(elap < time_limit) {
	timer.Start();
	dgemm(m, n, k, alpha, pA, lda, pB, ldb, beta, pC, ldc);
	double s = timer.getSec();
	elap += s;
	fp.setFlops(s);
    }

    double average, s;
    fp.getAverage(average, s);
    std::cout << average << " GFlops, " << s << std::endl;
//    std::cout << "elap " << elap << " time_limit " << time_limit << std::endl;

    destruct(pA);
    destruct(pB);
    destruct(pC);

    return ;
}
