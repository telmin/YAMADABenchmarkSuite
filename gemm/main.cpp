#include <iostream>
#include "dgemm.hpp"
#include "sgemm.hpp"
#include <omp.h>

int main(int argc, char** argv)
{
    // 300秒実行する
    double time_limit = 30.0;

//    int n = 40960;
//    int m = 40960;
//    int k = 40960;
    int n = 10240;
    int m = 10240;
    int k = 10240;

    double alpha = -1.0;
    double beta = 1.0;

    std::cout << "M " << m << " N " << n << " K " << k << " al " << alpha << " b " << beta << std::endl;

    //MKLVersion mkl_version;
    //mkl_get_version(&mkl_version);
    //std::cout << "MKL " << mkl_version.MajorVersion << "." << mkl_version.MinorVersion << std::endl;

    exec_dgemm(n, m, k, alpha, beta, time_limit);
    exec_sgemm(n, m, k, static_cast<float>(alpha), static_cast<float>(beta), time_limit);
}
