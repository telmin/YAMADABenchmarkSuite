#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <fftw3.h>
#include <omp.h>

void test3d(int nx, int ny, int nz) {
    fftw_complex *in, *out;
    fftw_plan pf, pb;

    int omp_num;
#pragma omp parallel
    omp_num = omp_get_num_threads();

    fftw_plan_with_nthreads(omp_num);

    int n = nx * ny * nz;

    in = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * n);
    out = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * n);

    pf = fftw_plan_dft_3d(nx, ny, nz, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
    pb = fftw_plan_dft_3d(nx, ny, nz, in, out, FFTW_BACKWARD, FFTW_ESTIMATE);

    double t1 = omp_get_wtime();
    for(int i = 0; i < 10; ++i) {
	fftw_execute(pf);
    }
    double t2 = omp_get_wtime();
    for(int i = 0; i < 10; ++i) {
	fftw_execute(pb);
    }
    double t3 = omp_get_wtime();

    double gflops = 5.0 * (double)nx * (double)ny * (double)nz *
		    (log((double)nx) + log((double)ny) + log((double)nz))
		    / log(2.0) / 1e9;

    printf("3-D FFT: %d x %d x %d\n", nx, ny, nz);
    double t = (t2 - t1) / 10.0;
    printf("On-board: %f msec, %lf GFLOPS.\n", t, gflops / t);
    t = (t3 - t2) / 10.0;
    printf("On-board: %f msec, %lf GFLOPS.\n", t, gflops / t);

    fftw_destroy_plan(pf);
    fftw_destroy_plan(pb);
    fftw_free(in);
    fftw_free(out);
}

int main(int argc, char **argv) {
    test3d(128, 128, 128);
    test3d(256, 256, 256);
}
