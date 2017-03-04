/*-----------------------------------------------------------------------*/
/* Program: STREAM							 */
/* Revision: $Id: stream.c,v 5.10 2013/01/17 16:01:06 mccalpin Exp mccalpin $ */
/* Original code developed by John D. McCalpin				 */
/* Programmers: John D. McCalpin					 */
/*		Joe R. Zagar						 */
/*									 */
/* This program measures memory transfer rates in MB/s for simple	 */
/* computational kernels coded in C.					 */
/*-----------------------------------------------------------------------*/
/* Copyright 1991-2013: John D. McCalpin				 */
/*-----------------------------------------------------------------------*/
/* License:								 */
/*  1. You are free to use this program and/or to redistribute		 */
/*     this program.							 */
/*  2. You are free to modify this program for your own use,		 */
/*     including commercial use, subject to the publication		 */
/*     restrictions in item 3.						 */
/*  3. You are free to publish results obtained from running this	 */
/*     program, or from works that you derive from this program,	 */
/*     with the following limitations:					 */
/*     3a. In order to be referred to as "STREAM benchmark results",	 */
/*	   published results must be in conformance to the STREAM	 */
/*	   Run Rules, (briefly reviewed below) published at		 */
/*	   http://www.cs.virginia.edu/stream/ref.html			 */
/*	   and incorporated herein by reference.			 */
/*	   As the copyright holder, John McCalpin retains the		 */
/*	   right to determine conformity with the Run Rules.		 */
/*     3b. Results based on modified source code or on runs not in	 */
/*	   accordance with the STREAM Run Rules must be clearly		 */
/*	   labelled whenever they are published.  Examples of		 */
/*	   proper labelling include:					 */
/*	     "tuned STREAM benchmark results"				 */
/*	     "based on a variant of the STREAM benchmark code"		 */
/*	   Other comparable, clear, and reasonable labelling is		 */
/*	   acceptable.							 */
/*     3c. Submission of results to the STREAM benchmark web site	 */
/*	   is encouraged, but not required.				 */
/*  4. Use of this program or creation of derived works based on this	 */
/*     program constitutes acceptance of these licensing restrictions.	 */
/*  5. Absolutely no warranty is expressed or implied.			 */
/*-----------------------------------------------------------------------*/

#include <iostream>
#include <limits>
#include <string>
#include <algorithm>
#include <vector>
#include <chrono>
#include <iomanip>
#include <cstdio>

#include "parsearg.hpp"

static void PrintMessages(size_t STREAM_ARRAY_SIZE, size_t NTIMES, size_t OFFSET, size_t SizeOfStream)
{
    const std::string HLINE = "-------------------------------------------------------------";

    std::cout << HLINE << std::endl;
    std::cout << "STREAM version $Revision : 5.10 $" << std::endl;
    std::cout << HLINE << std::endl;

    std::cout << "This system uses " << SizeOfStream << " bytes per array element." << std::endl;
    std::cout << HLINE << std::endl;

    std::cout << "Array Size = " << STREAM_ARRAY_SIZE << " (elements), Offset = " << OFFSET << " (elements) " << std::endl;
    std::cout << "Memory per array = " << SizeOfStream * ((double)STREAM_ARRAY_SIZE / 1024.0 / 1024.0) << " MiB (= " <<
	SizeOfStream * ((double)STREAM_ARRAY_SIZE / 1024.0 / 1024.0 / 1024.0) << " GiB)." << std::endl;
    std::cout << "Total Memory required = " << (3.0 * SizeOfStream) * ((double)STREAM_ARRAY_SIZE / 1024.0 / 1024.0) << " MiB (= " <<
	(3.0 * SizeOfStream) * ((double)STREAM_ARRAY_SIZE / 1024.0 / 1024.0 / 1024.0) << " GiB). " << std::endl;

    std::cout << "Each kernel will be executed " << NTIMES << " times." << std::endl;
}

template <typename T>
static double Copy(std::vector<T>& c, const std::vector<T>& a)
{
    const auto SIZE = c.size();
    auto start = std::chrono::system_clock::now();
#pragma omp parallel for
    for(auto j = decltype(SIZE)(0); j < SIZE; ++j) {
	c[j] = a[j];
    }
    auto end = std::chrono::system_clock::now();
    std::chrono::duration<double, std::milli> dur = (end - start);

    return dur.count() / 1000.0;
}

template <typename T>
static double Scale(std::vector<T>& b, const std::vector<T>& c, T scalar)
{
    const auto SIZE = b.size();
    auto start = std::chrono::system_clock::now();
#pragma omp parallel for
    for(auto j = decltype(SIZE)(0); j < SIZE; ++j) {
	b[j] = scalar * c[j];
    }
    auto end = std::chrono::system_clock::now();
    std::chrono::duration<double, std::milli> dur = (end - start);

    return dur.count() / 1000.0;
}

template <typename T>
static double Add(std::vector<T>& c, const std::vector<T>& a, const std::vector<T>& b)
{
    const auto SIZE = c.size();
    auto start = std::chrono::system_clock::now();
#pragma omp parallel for
    for(auto j = decltype(SIZE)(0); j < SIZE; ++j) {
	c[j] = a[j] + b[j];
    }
    auto end = std::chrono::system_clock::now();
    std::chrono::duration<double, std::milli> dur = (end - start);

    return dur.count() / 1000.0;
}

template <typename T>
static double Triad(std::vector<T>& a, const std::vector<T>& b, const std::vector<T>& c, T scalar)
{
    const auto SIZE = a.size();
    auto start = std::chrono::system_clock::now();
#pragma omp parallel for
    for(auto j = decltype(SIZE)(0); j < SIZE; ++j) {
	a[j] = b[j] + scalar * c[j];
    }
    auto end = std::chrono::system_clock::now();
    std::chrono::duration<double, std::milli> dur = (end - start);

    return dur.count() / 1000.0;
}

template <typename STREAM_TYPE>
static void ShowSummary(const std::vector<std::vector<double>>& times, size_t STREAM_ARRAY_SIZE, size_t NTIMES)
{
    // Summary
    std::vector<double> avgtime(4, 0);
    std::vector<double> mintime(4, std::numeric_limits<double>::max());
    std::vector<double> maxtime(4, std::numeric_limits<double>::min());

    for(auto k = decltype(NTIMES)(1); k < NTIMES; ++k) {
	const auto& cur_time = times[k];

	for(size_t j = 0; j < 4; ++j) {
	    avgtime[j] += cur_time[j];
	    mintime[j] = std::min(mintime[j], cur_time[j]);
	    maxtime[j] = std::max(maxtime[j], cur_time[j]);
	}
    }

    const std::string label[4] = {"Copy:      ", "Scale:     ",
				  "Add:	   ", "Triad:     "};

    const double bytes[4] = {    2.0 * sizeof(STREAM_TYPE) * STREAM_ARRAY_SIZE,
				 2.0 * sizeof(STREAM_TYPE) * STREAM_ARRAY_SIZE,
				 3.0 * sizeof(STREAM_TYPE) * STREAM_ARRAY_SIZE,
				 3.0 * sizeof(STREAM_TYPE) * STREAM_ARRAY_SIZE};

    // 諦めた
    printf("Function    Best Rate MB/s  Avg time     Min time     Max time\n");
    for(size_t j = 0; j < 4; ++j) {
	auto avg = avgtime[j] / static_cast<double>(NTIMES - 1);

	printf("%s%12.1f  %11.6f  %11.6f  %11.6f\n", label[j].c_str(),
	       1.0e-6 * bytes[j] / mintime[j],
	       avg, mintime[j], maxtime[j]);
    }
    const std::string HLINE = "-------------------------------------------------------------";
    std::cout << HLINE << std::endl;
}

template <typename STREAM_TYPE>
void checkSTREAMresults (const std::vector<STREAM_TYPE>& a, const std::vector<STREAM_TYPE>& b, const std::vector<STREAM_TYPE>& c,
			 size_t STREAM_ARRAY_SIZE, size_t NTIMES)
{
    STREAM_TYPE aj,bj,cj,scalar;
    STREAM_TYPE aSumErr,bSumErr,cSumErr;
    STREAM_TYPE aAvgErr,bAvgErr,cAvgErr;
    double epsilon;
    ssize_t	j;
    int	k,ierr,err;

    /* reproduce initialization */
    aj = 1.0;
    bj = 2.0;
    cj = 0.0;
    /* a[] is modified during timing check */
    aj = 2.0E0 * aj;
    /* now execute timing loop */
    scalar = 3.0;
    for (k=0; k<NTIMES; k++)
    {
	cj = aj;
	bj = scalar*cj;
	cj = aj+bj;
	aj = bj+scalar*cj;
    }

    /* accumulate deltas between observed and expected results */
    aSumErr = 0.0;
    bSumErr = 0.0;
    cSumErr = 0.0;
    for (j=0; j<STREAM_ARRAY_SIZE; j++) {
	aSumErr += abs(a[j] - aj);
	bSumErr += abs(b[j] - bj);
	cSumErr += abs(c[j] - cj);
	// if (j == 417) printf("Index 417: c[j]: %f, cj: %f\n",c[j],cj);	// MCCALPIN
    }
    aAvgErr = aSumErr / (STREAM_TYPE) STREAM_ARRAY_SIZE;
    bAvgErr = bSumErr / (STREAM_TYPE) STREAM_ARRAY_SIZE;
    cAvgErr = cSumErr / (STREAM_TYPE) STREAM_ARRAY_SIZE;

    if (sizeof(STREAM_TYPE) == 4) {
	epsilon = 1.e-6;
    }
    else if (sizeof(STREAM_TYPE) == 8) {
	epsilon = 1.e-13;
    }
    else {
	printf("WEIRD: sizeof(STREAM_TYPE) = %lu\n",sizeof(STREAM_TYPE));
	epsilon = 1.e-6;
    }

    err = 0;
    if (abs(aAvgErr/aj) > epsilon) {
	err++;
	printf ("Failed Validation on array a[], AvgRelAbsErr > epsilon (%e)\n",epsilon);
	printf ("     Expected Value: %e, AvgAbsErr: %e, AvgRelAbsErr: %e\n",aj,aAvgErr,abs(aAvgErr)/aj);
	ierr = 0;
	for (j=0; j<STREAM_ARRAY_SIZE; j++) {
	    if (abs(a[j]/aj-1.0) > epsilon) {
		ierr++;
#ifdef VERBOSE
		if (ierr < 10) {
		    printf("         array a: index: %ld, expected: %e, observed: %e, relative error: %e\n",
			   j,aj,a[j],abs((aj-a[j])/aAvgErr));
		}
#endif
	    }
	}
	printf("     For array a[], %d errors were found.\n",ierr);
    }
    if (abs(bAvgErr/bj) > epsilon) {
	err++;
	printf ("Failed Validation on array b[], AvgRelAbsErr > epsilon (%e)\n",epsilon);
	printf ("     Expected Value: %e, AvgAbsErr: %e, AvgRelAbsErr: %e\n",bj,bAvgErr,abs(bAvgErr)/bj);
	printf ("     AvgRelAbsErr > Epsilon (%e)\n",epsilon);
	ierr = 0;
	for (j=0; j<STREAM_ARRAY_SIZE; j++) {
	    if (abs(b[j]/bj-1.0) > epsilon) {
		ierr++;
#ifdef VERBOSE
		if (ierr < 10) {
		    printf("         array b: index: %ld, expected: %e, observed: %e, relative error: %e\n",
			   j,bj,b[j],abs((bj-b[j])/bAvgErr));
		}
#endif
	    }
	}
	printf("     For array b[], %d errors were found.\n",ierr);
    }
    if (abs(cAvgErr/cj) > epsilon) {
	err++;
	printf ("Failed Validation on array c[], AvgRelAbsErr > epsilon (%e)\n",epsilon);
	printf ("     Expected Value: %e, AvgAbsErr: %e, AvgRelAbsErr: %e\n",cj,cAvgErr,abs(cAvgErr)/cj);
	printf ("     AvgRelAbsErr > Epsilon (%e)\n",epsilon);
	ierr = 0;
	for (j=0; j<STREAM_ARRAY_SIZE; j++) {
	    if (abs(c[j]/cj-1.0) > epsilon) {
		ierr++;
#ifdef VERBOSE
		if (ierr < 10) {
		    printf("         array c: index: %ld, expected: %e, observed: %e, relative error: %e\n",
			   j,cj,c[j],abs((cj-c[j])/cAvgErr));
		}
#endif
	    }
	}
	printf("     For array c[], %d errors were found.\n",ierr);
    }
    if (err == 0) {
	printf ("Solution Validates: avg error less than %e on all three arrays\n",epsilon);
    }
#ifdef VERBOSE
    printf ("Results Validation Verbose Results: \n");
    printf ("    Expected a(1), b(1), c(1): %f %f %f \n",aj,bj,cj);
    printf ("    Observed a(1), b(1), c(1): %f %f %f \n",a[1],b[1],c[1]);
    printf ("    Rel Errors on a, b, c:     %e %e %e \n",abs(aAvgErr/aj),abs(bAvgErr/bj),abs(cAvgErr/cj));
#endif
}

int main(int argc, char** argv)
{
    ParseArg parser(argc, argv);

    const size_t STREAM_ARRAY_SIZE = parser.getArraySize();
    const size_t NTIMES = parser.getTimes();

    constexpr size_t OFFSET = 0;

    using STREAM_TYPE = double;

    auto a = std::vector<STREAM_TYPE>(STREAM_ARRAY_SIZE + OFFSET);
    auto b = std::vector<STREAM_TYPE>(STREAM_ARRAY_SIZE + OFFSET);
    auto c = std::vector<STREAM_TYPE>(STREAM_ARRAY_SIZE + OFFSET);

    PrintMessages(STREAM_ARRAY_SIZE, NTIMES, OFFSET, sizeof(STREAM_TYPE));

    // initial
    std::fill(a.begin(), a.end(), 1.0);
    std::fill(b.begin(), b.end(), 2.0);
    std::fill(c.begin(), c.end(), 0.0);

    for(auto& ia : a) {
	ia *= 2.0E0;
    }

    std::vector<std::vector<double>> times(NTIMES);
    for(auto& it : times) {
	it.resize(4);
    }

    STREAM_TYPE scalar = 3.0;
    for(auto k = decltype(NTIMES)(0); k < NTIMES; ++k) {
	// copy
	times[k][0] = Copy(c, a);

	// scale
	times[k][1] = Scale(b, c, scalar);

	// add
	times[k][2] = Add(c, a, b);

	// triad
	times[k][3] = Triad(a, b, c, scalar);
    }

    ShowSummary<STREAM_TYPE>(times, STREAM_ARRAY_SIZE, NTIMES);
    checkSTREAMresults<STREAM_TYPE>(a, b, c, STREAM_ARRAY_SIZE, NTIMES);

    return 0;
}
