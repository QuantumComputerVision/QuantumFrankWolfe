#include "mex.h"
#include "matrix.h"
#include<math.h>
#include <Eigen/Core>
#include <Eigen/Dense>
#include<vector>
#include<bitset>
#include<iostream>
#include <stdlib.h>
//#include<omp.h>


double solve_exhaustive_binary(double* Q, const mwSize n, double* q_result)
{
	mwSize N = (mwSize)pow(2, n);

	// Don't do the following in your C++ MEX main file
	Eigen::MatrixXd Q_eig = Eigen::Map<Eigen::MatrixXd>(Q, n, n);
	//Eigen::VectorXd q_result = Eigen::Map<Eigen::VectorXd>(qResult, n);
	Eigen::MatrixXd q(n, 1);

	double minErr = 99999999;

// #pragma omp parallel for
	for (mwSize i = 0; i < N; i++)
	{
		size_t j;
		auto qbits = std::bitset<512>(i);
		
		// qbits; 
		for (j = 0; j < n; j++)
			q(n-j-1,0) = qbits[j];

		double err = (q.transpose() * Q_eig * q)(0, 0);
		// #pragma omp critical {
		if (err < minErr) {
			for (j = 0; j < n; j++)
				q_result[j] = q(j);
			minErr = err;
		}
		// 		}
	}

	return minErr;
}

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
	// prhs[0]: a cell array of length T, each element is a vector with different lengths
	// mexPrintf("getting dims...");

	const mwSize* dims = mxGetDimensions(prhs[0]);
	const mwSize M = dims[0];
	const mwSize N = dims[1];

	if (M != N)
		mexPrintf("Matrix not square: %d x %d", M, N);

	double* Qptr = mxGetPr(prhs[0]);
	// Eigen::Map<const Eigen::MatrixXd> Q(mxGetPr(plhs[0]), lhsSize);

	// create a cell matrix with T rows and one columns
	// mexPrintf("mxCreateDoubleMatrix...");
	plhs[0] = mxCreateDoubleMatrix(M, 1, mxREAL);
	double* qoutPtr = mxGetPr(plhs[0]);

	// mexPrintf("running solve_exhaustive_binary...");
	solve_exhaustive_binary(Qptr, M, qoutPtr);
	
}
/*
void main()
{
	double sampleQ[81] = { -1.210090, 0.452396, -0.160053, 0.068573, 0.777775, -0.008967, 0.274021, 0.228899, -0.352809, 0.452396, 0.188618, -0.206941, 0.848165, -0.237047, -0.072873, -0.032667, 0.045315, 0.204572, -0.160053, -0.206941, -0.417389, 0.157309, -0.772926, 0.852817, 0.589696, -0.003636, -0.504836, 0.068573, 0.848165, 0.157309, -0.277294, 0.696767, -0.279095, 0.404297, 0.954596, -0.351223, 0.777775, -0.237047, -0.772926, 0.696767, -0.487711, -0.132707, 0.209891, 0.536868, -0.085013, -0.008967, -0.072873, 0.852817, -0.279095, -0.132707, 0.104150, -0.306324, -0.732855, -0.439518, 0.274021, -0.032667, 0.589696, 0.404297, 0.209891, -0.306324, 1.454526, -0.810736, -0.494784, 0.228899, 0.045315, -0.003636, 0.954596, 0.536868, -0.732855, -0.810736, -0.285088, 0.444851, -0.352809, 0.204572, -0.504836, -0.351223, -0.085013, -0.439518, -0.494784, 0.444851, -1.022690 };
	mwSize sampleN = 9;
	double qResult[9];

	double minErr = solve_exhaustive_binary(sampleQ, sampleN, qResult);
	for (int i = 0; i < 9; i++)
		mexPrintf("%f, \n", qResult[i]);
	mexPrintf("\n");
	mexPrintf("success %g\n", minErr);
}*/