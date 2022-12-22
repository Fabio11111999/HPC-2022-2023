#include "walltime.h"
#include <math.h>
#include <stdlib.h>
#include <omp.h>

double fastexp(double base, int esp) {
	double res = 1.0;
	while (esp) {
		if (esp % 2) {
			res *= base;
		}
		esp /= 2;
		base *= base;
	}
	return res;
}
int main(int argc, char *argv[]) {
  long long N = 2000000000;
  double up = 1.00000001;
  double Sn = 1.00000001;
  int n;
  /* allocate memory for the recursion */
  double *opt = (double *)malloc((N + 1) * sizeof(double));

  if (opt == NULL)
    die("failed to allocate problem size");

  double time_start = wall_time();
  // TODO: YOU NEED TO PARALLELIZE THIS LOOP

#pragma omp	parallel 
{
	int i, id, nths, istart, iend;
	id = omp_get_thread_num();
	nths = omp_get_num_threads();
	istart = id * N / nths;
	iend = (id + 1) * N / nths;
	double factor = 1.0;
	if (id == nths - 1) {
		iend = N;
	}
	opt[istart] = Sn * fastexp(up, istart);
	factor *= up;
	for (i = istart + 1; i < iend; i++) {
		opt[i] = opt[i - 1] * up;
		factor *= up;
	}
#pragma omp critical
	Sn *= factor;
}

  printf("Parallel RunTime   :  %f seconds\n", wall_time() - time_start);
  printf("Final Result Sn    :  %.17g \n", Sn);

  double temp = 0.0;
  for (n = 0; n <= N; ++n) {
    temp += opt[n] * opt[n];
  }
  printf("Result ||opt||^2_2 :  %f\n", temp / (double)N);
  printf("\n");

  return 0;
}
