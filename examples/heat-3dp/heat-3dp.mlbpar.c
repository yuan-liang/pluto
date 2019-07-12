#include <omp.h>
#include <math.h>
#define ceild(n,d)  ceil(((double)(n))/((double)(d)))
#define floord(n,d) floor(((double)(n))/((double)(d)))
#define max(x,y)    ((x) > (y)? (x) : (y))
#define min(x,y)    ((x) < (y)? (x) : (y))

/*
 * Discretized 3D heat equation stencil with non periodic boundary conditions
 * Adapted from Pochoir test bench
 *
 * Irshad Pananilath: irshad@csa.iisc.ernet.in
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>

/*
 * N is the number of points
 * T is the number of timesteps
 */
#ifdef HAS_DECLS
#include "decls.h"
#else
#define N 300L
#define T 500L
#endif

#define NUM_FP_OPS 15

/* Define our arrays */

double A[2][N][N][N];
double total;
double sum_err_sqr = 0;
int chtotal = 0;

/* Subtract the `struct timeval' values X and Y,
 * storing the result in RESULT.
 *
 * Return 1 if the difference is negative, otherwise 0.
 */
int timeval_subtract(struct timeval *result, struct timeval *x,
                     struct timeval *y) {
  /* Perform the carry for the later subtraction by updating y. */
  if (x->tv_usec < y->tv_usec) {
    int nsec = (y->tv_usec - x->tv_usec) / 1000000 + 1;

    y->tv_usec -= 1000000 * nsec;
    y->tv_sec += nsec;
  }

  if (x->tv_usec - y->tv_usec > 1000000) {
    int nsec = (x->tv_usec - y->tv_usec) / 1000000;

    y->tv_usec += 1000000 * nsec;
    y->tv_sec -= nsec;
  }

  /* Compute the time remaining to wait.
   * tv_usec is certainly positive.
   */
  result->tv_sec = x->tv_sec - y->tv_sec;
  result->tv_usec = x->tv_usec - y->tv_usec;

  /* Return 1 if result is negative. */
  return x->tv_sec < y->tv_sec;
}

int main(int argc, char *argv[]) {
  long int t, i, j, k;
  const int BASE = 1024;
  long count = 0;
  // for timekeeping
  int ts_return = -1;
  struct timeval start, end, result;
  double tdiff = 0.0;

  printf("Number of points = %ld\t|Number of timesteps = %ld\t", N, T);

  /* Initialization */
  srand(42); // seed with a constant value to verify results

  for (i = 0; i < N; i++) {
    for (j = 0; j < N; j++) {
      for (k = 0; k < N; k++) {
        A[0][i][j][k] = 1.0 * (rand() % BASE);
      }
    }
  }

#ifdef TIME
  gettimeofday(&start, 0);
#endif

  short _N = N - 1;
  short _T = T;

  int t1, t2, t3, t4, t5, t6, t7;
 int lb, ub, lbp, ubp, lb2, ub2;
 register int lbv, ubv;
if (_T >= 1) {
  for (t1=-11;t1<=floord(_T-1,8);t1++) {
    lbp=max(0,ceild(t1,2));
    ubp=min(floord(t1+28,2),floord(_T+148,16));
#pragma omp parallel for private(lbv,ubv,t3,t4,t5,t6,t7)
    for (t2=lbp;t2<=ubp;t2++) {
      for (t3=max(max(max(0,ceild(t1-1,2)),-t1-2),-t1+3*t2-29);t3<=min(min(min(min(floord(t1+29,2),floord(-8*t1+8*t2+_T-1,8)),floord(_T+148,16)),floord(-16*t1+3*_T-3,16)),-t1+3*t2+1);t3++) {
        for (t4=max(max(max(max(max(0,ceild(16*t1+16*t3,3)),8*t1),16*t2-149),16*t3-149),8*t1-8*t2+8*t3);t4<=min(min(min(min(min(floord(16*t1+16*t3+45,3),_T-1),8*t1+89),16*t2+15),16*t3+15),8*t1-8*t2+8*t3+89);t4++) {
          for (t5=max(max(max(16*t2,t4),-16*t1+16*t2+2*t4-164),-16*t1+16*t2-16*t3+3*t4-30);t5<=min(min(min(16*t2+15,t4+149),-16*t1+16*t2+2*t4),-16*t1+16*t2-16*t3+3*t4);t5++) {
            for (t6=max(max(16*t3,t4),-16*t1+16*t2+3*t4-t5-15);t6<=min(min(16*t3+15,t4+149),-16*t1+16*t2+3*t4-t5);t6++) {
              lbv=0;
              ubv=299;
#pragma ivdep
#pragma vector always
              for (t7=lbv;t7<=ubv;t7++) {
                A[(t4 + 1) % 2][(-t4+t5)][(-t4+t6)][t7] = ((((0.125 * ((A[t4 % 2][(-t4+t5) == _N ? 0 : (-t4+t5) + 1][(-t4+t6)][t7] - (2.0 * A[t4 % 2][(-t4+t5)][(-t4+t6)][t7])) + A[t4 % 2][(-t4+t5) == 0 ? _N : (-t4+t5) - 1][(-t4+t6)][t7])) + (0.125 * ((A[t4 % 2][(-t4+t5)][(-t4+t6) == _N ? 0 : (-t4+t6) + 1][t7] - (2.0 * A[t4 % 2][(-t4+t5)][(-t4+t6)][t7])) + A[t4 % 2][(-t4+t5)][(-t4+t6) == 0 ? _N : (-t4+t6) - 1][t7]))) + (0.125 * ((A[t4 % 2][(-t4+t5)][(-t4+t6)][t7 == 0 ? _N : t7 - 1] - (2.0 * A[t4 % 2][(-t4+t5)][(-t4+t6)][t7])) + A[t4 % 2][(-t4+t5)][(-t4+t6)][t7 == _N ? 0 : t7 + 1]))) + A[t4 % 2][(-t4+t5)][(-t4+t6)][t7]);;
                A[(t4 + 1) % 2][(t4-t5+299)][(-t4+t6)][t7] = ((((0.125 * ((A[t4 % 2][(t4-t5+299) == _N ? 0 : (t4-t5+299) + 1][(-t4+t6)][t7] - (2.0 * A[t4 % 2][(t4-t5+299)][(-t4+t6)][t7])) + A[t4 % 2][(t4-t5+299) == 0 ? _N : (t4-t5+299) - 1][(-t4+t6)][t7])) + (0.125 * ((A[t4 % 2][(t4-t5+299)][(-t4+t6) == _N ? 0 : (-t4+t6) + 1][t7] - (2.0 * A[t4 % 2][(t4-t5+299)][(-t4+t6)][t7])) + A[t4 % 2][(t4-t5+299)][(-t4+t6) == 0 ? _N : (-t4+t6) - 1][t7]))) + (0.125 * ((A[t4 % 2][(t4-t5+299)][(-t4+t6)][t7 == 0 ? _N : t7 - 1] - (2.0 * A[t4 % 2][(t4-t5+299)][(-t4+t6)][t7])) + A[t4 % 2][(t4-t5+299)][(-t4+t6)][t7 == _N ? 0 : t7 + 1]))) + A[t4 % 2][(t4-t5+299)][(-t4+t6)][t7]);;
                A[(t4 + 1) % 2][(-t4+t5)][(t4-t6+299)][t7] = ((((0.125 * ((A[t4 % 2][(-t4+t5) == _N ? 0 : (-t4+t5) + 1][(t4-t6+299)][t7] - (2.0 * A[t4 % 2][(-t4+t5)][(t4-t6+299)][t7])) + A[t4 % 2][(-t4+t5) == 0 ? _N : (-t4+t5) - 1][(t4-t6+299)][t7])) + (0.125 * ((A[t4 % 2][(-t4+t5)][(t4-t6+299) == _N ? 0 : (t4-t6+299) + 1][t7] - (2.0 * A[t4 % 2][(-t4+t5)][(t4-t6+299)][t7])) + A[t4 % 2][(-t4+t5)][(t4-t6+299) == 0 ? _N : (t4-t6+299) - 1][t7]))) + (0.125 * ((A[t4 % 2][(-t4+t5)][(t4-t6+299)][t7 == 0 ? _N : t7 - 1] - (2.0 * A[t4 % 2][(-t4+t5)][(t4-t6+299)][t7])) + A[t4 % 2][(-t4+t5)][(t4-t6+299)][t7 == _N ? 0 : t7 + 1]))) + A[t4 % 2][(-t4+t5)][(t4-t6+299)][t7]);;
                A[(t4 + 1) % 2][(t4-t5+299)][(t4-t6+299)][t7] = ((((0.125 * ((A[t4 % 2][(t4-t5+299) == _N ? 0 : (t4-t5+299) + 1][(t4-t6+299)][t7] - (2.0 * A[t4 % 2][(t4-t5+299)][(t4-t6+299)][t7])) + A[t4 % 2][(t4-t5+299) == 0 ? _N : (t4-t5+299) - 1][(t4-t6+299)][t7])) + (0.125 * ((A[t4 % 2][(t4-t5+299)][(t4-t6+299) == _N ? 0 : (t4-t6+299) + 1][t7] - (2.0 * A[t4 % 2][(t4-t5+299)][(t4-t6+299)][t7])) + A[t4 % 2][(t4-t5+299)][(t4-t6+299) == 0 ? _N : (t4-t6+299) - 1][t7]))) + (0.125 * ((A[t4 % 2][(t4-t5+299)][(t4-t6+299)][t7 == 0 ? _N : t7 - 1] - (2.0 * A[t4 % 2][(t4-t5+299)][(t4-t6+299)][t7])) + A[t4 % 2][(t4-t5+299)][(t4-t6+299)][t7 == _N ? 0 : t7 + 1]))) + A[t4 % 2][(t4-t5+299)][(t4-t6+299)][t7]);;
              }
            }
          }
        }
      }
    }
  }
}
#ifdef TIME
  gettimeofday(&end, 0);

  ts_return = timeval_subtract(&result, &end, &start);
  tdiff = (double)(result.tv_sec + result.tv_usec * 1.0e-6);

  printf("|Time taken: %7.5lfms\t", tdiff * 1.0e3);
  printf("|MFLOPS: %f\n",
         ((((double)NUM_FP_OPS * N * N * N * (T - 1)) / tdiff) / 1000000L));
#endif

  if (fopen(".test", "r")) {
    total = 0;
    for (i = 0; i < N; i++) {
      for (j = 0; j < N; j++) {
        for (k = 0; k < N; k++) {
          total += A[T % 2][i][j][k];
        }
      }
    }
    fprintf(stderr, "|sum: %e\t", total);
    for (i = 0; i < N; i++) {
      for (j = 0; j < N; j++) {
        for (k = 0; k < N; k++) {
          sum_err_sqr += (A[T % 2][i][j][k] - (total / N)) *
                         (A[T % 2][i][j][k] - (total / N));
        }
      }
    }
    fprintf(stderr, "|rms(A) = %7.2f\t", sqrt(sum_err_sqr));
    for (i = 0; i < N; i++) {
      for (j = 0; j < N; j++) {
        for (k = 0; k < N; k++) {
          chtotal += ((char *)A[T % 2][i][j])[k];
        }
      }
    }
    fprintf(stderr, "|sum(rep(A)) = %d\n", chtotal);
  }

  return 0;
}

// icc -O3 -fp-model precise heat_1d_np.c -o op-heat-1d-np -lm
// /* @ begin PrimeTile (num_tiling_levels=1; first_depth=1; last_depth=-1;
// boundary_tiling_level=-1;) @*/
// /* @ begin PrimeRegTile (scalar_replacement=0; T1t5=4; T1t6=4; T1t7=4;
// T1t8=4; ) @*/
// /* @ end @*/
