#include <math.h>
#define ceild(n,d)  ceil(((double)(n))/((double)(d)))
#define floord(n,d) floor(((double)(n))/((double)(d)))
#define max(x,y)    ((x) > (y)? (x) : (y))
#define min(x,y)    ((x) < (y)? (x) : (y))

/*
 * Discretized 2D heat equation stencil with non periodic boundary conditions
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
#define N 1600L
#define T 500L
#endif

#define NUM_FP_OPS 10

/* Define our arrays */
double A[2][N][N];
double sum_err_sqr = 0;
int chtotal = 0;

int timeval_subtract(struct timeval *result, struct timeval *x,
                     struct timeval *y) {
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

  result->tv_sec = x->tv_sec - y->tv_sec;
  result->tv_usec = x->tv_usec - y->tv_usec;

  return x->tv_sec < y->tv_sec;
}

int main(int argc, char *argv[]) {
  long int t, i, j, k;
  const int BASE = 1024;

  // for timekeeping
  int ts_return = -1;
  struct timeval start, end, result;
  double tdiff = 0.0;
  double total;

  printf("Number of points = %ld\t|Number of timesteps = %ld\t", N * N, T);

  /* Initialization */
  srand(42); // seed with a constant value to verify results

  for (i = 0; i < N; i++) {
    for (j = 0; j < N; j++) {
      A[0][i][j] = 1.0 * (rand() % BASE);
    }
  }

#ifdef TIME
  gettimeofday(&start, 0);
#endif

  short _N = N - 1;
  int _T = T;

  int t1, t2, t3, t4, t5, t6;
 register int lbv, ubv;
if (_T >= 1) {
  for (t1=-13;t1<=floord(_T-1,64);t1++) {
    for (t2=max(t1,-t1-1);t2<=min(min(floord(-32*t1+_T-1,32),floord(_T+798,64)),t1+25);t2++) {
      for (t3=max(0,ceild(t1+t2-1,2));t3<=min(floord(t1+t2+26,2),floord(_T+798,64));t3++) {
        for (t4=max(max(max(0,32*t1+32*t2),64*t2-799),64*t3-799);t4<=min(min(min(_T-1,64*t1+862),64*t3+63),32*t1+32*t2+63);t4++) {
          for (t5=max(max(64*t2,t4),-64*t1+2*t4-63);t5<=min(min(-64*t1+2*t4,64*t2+63),t4+799);t5++) {
            lbv=max(64*t3,t4);
            ubv=min(64*t3+63,t4+799);
#pragma ivdep
#pragma vector always
            for (t6=lbv;t6<=ubv;t6++) {
              A[(t4 + 1) % 2][(-t4+t5)][(-t4+t6)] = (((0.125 * ((A[t4 % 2][(-t4+t5) == _N ? 0 : (-t4+t5) + 1][(-t4+t6)] - (2.0 * A[t4 % 2][(-t4+t5)][(-t4+t6)])) + A[t4 % 2][(-t4+t5) == 0 ? _N : (-t4+t5) - 1][(-t4+t6)])) + (0.125 * ((A[t4 % 2][(-t4+t5)][(-t4+t6) == _N ? 0 : (-t4+t6) + 1] - (2.0 * A[t4 % 2][(-t4+t5)][(-t4+t6)])) + A[t4 % 2][(-t4+t5)][(-t4+t6) == 0 ? _N : (-t4+t6) - 1]))) + A[t4 % 2][(-t4+t5)][(-t4+t6)]);;
              A[(t4 + 1) % 2][(t4-t5+1599)][(-t4+t6)] = (((0.125 * ((A[t4 % 2][(t4-t5+1599) == _N ? 0 : (t4-t5+1599) + 1][(-t4+t6)] - (2.0 * A[t4 % 2][(t4-t5+1599)][(-t4+t6)])) + A[t4 % 2][(t4-t5+1599) == 0 ? _N : (t4-t5+1599) - 1][(-t4+t6)])) + (0.125 * ((A[t4 % 2][(t4-t5+1599)][(-t4+t6) == _N ? 0 : (-t4+t6) + 1] - (2.0 * A[t4 % 2][(t4-t5+1599)][(-t4+t6)])) + A[t4 % 2][(t4-t5+1599)][(-t4+t6) == 0 ? _N : (-t4+t6) - 1]))) + A[t4 % 2][(t4-t5+1599)][(-t4+t6)]);;
              A[(t4 + 1) % 2][(-t4+t5)][(t4-t6+1599)] = (((0.125 * ((A[t4 % 2][(-t4+t5) == _N ? 0 : (-t4+t5) + 1][(t4-t6+1599)] - (2.0 * A[t4 % 2][(-t4+t5)][(t4-t6+1599)])) + A[t4 % 2][(-t4+t5) == 0 ? _N : (-t4+t5) - 1][(t4-t6+1599)])) + (0.125 * ((A[t4 % 2][(-t4+t5)][(t4-t6+1599) == _N ? 0 : (t4-t6+1599) + 1] - (2.0 * A[t4 % 2][(-t4+t5)][(t4-t6+1599)])) + A[t4 % 2][(-t4+t5)][(t4-t6+1599) == 0 ? _N : (t4-t6+1599) - 1]))) + A[t4 % 2][(-t4+t5)][(t4-t6+1599)]);;
              A[(t4 + 1) % 2][(t4-t5+1599)][(t4-t6+1599)] = (((0.125 * ((A[t4 % 2][(t4-t5+1599) == _N ? 0 : (t4-t5+1599) + 1][(t4-t6+1599)] - (2.0 * A[t4 % 2][(t4-t5+1599)][(t4-t6+1599)])) + A[t4 % 2][(t4-t5+1599) == 0 ? _N : (t4-t5+1599) - 1][(t4-t6+1599)])) + (0.125 * ((A[t4 % 2][(t4-t5+1599)][(t4-t6+1599) == _N ? 0 : (t4-t6+1599) + 1] - (2.0 * A[t4 % 2][(t4-t5+1599)][(t4-t6+1599)])) + A[t4 % 2][(t4-t5+1599)][(t4-t6+1599) == 0 ? _N : (t4-t6+1599) - 1]))) + A[t4 % 2][(t4-t5+1599)][(t4-t6+1599)]);;
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

  printf("|Time taken =  %7.5lfms\t", tdiff * 1.0e3);
  printf("|MFLOPS =  %f\n",
         ((((double)NUM_FP_OPS * N * N * T) / tdiff) / 1000000L));
#endif

  if (fopen(".test", "r")) {
    total = 0;
    for (i = 0; i < N; i++) {
      for (j = 0; j < N; j++) {
        total += A[T % 2][i][j];
      }
    }
    fprintf(stderr, "|sum: %e\t", total);
    for (i = 0; i < N; i++) {
      for (j = 0; j < N; j++) {
        sum_err_sqr +=
            (A[T % 2][i][j] - (total / N)) * (A[T % 2][i][j] - (total / N));
      }
    }
    fprintf(stderr, "|rms(A) = %7.2f\t", sqrt(sum_err_sqr));
    for (i = 0; i < N; i++) {
      for (j = 0; j < N; j++) {
        chtotal += ((char *)A[T % 2][i])[j];
      }
    }
    fprintf(stderr, "|sum(rep(A)) = %d\n", chtotal);
  }
  return 0;
}

// icc -O3 -fp-model precise heat_1d_np.c -o op-heat-1d-np -lm
// /* @ begin PrimeTile (num_tiling_levels=1; first_depth=1; last_depth=-1;
// boundary_tiling_level=-1;) @*/
// /* @ begin PrimeRegTile (scalar_replacement=0; T1t3=8; T1t4=8; ) @*/
// /* @ end @*/
