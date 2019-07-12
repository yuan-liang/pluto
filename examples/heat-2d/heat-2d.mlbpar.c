#include <omp.h>

#define S1(zT3,zT4,zT5,t,i,j)	A[(t + 1) % 2][i][j] = (((0.125 * ((A[t % 2][i + 1][j] - (2.0 * A[t % 2][i][j])) + A[t % 2][i - 1][j])) + (0.125 * ((A[t % 2][i][j + 1] - (2.0 * A[t % 2][i][j])) + A[t % 2][i][j - 1]))) + A[t % 2][i][j]);

		int t1, t2, t3, t4, t5, t6;

	int lb, ub, lbp, ubp, lb2, ub2;
	register int lbv, ubv;

/* Start of CLooG code */
for (t1=-64;t1<=31;t1++) {
  lbp=max(0,ceild(t1,2));
  ubp=min(min(78,floord(t1+188,2)),t1+125);
#pragma omp parallel for private(lbv,ubv,t3,t4,t5,t6)
  for (t2=lbp;t2<=ubp;t2++) {
    for (t3=max(max(max(0,ceild(t1-1,2)),-t1-2),-t1+3*t2-189);t3<=min(min(floord(t1+189,2),-t1+46),-t1+3*t2+1);t3++) {
      for (t4=max(max(max(max(max(max(0,ceild(64*t1+64*t3,3)),32*t1+1),64*t2-4000),64*t3-4000),64*t1-64*t2+2),32*t1-32*t2+32*t3+1);t4<=min(min(min(min(min(min(999,floord(64*t1+64*t3+189,3)),32*t1+2063),64*t2+62),64*t3+62),64*t1-64*t2+8063),32*t1-32*t2+32*t3+2063);t4++) {
        for (t5=max(max(max(64*t2,t4+1),-64*t1+64*t2+2*t4-4063),-64*t1+64*t2-64*t3+3*t4-126);t5<=min(min(min(64*t2+63,t4+4000),-64*t1+64*t2-64*t3+3*t4),-64*t1+64*t2+2*t4-1);t5++) {
          lbv=max(max(64*t3,t4+1),-64*t1+64*t2+3*t4-t5-63);
          ubv=min(min(64*t3+63,t4+4000),-64*t1+64*t2+3*t4-t5);
          #pragma ivdep
          #pragma vector always
          for (t6=lbv;t6<=ubv;t6++) {
            S1((t1-t2),t2,t3,t4,(-t4+t5),(-t4+t6));
          }
        }
      }
    }
  }
}
/* End of CLooG code */
