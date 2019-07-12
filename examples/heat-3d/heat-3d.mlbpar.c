#include <omp.h>

#define S1(zT4,zT5,zT6,zT7,t,i,j,k)	A[(t + 1) % 2][i][j][k] = ((((0.125 * ((A[t % 2][i + 1][j][k] - (2.0 * A[t % 2][i][j][k])) + A[t % 2][i - 1][j][k])) + (0.125 * ((A[t % 2][i][j + 1][k] - (2.0 * A[t % 2][i][j][k])) + A[t % 2][i][j - 1][k]))) + (0.125 * ((A[t % 2][i][j][k - 1] - (2.0 * A[t % 2][i][j][k])) + A[t % 2][i][j][k + 1]))) + A[t % 2][i][j][k]);

		int t1, t2, t3, t4, t5, t6, t7, t8;

	int lb, ub, lbp, ubp, lb2, ub2;
	register int lbv, ubv;

/* Start of CLooG code */
for (t1=-39;t1<=62;t1++) {
  lbp=max(max(0,ceild(t1,2)),t1-30);
  ubp=min(49,floord(t1+75,2));
#pragma omp parallel for private(lbv,ubv,t3,t4,t5,t6,t7,t8)
  for (t2=lbp;t2<=ubp;t2++) {
    for (t3=max(max(max(max(0,ceild(t1-1,2)),-t1-21),t2-19),-t1+3*t2-76);t3<=min(min(min(min(49,floord(t1+76,2)),-t1+93),t2+19),-t1+3*t2+1);t3++) {
      for (t5=max(max(max(max(max(max(0,ceild(16*t1+16*t3+1,3)),8*t1+1),16*t2-300),16*t3-300),16*t1-16*t2+3),8*t1-8*t2+8*t3+1);t5<=min(min(min(min(min(498,floord(16*t1+16*t3+345,3)),8*t1+315),16*t2+14),16*t3+14),8*t1-8*t2+8*t3+315);t5++) {
        for (t6=max(max(max(16*t2,t5+1),-16*t1+16*t2+2*t5-615),-16*t1+16*t2-16*t3+3*t5-330);t6<=min(min(min(16*t2+15,t5+300),-16*t1+16*t2+2*t5-2),-16*t1+16*t2-16*t3+3*t5-1);t6++) {
          for (t7=max(max(16*t3,t5+1),-16*t1+16*t2+3*t5-t6-315);t7<=min(min(16*t3+15,t5+300),-16*t1+16*t2+3*t5-t6-1);t7++) {
            lbv=max(t5+1,-16*t1+16*t2+4*t5-t6-t7-15);
            ubv=min(t5+300,-16*t1+16*t2+4*t5-t6-t7);
            #pragma ivdep
            #pragma vector always
            for (t8=lbv;t8<=ubv;t8++) {
              S1((t1-t2),t2,t3,0,t5,(-t5+t6),(-t5+t7),(-t5+t8));
            }
          }
        }
      }
    }
  }
}
/* End of CLooG code */
