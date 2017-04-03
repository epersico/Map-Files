#include <math.h>
#define FACTOR 1.6
#define NTRY 50
#define INTERROR 123456
#define TOOBIG=0.1L/LDBL_EPSILON;
#define UPPERBOUND TOOBIG;
#define LOWERBOUND TOOBIG;

int zbrac(long double (*func)(long double), long double *,  long double *);

int zbrac(long double (*func)(long double), long double *x1,  long double *x2)
{
  
  int j;
  long double f1, f2;

  if (*x1 == *x2) (*x2 += 10.0L*LDBL_EPSILON);
  if (*x1 > *x2) { f1=*x1; *x2=*x1; *x1=f1; }
  f1=(*func)(*x1);
  f2=(*func)(*x2);
  for (j=1;j<=NTRY;j++) {
    if (f1*f2 < 0.0) return 1;
    if (fabsl(f1) < fabsl(f2))
      f1=(*func)(*x1 += FACTOR*(*x1 - *x2));
    else
      f2=(*func)(*x2 += FACTOR*(*x2 - *x1));
    if (fabsl(*x2)>TOOBIG || (*x1<LOWERBOUND && *x2>UPPERBOUND)) break;
  }
  return 0;
}
