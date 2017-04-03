#include <math.h>
#define JMAX 40
#define ERROR 123456789.0
long double rtbis(long double (*)(long double), long double, long double, long double);

long double rtbis(long double (*func)(long double), long double x1, long double x2, long double xxacc)
{
  int j;
  long double dx, f, fmid, xmid, rtb;

  f=(*func)(x1);
  fmid=(*func)(x2);
  if (f*fmid >= 0.0) return ERROR; /* root not bracketed */
  rtb = f < 0.0 ? (dx=x2-x1, x1) : (dx=x1-x2, x2); 
  for (j=1;j<=JMAX;j++) {
    fmid=(*func)(xmid=rtb+(dx *= 0.5L)); 
    if (fmid <= 0.0) rtb=xmid;
    if (fabsl(dx) < xxacc || fmid == 0.0) return rtb;
  }
  return ERROR; /* Too many bisections in rtbis */
}
