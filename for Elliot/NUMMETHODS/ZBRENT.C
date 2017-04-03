#include <math.h>
#define ERROR 123456789.0   
#define ITMAX 100   /* max number of "bisections" in zbrent and brent */
#define EPS LDBL_EPSILON /* machine precision .. used in ZBRENT */
#define ztol EPS /* accuracy of root returned by zbrent */ 

long double zbrent(long double (*)(long double, long double *, long double *), 
		   long double *, long double *, long double, long double);

long double zbrent(long double (*func)(long double,long double *,long double *)
		   , long double *y_m, long double *range,
		   long double x1, long double x2)
{
  int i;
  long double a=x1,b=x2,c=x2,d,e,min1,min2;
  long double fa, fb, fc,p,q,r,s,tol1,xm;

  fa = (*func)(a, y_m, range); fb = (*func)(b, y_m, range);
  if ((fa > 0.0 && fb > 0.0) || (fa < 0.0 && fb < 0.0)) return ERROR;
  fc=fb;
  for (i=1;i<=ITMAX;i++) {
    if ((fb > 0.0 && fc > 0.0) || (fb < 0.0 && fc < 0.0)) {
      c=a;
      fc=fa;
      e=d=b-a;
    }
    if (fabsl(fc) < fabsl(fb)) {
      a=b;
      b=c;
      c=a;
      fa=fb;
      fb=fc;
      fc=fa;
    }
    tol1=2.0*EPS*fabsl(b)+0.5*ztol;
    xm=0.5*(c-b);
    if (fabsl(xm) <= tol1 || fb == 0.0) return b;
    if (fabsl(e) >= tol1 && fabsl(fa) > fabsl(fb)) {
      s=fb/fa;
      if (a == c) {
	p=2.0*xm*s;
	q=1.0-s;
      } else {
	q=fa/fc;
	r=fb/fc;
	p=s*(2.0*xm*q*(q-r)-(b-a)*(r-1.0));
	q=(q-1.0)*(r-1.0)*(s-1.0);
      }
      if (p > 0.0) q = -q;
      p=fabsl(p);
      min1=3.0*xm*q-fabsl(tol1*q);
      min2=fabsl(e*q);
      if (2.0*p < (min1 <min2 ? min1 : min2)) {
	e=d;
	d=p/q;
      } else {
	d=xm;
	e=d;
      }
    } else {
      d=xm;
      e=d;
    }
    a=b;
    fa=fb;
    if (fabsl(d) > tol1)
      b += d;
    else
      b += (xm > 0.0 ? fabsl(tol1) : -fabsl(tol1));
    fb=(*func)(b,y_m,range);
  }
  return ERROR;
}
