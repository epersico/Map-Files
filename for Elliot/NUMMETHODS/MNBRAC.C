#include <math.h>
#define SHFT(a,b,c,d) (a)=(b);(b)=(c);(c)=(d);
#define GLIMIT 100.0
#define EPS LDBL_EPSILON /* machine precision .. used in ZBRENT */
#define TINY EPS /* really small number; used in mnbrac */
#define GOLD  1.6180339887498948482045868343656381177203091798058

void mnbrac(long double *, long double *, long double *,
	    long double (*)(long double));

void mnbrac(long double *ax, long double *bx, long double *cx, 
	    long double (*func)(long double))
{
  long double ulim, u, r, q, fu, dum;
  long double fa, fb, fc;

  fa = (*func)(*ax); fb = (*func)(*bx);
  if (fb > fa) {
    SHFT(dum, *ax, *bx, dum);
    SHFT(dum, fb, fa, dum);
  }
  *cx = (*bx) + GOLD*(*bx-*ax);
  fc = (*func)(*cx);
  while (fb > fc) {
    r = (*bx - *ax)*(fb-fc);
    q = (*bx - *cx)*(fb-fa);
    u = (*bx) - 
      ((*bx - *cx)*q - (*bx - *ax)*r)/(2.0*SIGN(FMAX(fabsl(q-r), TINY),q-r));
    ulim = (*bx) + GLIMIT*(*cx-*bx);
    if ((*bx-u)*(u-*cx) > 0.0) {
      fu = (*func)(u);
      if (fu < fc) {
	*ax = (*bx);
	*bx = u;
	fa = fb;
	return;
      } else if (fu > fb) {
	*cx = u;
	fc = fu;
	return;
      }
      u = (*cx) + GOLD*(*cx-*bx);
      fu = (*func)(u);
    } else if ((*cx-u)*(u-ulim) > 0.0) {
      fu = (*func)(u);
      if (fu < fc) {
	SHFT(*bx, *cx, u, *cx+GOLD*(*cx-*bx));
	SHFT(fb, fc, fu, (*func)(u));
      }
    } else if ((u-ulim)*(ulim-*cx) >= 0.0) {
      u = ulim;
      fu = (*func)(u);
    } else {
      u = (*cx) + GOLD*(*cx-*bx);
      fu = (*func)(u);
    }
    SHFT(*ax, *bx, *cx, u);
    SHFT(fa, fb, fc, fu);
  }
}
