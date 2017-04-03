
/* --**-- per-orbs-on-bif-curve.v1.c --**-- */
/* this is identical to per-orbs-on-bif-curve.c escept the follwing:
   i for s1 symmetry line, i find all the orbits and hence the yrange
   bounded by yup and ydown decreases with each period. but for the other
   symmetry lines, i search s3 for m=odd and s2 for even/odd.
   which means that for let's say [35] on s2, i'm not using the yrange of
   [33] on s2 but i'm using the yrange of [30] on s2. and thus it's larger
   than it needs to be. hence, this program adjusts that by assigning the
   yrange correctly for the next symmetry line... does it make sense???
   probably not....
*/

#include <stdio.h>
#include <math.h>
#include <time.h> /* **** ...ONLY FOR TIMING PURPOSES... **** */
#include <limits.h> /* for limits of int and char types like MAX_INT etc. */
#include <float.h> /* for limits of floats like FLT_MAX etc. */

#define pi    3.1415926535897932384626433832795028841971693993748L
#define f(y) a*(1-(y)*(y))
#define df(y) (-2.0L)*a*(y)
#define symm_line(y) dreml((0.5*f(y)*((symm_line_type-1)/2) + 0.5*((symm_line_type-1)%2)),1.0L)
#define up 1
#define down 0
#define ERROR 11111111111111.0   
#define ITMAX 1000
#define EPS LDBL_EPSILON /* machine precision .. used in ZBRENT */
#define ztol EPS /* accuracy of root returned by zbrent */ 

#define mod4(x) ( (x) <= 4 ? (x) : ((x)-4) )

long double zbrent(long double (*)(long double), long double, long double);
long double F_iter(long double);
unsigned int getorbit(long double, long double *, long double *, long double *);
int find_mid_point(long double, long double *, long double *,
		    long double *, long double *, long double *, int *);

unsigned int m, iter;
int period, symm_line_type;
long double a, k;
unsigned int numerator[50]   = {0};
unsigned int denominator[50] = {0};
int sml[2][43] = {{0},{0}};

int main(int argc, char* argv[])
{
  int i, max_i, j, found, l, tmp1,tmp2, tmp0;
  unsigned int sml_mid, m_mid;
  long double ystep, ymin, ymax, y0, y1, f0, f1, root, xerr, yerr, Res_mid;
  long double res, yfound[2], x_err[2], y_err[2], Res[2], xmid, ymid; 
  long double xorb, xerrmid, yerrmid, dummy;
  long double yorb[2][3]={{-1.1,-1.1,-1.1},{1.1,1.1,1.1}}, ystep_min=1.0e-14;
  FILE *fptr, *fout;

  scanf("%Lf%Lf%d", &a, &k, &max_i);
  if ( (fptr=fopen(argv[1], "r")) == NULL) {
    fprintf(stderr, "*** error opering in file  %s  ***\n", argv[1]);
    exit(1);
  }

  if ( (fout=fopen(argv[2], "w")) == NULL) {
    fprintf(stderr, "*** error opering out file  %s  ***\n", argv[2]);
    exit(1);
  }

  nice(15);
  for (i = 0; i <= max_i ; i++) {
    fscanf(fptr, "%d%d%d%Lf\n", &tmp0,&tmp1,&tmp2,&dummy);
    numerator[2*i]=tmp1; denominator[2*i]=tmp2; 
    if (numerator[2*i]%2==0) sml[1][2*i]=2; else sml[1][2*i]=3;
    fscanf(fptr, "%d%d%d%Lf\n", &tmp0,&tmp1,&tmp2,&dummy);
    m = numerator[2*i]; period = denominator[2*i]; 
    for (j=0; j<2; j++) {
      symm_line_type = sml[j][2*i] + (j+1)%2; found = 0;
      ymin = yorb[down][symm_line_type-1]; ymax = yorb[up][symm_line_type-1]; 
      ystep = (ymax-ymin)/100.0L;
      y0 = ymin; y1 = y0+ystep;
      f0 = F_iter(y0); f1 = F_iter(y1);
      while (found < 2) {
	if (f0*f1 <= 0) {
	  root = zbrent(F_iter, y0, y1);
	  if (getorbit(root, &xerr, &yerr, &res) == m) {
	    yfound[found] = root; x_err[found] = xerr;
	    y_err[found] = yerr; Res[found] = res; found += 1;
	  }
	}
	y0 = y1; f0 = f1; y1 += ystep; 
	if (y1 > ymax && found < 2) {
	  ystep = ystep/10.0L;
	  if (ystep < ystep_min) exit(0);
	  y0 = ymin; y1 = y0 + ystep;
	  f0 = F_iter(y0); f1 = F_iter(y1);
	} else f1 = F_iter(y1);
      }
      for (l=0; l<2; l++) {
	yorb[l][symm_line_type-1] = yfound[l]; xorb = symm_line(yfound[l]);
	fprintf(fout, "%9d: %9d: %25.18Le: %25.18Le: %d: %d: %25.18Le: %7.3Lf: %25.18Le: %7.3Lf: %27.18Le: %8.3Lf\n", 
	       m, period, a, k, symm_line_type, l, xorb, 
	       log10l(fabsl(x_err[l])), yorb[l][symm_line_type-1], 
	       log10l(fabsl(y_err[l])), Res[l], log10l(fabsl(Res[l])));
	sml_mid = find_mid_point(yorb[l][symm_line_type-1], &xmid, &ymid,
				 &xerrmid, &yerrmid, &Res_mid, &m_mid);
	if (period%2==0 && j==0) yorb[l][1]=ymid;
	if (m%2==0 && j==0) yorb[l][2]=ymid;
	fprintf(fout, "%9d: %9d: %25.18Le: %25.18Le: %d: %d: %25.18Le: %7.3Lf: %25.18Le: %7.3Lf: %27.18Le: %8.3Lf\n", 
	       m_mid, period, a, k, sml_mid, l, xmid, log10l(fabsl(xerrmid)),
	       ymid, log10l(fabsl(yerrmid)), Res_mid, log10l(fabsl(Res_mid)));
	fflush(fout);
      }
    }
  }
  fclose(fptr);
  fclose(fout);
}

long double F_iter(long double y)
{
  long double xold, yold, xnew, ynew, tmp1;
  int i;

  if ((period % 2) == 1) {
    if ((symm_line_type == 1) || (symm_line_type == 2))
      iter = (period + 1)/2;
    else iter = (period - 1)/2;
  } else iter = period/2;
  yold = y; xold = symm_line(y);

  for (i=1; i <= iter; i++) {
    ynew = yold - k*sinl(2*pi*xold); xnew = dreml((xold + f(ynew)),1.0L);
    xold = xnew; yold = ynew;
  }

  if ((period % 2) == 0) {
    if ((symm_line_type == 1) || (symm_line_type == 2))
      return sinl(2*pi*xnew);
    else return sinl(2*pi*(xnew - 0.5*f(ynew)));
  } else {
    if ((symm_line_type == 1) || (symm_line_type == 2))
      return sinl(2*pi*(xnew - 0.5*f(ynew)));
    else return sinl(2*pi*xnew);
  }
}

unsigned int getorbit(long double y, long double *xerr, 
		      long double *yerr, long double *res)
{
  int i; unsigned int tmp;
  long double x0, xold, yold, xnew, ynew, xold_u, xnew_u, wind_num;
  long double M11, M12, M21, M22, m11, m12, m21, m22;

  x0 = symm_line(y);
  yold = y; xold = x0; xold_u = xold;
  M11 = 1.0 - 2*pi*k*cosl(2*pi*xold)*df(yold);
  M12 = df(yold);
  M21 = (-2.0)*pi*k*cosl(2*pi*xold);
  M22 = 1;
  for (i=1; i <= period; i++) {
    ynew = yold - (k)*sinl(2*pi*xold); xnew = dreml((xold + f(ynew)),1.0L);
    xnew_u = xold_u + f(ynew);
    xold = xnew; xold_u = xnew_u; yold = ynew;

    if (i != period) {
      m11 = M11*(1.0 - 2*pi*k*cos(2*pi*xold)*df(yold)) + 
	M12*((-2.0)*pi*k*cos(2*pi*xold));
      m12 = M11*df(yold) + M12;
      m21 = M21*(1.0 - 2*pi*k*cos(2*pi*xold)*df(yold)) +
	M22*((-2.0)*pi*k*cos(2*pi*xold));
      m22 = M21*df(yold) + M22;
      M11 = m11; M12 = m12; M21 = m21; M22 = m22;
    }
  }

  *res = (1.0/4)*(2.0-(M11+M22)); *yerr = ynew - y;
  if ((xnew - x0) > 0.5) *xerr = xnew - x0 - 1.0L;
  else if ((xnew - x0) <= -0.5) *xerr = xnew - x0 + 1.0L;
  else *xerr = xnew - x0;

  if ((xnew_u-x0) >= 0.0) {
    tmp = (unsigned int) (xnew_u-x0);
    if (((xnew_u-x0)-tmp) >= 0.5)
      return (tmp+1);
    else
      return tmp;
  } else {
    tmp = (unsigned int) -(xnew_u-x0);
    if ((-(xnew_u-x0)-tmp) >= 0.5)
      return -(tmp+1);
    else
      return -tmp;
  }
}

int find_mid_point(long double y, long double *xm, long double *ym,
		   long double *xme, long double *yme, 
		   long double *resm, int *mm)
{
  int outsml, tmpsml, iter, i;
  long double xold, yold, xnew, ynew;

  if ((m % 2 == 0) && (period % 2 == 1))
    outsml = mod4(symm_line_type + 2);
  else if (m % 2 == 1)
    if (period % 2 == 1)
      outsml = mod4(symm_line_type + 2*(symm_line_type%2) + 1);
    else 
      outsml = mod4(symm_line_type + 2*(symm_line_type%2) - 1);
  
  if ((period % 2) == 1) {
    if ((symm_line_type == 1) || (symm_line_type == 2))
      iter = (period + 1)/2;
    else iter = (period - 1)/2;
  } else iter = period/2;

  yold = y; xold = symm_line(y);

  for (i=1; i <= iter; i++) {
    ynew = yold - (k)*sinl(2*pi*xold); xnew = dreml((xold + f(ynew)),1.0L);
    xold = xnew; yold = ynew;
  }

  *xm = xnew; *ym = ynew;
  tmpsml = symm_line_type;
  symm_line_type = outsml;
  *mm = getorbit(ynew, xme, yme, resm);
  symm_line_type = tmpsml;
  return outsml;
}

long double zbrent(long double (*func)(long double), long double x1, 
		   long double x2)
{
  int i;
  long double a=x1,b=x2,c=x2,d,e,min1,min2;
  long double fa=(*func)(a),fb=(*func)(b),fc,p,q,r,s,tol1,xm;

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
    fb=(*func)(b);
  }
  return ERROR;
}
