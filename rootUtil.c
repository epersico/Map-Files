#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#define FACTOR 1.6
#define NTRY 50
#define ITMAX 100
#define EPS LDBL_EPSILON
#define SHFT(a,b,c,d) (a)=(b);(b)=(c);(c)=(d);
#define FMAX(a,b) (((a)>(b))?(a):(b))
#define SHFT2(a,b,c) (a)=(b);(b)=(c);
#define GOLD  1.6180339887498948482045868343656381177203091798058
#define CGOLD 0.3819660112501051517954131656343618822796908201942
#define DABS(d) (((d)<0)?-1.0*(d):(d))
#define SIGN(a,b) ((b)>=0.0 ? DABS(a):-1.0*DABS(b))
#define TINY EPS /* really small number; used in mnbrac */
#define ZEPS EPS
#define GLIMIT 50.0
#define ERROR 12345678

/*
Using Brent's method, find the root of a function func known to lie between x1 and x2. The
root, returned as zbrent, will be refined until its accuracy is tol
*/
//Fixed a b c and orbit find the y root
long double zbrentY(long double (*func)(long double, long double, long double, long double, int),
				   long double x1, long double x2, long double tol,
				   long double param1, long double param2, long double param3, int param4)
{
	int iter;
	long double a=x1,b=x2,c=x2,d,e,min1,min2;
	long double fa=(*func)(a, param1, param2, param3, param4),fb=(*func)(b, param1, param2, param3, param4),fc,p,q,r,s,tol1,xm;

	if((fa>0.0 && fb >0.0)||(fa<0.0 && fb<0.0)){
		fprintf(stderr,"Root must be bracketed : zbrentY\n");
		//exit(EXIT_FAILURE);
		return ERROR;
	}
	
	fc=fb;
	for(iter=1; iter<=ITMAX; iter++){
		if((fb>0.0 && fc>0.0)||(fb<0.0 && fc<0.0)){
			c=a;
			fc=fa;
			e=d=b-a;
		}
		if (DABS(fc)<DABS(fb)){
			a=b;
			b=c;
			c=a;
			fa=fb;
			fb=fc;
			fc=fa;
		}
		//Convergence check
		tol1=2*EPS*DABS(b)+.5*tol;
		xm=.5*(c-b);
		if(DABS(xm) <= tol1 || fb == 0.0) return b;
		if(DABS(e) >= tol1 && DABS(fa) > DABS(fb)){
			// attempt inverse quadratic interpolation
			s=fb/fa;
			if (a == c){
				p=2.0*xm*s;
				q=1.0-s;
			}else{
				q=fa/fc;
				r=fb/fc;
				p=s*(2.0*xm*q*(q-r)-(b-a)*(r-1));
				q=(q-1.0)*(r-1.0)*(s-1.0);
			}
			//Check whether in bounds
			if(p>0.0) q=-q;
			p = DABS(p);
			min1=3.0*xm*q-DABS(tol1*q);
			min2=DABS(e*q);
			//Accept interpolation
			if(2.0*p < (min1<min2 ? min1 : min2)){
				e=d;
				d=p/q;
			//Interpolation failed, use bisection
			}else{
				d=xm;
				e=d;
			}
		//bounds decreasing too slowly, use bisection
		}else{
			d=xm;
			e=d;
		}
		a=b;
		fa=fb;
		if(DABS(d) > tol1){
			b+= d;
		}else{
			b+= SIGN(tol1,xm);
		}
		fb=(*func)(b, param1, param2, param3, param4);
	}
	fprintf(stderr, "Maximum number of iterations exceeded in zbrentY\n");
	//exit(EXIT_FAILURE);
	return ERROR;
}

/*
Using Brent's method, find the root of a function func known to lie between x1 and x2. The
root, returned as zbrent, will be refined until its accuracy is tol
*/
//Fixed y a c and orbit find the b root
long double zbrentB(long double (*func)(long double, long double, long double, long double, int),
					long double x1, long double x2, long double tol,
					long double param1, long double param2, long double param3, int param4)
{
	int iter;
	long double a=x1,b=x2,c=x2,d,e,min1,min2;
	long double fa=(*func)(param1, param2, a, param3, param4),fb=(*func)(param1, param2, b, param3, param4),fc,p,q,r,s,tol1,xm;

	if((fa>0.0 && fb >0.0)||(fa<0.0 && fb<0.0)){
		fprintf(stderr,"Root must be bracketed: zbrentB\n");
		//exit(EXIT_FAILURE);
		return ERROR;
	}
	
	fc=fb;
	for(iter=1; iter<=ITMAX; iter++){
		if((fb>0.0 && fc>0.0)||(fb<0.0 && fc<0.0)){
			c=a;
			fc=fa;
			e=d=b-a;
		}
		if (DABS(fc)<DABS(fb)){
			a=b;
			b=c;
			c=a;
			fa=fb;
			fb=fc;
			fc=fa;
		}
		//Convergence check
		tol1=2*EPS*DABS(b)+.5*tol;
		xm=.5*(c-b);
		if(DABS(xm) <= tol1 || fb == 0.0) return b;
		if(DABS(e) >= tol1 && DABS(fa) > DABS(fb)){
			// attempt inverse quadratic interpolation
			s=fb/fa;
			if (a == c){
				p=2.0*xm*s;
				q=1.0-s;
			}else{
				q=fa/fc;
				r=fb/fc;
				p=s*(2.0*xm*q*(q-r)-(b-a)*(r-1));
				q=(q-1.0)*(r-1.0)*(s-1.0);
			}
			//Check whether in bounds
			if(p>0.0) q=-q;
			p = DABS(p);
			min1=3.0*xm*q-DABS(tol1*q);
			min2=DABS(e*q);
			//Accept interpolation
			if(2.0*p < (min1<min2 ? min1 : min2)){
				e=d;
				d=p/q;
			//Interpolation failed, use bisection
			}else{
				d=xm;
				e=d;
			}
		//bounds decreasing too slowly, use bisection
		}else{
			d=xm;
			e=d;
		}
		a=b;
		fa=fb;
		if(DABS(d) > tol1){
			b+= d;
		}else{
			b+= SIGN(tol1,xm);
		}
		fb=(*func)(param1, param2,b,param3,param4);
	}
	fprintf(stderr, "Maximum number of iterations exceeded in zbrentB\n");
	//exit(EXIT_FAILURE);
	return ERROR;
}


/*
Given a function func, and given distinct initial points ax and bx, this routine searches in
the downhill direction (defined by the function as evaluated at the intial points) and returns
new points ax, bx, and cx that bracket a minimum of the function.
*/
//Fixed a b c and orbit find the y min brac
void mnbracY(long double *ax, long double *bx, long double *cx,
			long double (*func)(long double, long double, long double, long double, int),
			long double param1, long double param2, long double param3, int param4)
{
  long double ulim, u, r, q, fu, dum;
  long double fa, fb, fc;

  fa = (*func)(*ax, param1, param2, param3, param4); fb = (*func)(*bx, param1, param2, param3, param4);
  //switch roles of a and b so that wwe can go downhill in the direction from a to b
  if (fb > fa) {
    SHFT(dum, *ax, *bx, dum);
    SHFT(dum, fb, fa, dum);
  }
  //first guess for c
  *cx = (*bx) + GOLD*(*bx-*ax);
  fc = (*func)(*cx, param1, param2, param3, param4);
  //keep returning here until we bracket
  while (fb > fc) {
    //compute u by parabolic extrapolation from a, b, c.
    //TINY is used to prevent any possible division by zero
    r = (*bx - *ax)*(fb-fc);
    q = (*bx - *cx)*(fb-fa);
    u = (*bx) - 
      ((*bx - *cx)*q - (*bx - *ax)*r)/(2.0*SIGN(FMAX(DABS(q-r), TINY),q-r));
    ulim = (*bx) + GLIMIT*(*cx-*bx);
    //We will not go farther than this. Test possibities
	if ((*bx-u)*(u-*cx) > 0.0) {
	  //Parabolic u is between b and c
      fu = (*func)(u,param1, param2, param3, param4);
	  //Got a minimum between b and c
      if (fu < fc) {
		*ax = (*bx);
		*bx = u;
		fa = fb;
		fb = fu;
		return;
	  //Got a minimum between a and u
      } else if (fu > fb) {
		*cx = u;
		fc = fu;
		return;
      }
	  //Parabolic fit was no use. Use default magnification.
      u = (*cx) + GOLD*(*cx-*bx);
      fu = (*func)(u, param1, param2, param3, param4);
    //prabolicfit is between c and its allowed limit
    } else if ((*cx-u)*(u-ulim) > 0.0) {
      fu = (*func)(u, param1, param2, param3, param4);
      if (fu < fc) {
		SHFT(*bx, *cx, u, *cx+GOLD*(*cx-*bx));
		SHFT(fb, fc, fu, (*func)(u, param1, param2, param3, param4));
      }
	  //limit prabolic u to maximum allowed value.
    } else if ((u-ulim)*(ulim-*cx) >= 0.0) {
      u = ulim;
      fu = (*func)(u, param1, param2, param3, param4);
	  //reject parabolic u, use default magnification.
    } else {
      u = (*cx) + GOLD*(*cx-*bx);
      fu = (*func)(u, param1, param2, param3, param4);
    }
	//eliminate oldest point and continue
    SHFT(*ax, *bx, *cx, u);
    SHFT(fa, fb, fc, fu);
  }
}

/*
Given a function func, and given distinct initial points ax and bx, this routine searches in
the downhill direction (defined by the function as evaluated at the intial points) and returns
new points ax, bx, and cx that bracket a minimum of the function.
*/
//Fixed y a c and orbit find the b min brac
void mnbracB(long double *ax, long double *bx, long double *cx,
			long double (*func)(long double, long double, long double, long double, int),
			long double param1, long double param2, long double param3, int param4)
{
  long double ulim, u, r, q, fu, dum;
  long double fa, fb, fc;

  fa = (*func)(param1, param2, *ax, param3, param4); fb = (*func)(param1, param2, *bx, param3, param4);
  //switch roles of a and b so that wwe can go downhill in the direction from a to b
  if (fb > fa) {
    SHFT(dum, *ax, *bx, dum);
    SHFT(dum, fb, fa, dum);
  }
  //first guess for c
  *cx = (*bx) + GOLD*(*bx-*ax);
  fc = (*func)(param1, param2, *cx, param3, param4);
  //keep returning here until we bracket
  while (fb > fc) {
    //compute u by parabolic extrapolation from a, b, c.
    //TINY is used to prevent any possible division by zero
    r = (*bx - *ax)*(fb-fc);
    q = (*bx - *cx)*(fb-fa);
    u = (*bx) - 
      ((*bx - *cx)*q - (*bx - *ax)*r)/(2.0*SIGN(FMAX(DABS(q-r), TINY),q-r));
    ulim = (*bx) + GLIMIT*(*cx-*bx);
    //We will not go farther than this. Tes possibities
	if ((*bx-u)*(u-*cx) > 0.0) {
	  //Parabolic u is between b and c
      fu = (*func)(param1, param2,u,param3, param4);
	  //Got a minimum between b and c
      if (fu < fc) {
		*ax = (*bx);
		*bx = u;
		fa = fb;
		fb = fu;
		return;
	  //Got a minimum between a and u
      } else if (fu > fb) {
		*cx = u;
		fc = fu;
		return;
      }
	  //Parabolic fit was no use. Use default magnification.
      u = (*cx) + GOLD*(*cx-*bx);
      fu = (*func)(param1, param2,u,param3, param4);
    //prabolicfit is between c and its allowed limit
    } else if ((*cx-u)*(u-ulim) > 0.0) {
      fu = (*func)(param1, param2,u,param3,param4);
      if (fu < fc) {
		SHFT(*bx, *cx, u, *cx+GOLD*(*cx-*bx));
		SHFT(fb, fc, fu, (*func)(param1, param2,u,param3,param4));
      }
	  //limit prabolic u to maximum allowed value.
    } else if ((u-ulim)*(ulim-*cx) >= 0.0) {
      u = ulim;
      fu = (*func)(param1, param2,u,param3, param4);
	  //reject parabolic u, use default magnification.
    } else {
      u = (*cx) + GOLD*(*cx-*bx);
      fu = (*func)(param1, param2,u,param3,param4);
    }
	//eliminate oldest point and continue
    SHFT(*ax, *bx, *cx, u);
    SHFT(fa, fb, fc, fu);
  }
}

/*
Given a function f, and given a bracketing triplet of abscissas ax, bx, cx (such that bx is
between ax and cx, and f(bx) is less than both f(ax) and f(cx), this routine isolates the minimum
to a fraction precision of about tol using Brent's method. The abscissa of the minimum
is returned as xmin, and the minimum function value is returned as brent, the returned fuction value
*/
//Fixed a b c and orbit find the y min
long double brentY(long double ax, long double bx, long double cx,
				  long double (*f)(long double, long double, long double, long double, int),
				  long double tol, long double *xmin, long double param1, long double param2, long double param3, int param4)
{
  int iter;
  long double a,b,d,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm,e=0.0;
  //a and b mut be in ascending order but input abscissas need not be
  a=(ax < cx ? ax : cx);
  b=(ax > cx ? ax : cx);
  x=w=v=bx;
  fw=fv=fx=(*f)(x, param1, param2, param3, param4);
  //Main program loop
  for (iter=1; iter<=ITMAX; iter++) {
    xm=0.5*(a+b);
    tol2= 2.0*(tol1=tol*DABS(x)+ZEPS);
    //test for done here
	if (DABS(x-xm) <= (tol2-0.5*(b-a))) {
      *xmin = x;
      return fx;
    }
	// Construct a trial parabolic fit
    if (DABS(e) > tol1) {
      r = (x-w)*(fx-fv);
      q = (x-v)*(fx-fw);
      p = (x-v)*q - (x-w)*r;
      q = 2.0*(q-r);
      if (q > 0.0) p = -p;
      q = DABS(q);
      etemp=e;
      e=d;
	  //The condition determine tthe acceptability of the parabolic fit. 
	  //Here we take the golden section step into the larger of the two segments.
      if (DABS(p) >= DABS(0.5*q*etemp) || p <= q*(a-x) || p >= q*(b-x))
	d = CGOLD*(e=(x >= xm ? a-x : b-x));
      else {
		//Take the parabolic step
		d = p/q;
		u = x+d;
		if (u-a < tol2 || b-u < tol2)
			d = SIGN(tol1,xm-x);
	  }
	} else {
      d = CGOLD*(e=(x >= xm ? a-x : b-x));
    }
    u = (DABS(d) >= tol1 ? x+d : x+SIGN(tol1,d));    
	fu = (*f)(u, param1, param2, param3, param4); 
	//This is the one function evaluation per iteration
	//Now decide what to do with our function evaluation 
    if (fu <= fx) {
      if (u >= x) a = x; else b = x;
      SHFT(v,w,x,u);
      SHFT(fv,fw,fx,fu);
    } else {
      if (u < x) a = u; else b = u;
      if (fu <= fw || w == x) {
	v=w;
	w=u;
	fv=fw;
	fw=fu;
      } else if (fu <= fv || v == x || v == w) {
	v = u;
	fv = fu;
      }
    }
  }
  *xmin = x;
	fprintf(stderr, "Maximum number of iterations exceeded in brent");
	return ERROR;
}

/*
Given a function f, and given a bracketing triplet of abscissas ax, bx, cx (such that bx is
between ax and cx, and f(bx) is less than both f(ax) and f(cx), this routine isolates the minimum
to a fraction precision of about tol using Brent's method. The abscissa of the minimum
is returned as xmin, and the minimum function value is returned as brent, the returned fuction value
*/
//Fixed y a c and orbit find the b min
long double brentB(long double ax, long double bx, long double cx,
				  long double (*f)(long double, long double, long double, long double, int),
				  long double tol, long double *xmin, long double param1, long double param2, long double param3, int param4)
{
  int iter;
  long double a,b,d,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm,e=0.0;
  //a and b must be in ascending order but input abscissas need not be
  a=(ax < cx ? ax : cx);
  b=(ax > cx ? ax : cx);
  x=w=v=bx;
  fw=fv=fx=(*f)(param1, param2, x, param3, param4);
  //Main program loop
  for (iter=1; iter<=ITMAX; iter++) {
    xm=0.5*(a+b);
    tol2= 2.0*(tol1=tol*DABS(x)+ZEPS);
    //test for done here
	if (DABS(x-xm) <= (tol2-0.5*(b-a))) {
      *xmin = x;
      return fx;
    }
	// Construct a trial parabolic fit
    if (DABS(e) > tol1) {
      r = (x-w)*(fx-fv);
      q = (x-v)*(fx-fw);
      p = (x-v)*q - (x-w)*r;
      q = 2.0*(q-r);
      if (q > 0.0) p = -p;
      q = DABS(q);
      etemp=e;
      e=d;
	  //The condition determine tthe acceptability of the parabolic fit. 
	  //Here we take the golden section step into the larger of the two segments.
      if (DABS(p) >= DABS(0.5*q*etemp) || p <= q*(a-x) || p >= q*(b-x))
	d = CGOLD*(e=(x >= xm ? a-x : b-x));
      else {
		//Take the parabolic step
		d = p/q;
		u = x+d;
		if (u-a < tol2 || b-u < tol2)
			d = SIGN(tol1,xm-x);
	  }
	} else {
      d = CGOLD*(e=(x >= xm ? a-x : b-x));
    }
    u = (DABS(d) >= tol1 ? x+d : x+SIGN(tol1,d));    
	fu = (*f)(param1, param2, u, param3, param4); 
	//This is the one function evaluation per iteration
	//Now decide what to do with our function evaluation 
    if (fu <= fx) {
      if (u >= x) a = x; else b = x;
      SHFT(v,w,x,u);
      SHFT(fv,fw,fx,fu);
    } else {
      if (u < x) a = u; else b = u;
      if (fu <= fw || w == x) {
	v=w;
	w=u;
	fv=fw;
	fw=fu;
      } else if (fu <= fv || v == x || v == w) {
	v = u;
	fv = fu;
      }
    }
  }
  *xmin = x;
	fprintf(stderr, "Maximum number of iterations exceeded in brent");
	return ERROR;
}


/*
Given a function fx defined on the interval from x1-x2 subdivide the interval into n equally
spaced segments, and search for zero crossings of the function. nb is input as the maximum number
of roots sought, and is reset to the number of bracketing paris xb[1 .. nb], xb2[1 .. nb]
that are found
*/
void zbrakY(long double (*fx)(long double, long double, long double, long double, int), long double x1, long double x2, int n,
			long double xb1[], long double xb2[], int *nb, long double param1, long double param2, long double param3, int param4)
{
	int nbb, i;
	long double x, fp, fc, dx;

	nbb = 0;
	dx = (x2-x1)/n;
	fp=(*fx)(x=x1, param1, param2, param3, param4);
	for(i=1; i<=n; i++){
		fc = (*fx)(x+=dx, param1, param2, param3, param4);
		if(fc*fp < 0.0){
			xb1[++nbb]=x-dx;
			xb2[nbb]=x;
			if(*nb==nbb) return;
		}
		fp =fc;
	}
	*nb = nbb;
}

/*
Given a function fx defined on the interval from x1-x2 subdivide the interval into n equally
spaced segments, and search for zero crossings of the function. nb is input as the maximum number
of roots sought, and is reset to the number of bracketing paris xb[1 .. nb], xb2[1 .. nb]
that are found
*/
void zbrakB(long double (*fx)(long double, long double, long double, long double, int), long double x1, long double x2, int n,
			long double xb1[], long double xb2[], int *nb, long double param1, long double param2, long double param3, int param4)
{
	int nbb, i;
	long double x, fp, fc, dx;

	nbb = 0;
	dx = (x2-x1)/n;
	fp=(*fx)(param1, param2, x=x1, param3, param4);
	for(i=1; i<=n; i++){
		fc = (*fx)(param1, param2, x+=dx, param3, param4);
		if(fc*fp < 0.0){
			xb1[++nbb]=x-dx;
			xb2[nbb]=x;
			if(*nb==nbb) return;
		}
		fp =fc;
	}
	*nb = nbb;
}

/*
Given a function func and an initial guessed range x1 to x2, the routine expands the range
geometrically in one direction until a root is bracketed by the returned values x1 and x2 (in which case zbrac returns 1)
or until the range becomes unacceptably large (in which case zbrac returns 0)
*/
int zbracY(long double (*func)(long double, long double, long double, long double, int), long double *x1, long double *x2,
		   long double param1, long double param2, long double param3, int param4, long double dir)
{
	int j;
	long double f1, f2;

	if(*x1 == *x2) *x2 += TINY;
	f1 = (*func)(*x1, param1, param2, param3, param4);
	f2 = (*func)(*x2, param1, param2, param3, param4);
	for(j=1; j<=NTRY; j++){
		if(f1*f2 < 0.0) return 1;
		if(dir<0.0){
			f1=(*func)(*x1+= dir, param1, param2, param3, param4);
		}else{
			f2 = (*func)(*x2 += dir, param1, param2, param3, param4);
		}
	}
	return 0;
}

/*
Given a function func and an initial guessed range x1 to x2, the routine expands the range
geometrically until a root is bracketed by the returned values x1 and x2 (in which case zbrac returns 1)
or until the range becomes unacceptably large (in which case zbrac returns 0)
*/
int zbracB(long double (*func)(long double, long double, long double, long double, int), long double *x1, long double *x2,
		   long double param1, long double param2, long double param3, int param4)
{
	int j;
	long double f1, f2;

	if(*x1 == *x2) *x2 += TINY;
	f1 = (*func)(param1, param2,*x1, param3, param4);
	f2 = (*func)(param1, param2,*x2, param3, param4);
	for(j=1; j<=NTRY; j++){
		if(f1*f2 < 0.0) return 1;
		if(DABS(f1)<DABS(f2))
			f1=(*func)(param1, param2,*x1+= FACTOR*(*x1-*x2), param3, param4);
		else
			f2=(*func)(param1, param2,*x2 += FACTOR*(*x2-*x1),param3, param4);
	}
	return 0;
}
