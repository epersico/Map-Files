#include <stdio.h>
#include <math.h>
#define ROOT_NO 10000
#include <string.h>
#define NL(flnm,dumch) {do dumch=fgetc(flnm); while (dumch!='\n' && dumch!=EOF);}

const char *infilename="findallorbits.in", outfilename[4][180];

void makeinitials(void);
void prepoutfiles(void);
void ntmapstep(long double *, long double *, int *);
void invntmapstep(long double *, long double *, int *);
void ntdm(long double, long double, long double[2][2]);
long double findx0(long double);
void assignshino(int, long double *, long double *);
long double rootfunceven1(long double);
long double rootfunceven2(long double);
long double rootfuncodd1(long double);
long double rootfuncodd2(long double);
int wind(long double, long double);
long double residue(int, int *, int *, long double *);
 
void zbrak(long double (*)(long double), long double, long double, int, long double[], long double[], int *);
long double rtbis(long double (*)(long double), long double, long double, long double);

long double pi;
long double a0, b0, acc, eqd, brakl0, braku0;
long double xorbs[ROOT_NO], yorbs[ROOT_NO];
int  sln, n, m, brakn, outfiles;
int lift=0;
void (*map)(long double *, long double *, int *);
FILE *fl1, *fl2, *fl3;

main(int argc, char **argv)
{
  long double pery, phalf;
  long double brackl[ROOT_NO], bracku[ROOT_NO];
  int maxbrack=(ROOT_NO-4)/4;
  int i, j, meql, meqlold=0, idorb, idorbiter;
  long double (*rootfunc)(long double);

  /* read input */
  makeinitials();

  /* prepare output files */
  fl1=fopen(outfilename[0], "w");
  if (outfiles>1) fl2=fopen(outfilename[1], "w");
  if (outfiles>2) fl3=fopen(outfilename[2], "w");
  prepoutfiles();
    
  /* main calculation */
  meql=0; 
  fprintf(fl1,"# a=%15.11Lf, b=%15.11Lf\n", a0, b0);
  if (outfiles>1) fprintf(fl2,"# a=%15.11Lf, b=%15.11Lf\n", a0, b0);
  if (outfiles>1) fprintf(fl2,"# %12s %14s\n", "q0", "p0");
  if (outfiles>2) fprintf(fl3,"# a=%15.11Lf, b= %15.11Lf\n", a0, b0);
  if (outfiles>2) fprintf(fl3,"# %12s %14s\n", "q0", "p0");

  /* limits of shino point orbits */
  map=&ntmapstep;
  for (sln=1;sln<=4;sln++) {
      fprintf(fl1,"# shp=%2d, ", sln);  /* sln, number of roots found */
  }

  /*  limits of shino point inverse orbits */
  map=&invntmapstep;
  for (sln=1;sln<=4;sln++) {
      fprintf(fl1,"# shp=-%2d, ", sln);  /* sln, number of roots found */
  }
  


  map=&ntmapstep;
  for (sln=1;sln<=4;sln++) {
      meqlold=meql;
      /* find root number and brackets */
      maxbrack=(ROOT_NO-8)/4;
      if (2*(n/2)==n) {                       /* n even */
	  if (sln<=2) rootfunc=&rootfunceven1;      /* I0 -> I0 after n/2 */
	  else rootfunc=&rootfunceven2;             /* I1 -> I1 after n/2 */
      }
      else {                                  /* n odd */
	  if (sln<=2) rootfunc=&rootfuncodd1;       /* I0 -> I1 after n/2+1=(n+1)/2 */
	  else rootfunc=&rootfuncodd2;              /* I1 -> I0 after n/2=(n-1)/2 */
      }

      zbrak(rootfunc, brakl0, braku0, brakn, &brackl[0], &bracku[0], &maxbrack);
      fprintf(fl1,"# sln=%2d, maxbrack=%4d, ", sln, maxbrack);  /* sln, number of roots found */
  
      /* examine root brackets */
      for (i=1;i<=maxbrack;i++) {
	  
	  /* find roots */
	  pery=rtbis(rootfunc, brackl[i], bracku[i], acc);
	  
	  /* count and store ?/n roots with m/n winding numbers */
	  lift=wind(findx0(pery), pery); /* optional output to fl3 */
	  if (lift==m) {xorbs[++meql]=findx0(pery);yorbs[meql]=pery;}
      }
  
      fprintf(fl1,"meql=%4d\n", meql-meqlold);
      for (i=meqlold+1;i<=meql;i++) {
	  fprintf(fl1,"%15.11Lf %15.11Lf %10.3Lg ", 
		  xorbs[i], yorbs[i], residue(i, &idorb, &idorbiter, &phalf));
	  fprintf(fl1,"%10.5Lf %6d %6d %6d\n", phalf, i, idorb, idorbiter);
      }
  
      fprintf(fl1,"\n");
      if (outfiles>1) fprintf(fl2,"\n");
      if (outfiles>2) fprintf(fl3,"\n");
  }
    

  fclose(fl1);
  if (outfiles>1) fclose(fl2);
  if (outfiles>2) fclose(fl3);
  exit(0); 
} 


void makeinitials(void)      /* read and adjust initial values */
{
  FILE *fl;
  char dummy;
  int i;

  pi=4.0L*atanl(1.0L);

  /* read initials */
  fl=fopen(infilename, "r");
  
  fscanf(fl,"%Lf", &a0); NL(fl, dummy);       /* start a */
  fscanf(fl,"%Lf", &b0); NL(fl, dummy);       /* start b */
  fscanf(fl,"%ld", &m); NL(fl, dummy);        /* m */
  fscanf(fl,"%ld", &n); NL(fl, dummy);        /* period n */
  fscanf(fl,"%Le", &acc); NL(fl, dummy);      /* absolute root acc */
  fscanf(fl,"%Le", &eqd); NL(fl, dummy);      /* absolute root acc */
  fscanf(fl,"%Lf", &brakl0); NL(fl, dummy);   /* lower end of interval for "brak"ket */
  fscanf(fl,"%Lf", &braku0); NL(fl, dummy);   /* upper end of interval for "brak"ket */
  fscanf(fl,"%ld", &brakn); NL(fl, dummy);    /* divisions for "brak"ket interval */
  fscanf(fl,"%ld", &outfiles); NL(fl, dummy); /* number of output files */

  for (i=0;i<outfiles;i++) {
      fscanf(fl,"%s", &outfilename[i][0]); NL(fl, dummy);  /* bifc input filename */
  }
  
  fclose(fl);

  /* adjust initials */
  if (brakl0==braku0) {brakl0=-0.5L; braku0=0.5L;}
}


void prepoutfiles(void)
{
  fprintf(fl1,"# a:       %10.7Lf , b:      %10.7Lf \n", a0, b0); 
  fprintf(fl1,"# omega:   (m =%4d) / (n =%4d) \n", m, n); 
  fprintf(fl1,"# brakl0:  %10.7Lf , braku0:  %10.7Lf , brakn:    %10d\n", 
	  brakl0, braku0, brakn); 
  fprintf(fl1,"# acc:     %10.3Le , eqd:     %10.3Le\n", acc, eqd); 
  fprintf(fl1,"# %13s %15s %10s %10s %6s %6s %6s\n", 
	  "x0", "y0", "res", "yhalf", "orb#", "eq to", "@ i="); 

  if (outfiles>1) {
    fprintf(fl2,"# a:       %10.7Lf , b:      %10.7Lf \n", a0, b0); 
    fprintf(fl2,"# omega:   (m =%4d) / (n =%4d) \n", m, n); 
    fprintf(fl2,"# brakl0:  %10.7Lf , braku0:  %10.7Lf , brakn:    %10d\n", 
	    brakl0, braku0, brakn); 
    fprintf(fl2,"# acc:     %10.3Le , eqd:     %10.3Le\n", acc, eqd); 
  }
  
  if (outfiles>2) {
    fprintf(fl3,"# a:       %10.7Lf , b:      %10.7Lf \n", a0, b0); 
    fprintf(fl3,"# omega:   (m =%4d) / (n =%4d) \n", m, n); 
    fprintf(fl3,"# brakl0:  %10.7Lf , braku0:  %10.7Lf , brakn:    %10d\n", 
	    brakl0, braku0, brakn); 
    fprintf(fl3,"# acc:     %10.3Le , eqd:     %10.3Le\n", acc, eqd); 
  }
}

/* single iteration of NT-map */
void ntmapstep(long double *x, long double *y, int *lft)
{
  *y = *y - b0*sinl(2.0L*pi*(*x));
  *x = *x + a0*(1.0L-(*y)*(*y));
  
  /* x modulo 1 */
  while (*x> 0.55L) {*x=*x-1.0L; (*lft)++;}
  while (*x< -0.45L) {*x=*x+1.0L; (*lft)--;}
}

/* single iteration of inverse NT-map */
void invntmapstep(long double *x, long double *y, int *lft)
{
  *x = *x - a0*(1.0L-(*y)*(*y));
  *y = *y + b0*sinl(2.0L*pi*(*x));
    
  /* x modulo 1 */
  while (*x> 0.55L) {*x=*x-1.0L; (*lft)--;}
  while (*x< -0.45L) {*x=*x+1.0L; (*lft)++;}
}

/* functional determinant of NT-map */
void ntdm(long double x, long double y, long double dxy[2][2])
{
  long double tmp1, tmp2;

  tmp1=-2.0L*a0*(y-b0*sinl(2.0L*pi*x));
  tmp2=-2.0L*pi*b0*cosl(2.0L*pi*x);
  dxy[0][0]=1.0L+tmp1*tmp2;
  dxy[0][1]=tmp1;
  dxy[1][0]=tmp2;
  dxy[1][1]=1.0L;
}


/* find x0(y0) so that point lies on symmetry line */
long double findx0(long double y0)
{
  long double tmp;

  switch (sln) {
    case 1: {tmp=0.0L; break;}
    case 2: {tmp=0.5L; break;}
    case 3: {tmp=a0/2.0L*(1.0L-y0*y0); break;}
    case 4: {tmp=a0/2.0L*(1.0L-y0*y0)+0.5L; break;}
    default: {tmp=0.0L; printf("Uhoh - sln...\n");}
  }
   
  /* x modulo 1 */
  while (tmp> 0.55L) tmp=tmp-1.0L;
  while (tmp< -0.45L) tmp=tmp+1.0L;
  
  return tmp; 
}

/* assign shinohara points */
void assignshino(int ss, long double *xval, long double *yval)
{
  switch (sln) {
      case 1: {*xval=-0.25L; *yval=-0.5L*b0; break;}
      case 2: {*xval=a0/2.0L-0.25L; *yval=0.0L; break;}
      case 3: {*xval=0.25L; *yval=0.5L*b0;break;}
      case 4: {*xval=a0/2.0L-0.25L; *yval=0.0L; break;}
      default: {*xval=-0.25L; *yval=-0.5L*b0; printf("Uhoh - shinopoint...\n"); break;}
  }
}

    
/* sin(2 pi x) after n/2 map iterations for root search */
long double rootfunceven1(long double y0)
{
  int i;
  long double xx, yy;
  xx=findx0(y0); yy=y0; lift=0;
  for (i=1;i<=n/2;i++) map(&xx, &yy, &lift);
  return sinl(2.0L*pi*xx); /* I0 -> I0 after n/2 */
}

/* sin(2 pi x) after n/2 map iterations for root search */
long double rootfunceven2(long double y0)
{
  int i;
  long double xx, yy;

  xx=findx0(y0); yy=y0; lift=0;
  for (i=1;i<=n/2;i++) map(&xx, &yy, &lift);
  return sinl(2.0L*pi*(xx-a0/2.0L*(1.0L-yy*yy))); /* I1 -> I1 after n/2 */
}

/* sin(2 pi x) after n/2 map iterations for root search */
long double rootfuncodd1(long double y0)
{
  int i;
  long double xx, yy;

  xx=findx0(y0); yy=y0; lift=0;
  for (i=1;i<=n/2+1;i++) map(&xx, &yy, &lift);    /* n odd: n/2+1=(n+1)/2 */
  return sinl(2.0L*pi*(xx-a0/2.0L*(1.0L-yy*yy))); /* I0 -> I1 after n/2+1=(n+1)/2 */
}

/* sin(2 pi x) after n/2 map iterations for root search */
long double rootfuncodd2(long double y0)
{
  int i;
  long double xx, yy;

  xx=findx0(y0); yy=y0; lift=0;
  for (i=1;i<=n/2;i++) map(&xx, &yy, &lift); /* n odd: n/2=(n-1)/2 */
  return sinl(2.0L*pi*xx);                   /* I1 -> I0 after n/2=(n-1)/2 */
}


/* find limiting fixed point of shinohara orbit */
int findshinofp(int shp, long double *shx, long double *shy, long double *shlx, long double *shly)
{
    long double xx0, yy0, xx, yy, dxx, dyy;
    int i, j, lift;

    assignshino(shp, &xx0, &yy0);
    xx=xx0; yy=yy0;

    for (i=1;i<=9*n;i++) map(&xx,&yy,&lift);
    xx0=xx; yy0=yy;
    for (i=1;i<=n;i++) map(&xx,&yy,&lift);
    dxx=xx-xx0; dyy=yy-yy0;
    *shlx=exp(log(dxx)/((long double) (10)));
    *shly=exp(log(dyy)/((long double) (10)));
    xx0=xx; yy0=yy;
    for (i=1;i<=n;i++) map(&xx,&yy,&lift);
    j=11;
    
    while (fabs(xx-xx0)<=fabs(ddx)+acc && fabs(yy-yy0)<=fabs(dyy)+acc && fabs(yy-yy0)>acc && 
	fabs(*shlx-exp(log(xx-xx0)/((long double) (j))))>lacc &&
	fabs(*shly-exp(log(yy-yy0)/((long double) (j))))>lacc && j<n*n*n) {
	dxx=xx-xx0; dyy=yy-yy0;
	*shlx=exp(log(dxx)/((long double) (j)));
	*shly=exp(log(dyy)/((long double) (j)));
	xx0=xx; yy0=yy;
	for (i=1;i<=n;i++) map(&xx,&yy,&lift);
	j++;
    }
    *shx=xx; *shy=yy;
    if (fabs(yy-yy0)<acc) {
	*shx=xx; 
	*shy=yy;	
	*shlx=exp(log(dxx)/((long double) (j)));
	*shly=exp(log(dyy)/((long double) (j)));
	return 1;
    }
    else if (fabs(xx-xx0)<=fabs(ddx)+acc && fabs(yy-yy0)<=fabs(dyy)+acc && j<n*n*n) {
	*shlx=exp(log(dxx)/((long double) (j)));
	*shly=exp(log(dyy)/((long double) (j)));
	xx=xx0+pow(*shlx,((long double) (j)))/(1.0L-(*shlx));
	yy=yy0+pow(*shly,((long double) (j)))/(1.0L-(*shly));
	xx0=xx; yy0=yy;
	for (i=1;i<=n;i++) map(&xx,&yy,&lift);
	if (fabs(yy-yy0)<acc) {*shx=xx; *shy=yy; return 1;} else {
	    while (fabs(xx-xx0)<=fabs(ddx)+acc && fabs(yy-yy0)<=fabs(dyy)+acc 
		   && fabs(yy-yy0)>acc && j<n*n*n) {
		dxx=xx-xx0; dyy=yy-yy0;
		*shlx=exp(log(dxx)/((long double) (j)));
		*shly=exp(log(dyy)/((long double) (j)));
		xx0=xx; yy0=yy;
		for (i=1;i<=n;i++) map(&xx,&yy,&lift);
		j++;
	    }
	    if (fabs(yy-yy0)<acc) {
		*shx=xx; 
		*shy=yy;	
		*shlx=exp(log(dxx)/((long double) (j)));
		*shly=exp(log(dyy)/((long double) (j)));
		return 1;
	    }
	    else return 0;
	}
	
    }
    else return 0;    
}

/* find winding number of known periodic orbit, with optional output to fl3 */
int wind(long double q, long double p)
{
  long double qq=q, pp=p;
  int lft=0, i;

  lft=0;
  for (i=1;i<=n;i++) {if (outfiles>2) fprintf(fl3,"%14.8Lf %14.8Lf", qq, pp); map(&qq, &pp, &lft);}
  if (outfiles>2) fprintf(fl3,"%14.8Lf%14.8Lf%5d%10.2Le%10.2Le\n",qq,pp,lft,qq-q, pp-p);
  return lft;
}


/* find residues of a known periodic orbit, with opt. output to fl2 and *tmp on I0 */
long double residue(int ii, int *iieq, int *iiter, long double *tmp)
{
  long double qq, pp, dpq[2][2], lm[2][2], lm0[2][2];
  int i, j, k, lft=0;

  *iieq=ii;
  *iiter=0;
  qq=xorbs[ii];
  pp=yorbs[ii];
  ntdm(qq,pp,lm);

  for (k=1;k<n;k++) {

    if (outfiles>1) fprintf(fl2,"%14.8Lf %14.8Lf\n", qq, pp); 

    map(&qq, &pp, &lft);
    if (k==n/2) *tmp=pp;
    ntdm(qq,pp,dpq);
    for (i=0;i<=1;i++) for (j=0;j<=1;j++) lm0[i][j]=lm[i][j]; 
    for (i=0;i<=1;i++) for (j=0;j<=1;j++) lm[i][j]=dpq[i][0]*lm0[0][j]+dpq[i][1]*lm0[1][j];
    i=1; while (*iieq>=ii && i<ii) {
	if (fabs(pp-yorbs[i])<eqd && fabs(qq-xorbs[i])<eqd) {*iieq=i; *iiter=k;}
	i++;
    }    
  }

  if (outfiles>1) fprintf(fl2,"%14.8Lf %14.8Lf ", qq, pp);
  map(&qq, &pp, &lft);
  fprintf(fl2,"%14.4Le %14.4Le\n\n", qq-findx0(yorbs[i]), pp-yorbs[i]);

  return (2.0L-lm[0][0]-lm[1][1])/4.0L;
}

/* zbrak from Num. Rec. */
void zbrak(long double (*fx)(long double), long double x1, long double x2, int nn, 
	   long double xb1[], long double xb2[], int *nb)
{
  int nbb, i;
  long double x, fp, fc, dx;
  
  nbb=0;
  dx=(x2-x1)/nn;
  fp=(*fx)(x=x1);
  for (i=1;i<=nn;i++) {
    fc=(*fx)(x += dx);
    if (fc*fp <= 0.0) {
      xb1[++nbb]=x-dx;
      xb2[nbb]=x;
      if (*nb==nbb) return;
    }
    fp=fc;
  }
  *nb = nbb;
}


/* minimally modified version of rtbis from Num. Rec. */
long double rtbis(long double (*func)(long double), long double x1, long double x2, long double xacc)
{
  const int JMAX=40;
  int nrerror;
  int j;
  long double dx, f, fmid, xmid, rtb;
  
  f=(*func)(x1);
  fmid=(*func)(x2);
  if (f*fmid >= 0.0) nrerror=1; /* "Root must be bracketed for bisection in rtbis" */
  rtb = f < 0.0 ? (dx=x2-x1, x1) : (dx=x1-x2, x2); 
  for (j=1;j<=JMAX;j++) {
    fmid=(*func)(xmid=rtb+(dx *= 0.5L)); 
    if (fmid <= 0.0) rtb=xmid;
    if (fabs(dx) < xacc || fmid == 0.0) return rtb;
  }
  nrerror=2; /* "Too many bisections in rtbis" */
  return 0.0L;
}
