#include <stdio.h>
#include <math.h>
#define UNDEFINED -123456789.0

const char *infilename="windno-profile.in", 
  *outfilename="/tmp/windno-profile.out", 
  *extfilename="/tmp/windno-extrema.out";

long double a0[100], b0[100], a, b, dslmax, dslmin, dslstep,
  ymin, ymax, ystep, cutoff;
int n, dn, sl, nab, speedup;

void makeinitials(void);
void find_om(long double, long double, long double *, long double *, int *, long double *, int *, 
	     long double *, int *, long double *, int *);
long double symmline(long double, int);

main(int argc, char **argv)
{
  const long double pi=4.0L*atanl(1.0L);
  int abloop, xcount, nomax, nomin, ndom, ncutoff, skip;
  long double y0, yf, omx, omold, dir, dsl, dom, omega, omax, omin;
  FILE *fl, *flx;
  char dummy[200];

  makeinitials();
  
  /* preparing output */
  fl=fopen(outfilename, "w");
  fprintf(fl, "# n     = %20d \n", n);
  fprintf(fl, "# dn    = %20d \n", dn);
  fprintf(fl, "# speedy= %20d \n", speedup);
  fprintf(fl, "# cutoff= %20.15Lf \n", cutoff);
  fprintf(fl, "# sl    = %20d \n", sl);
  fprintf(fl, "# dslmin= %20.15Lf \n", dslmin);
  fprintf(fl, "# dslmax= %20.15Lf \n", dslmax);
  fprintf(fl, "# dslst = %20.15Lf \n", dslstep);
  fprintf(fl, "# ymin  = %20.15Lf \n", ymin);
  fprintf(fl, "# ymax  = %20.15Lf \n", ymax);
  fprintf(fl, "# ystep = %20.15Lf \n", ystep);

  flx=fopen(extfilename, "w");
  fprintf(flx, "# n     = %20d \n", n);
  fprintf(flx, "# dn    = %20d \n", dn);
  fprintf(flx, "# speedy= %20d \n", speedup);
  fprintf(flx, "# cutoff= %20.15Lf \n", cutoff);
  fprintf(flx, "# sl    = %20d \n", sl);
  fprintf(flx, "# dslmin= %20.15Lf \n", dslmin);
  fprintf(flx, "# dslmax= %20.15Lf \n", dslmax);
  fprintf(flx, "# dslst = %20.15Lf \n", dslstep);
  fprintf(flx, "# ymin  = %20.15Lf \n", ymin);
  fprintf(flx, "# ymax  = %20.15Lf \n", ymax);
  fprintf(flx, "# ystep = %20.15Lf \n#\n", ystep);

  for (abloop=0;abloop<nab;abloop++)
    {
      a=a0[abloop];
      b=b0[abloop];
      
      /* one of Shinohara's meander points */
      find_om(a/2.0L+0.25L, 0.0, &yf, &omega, &ncutoff, &dom, 
	      &ndom, &omax, &nomax, &omin, &nomin);

      fprintf(fl, "# a = %20.15Lf, b = %20.15Lf \n", a, b);
      fprintf(fl, "# meander has om = %20.15Lf after %8d iterations. \n", 
	      omega, ncutoff);
      fprintf(fl, "# %8s %20s %20s %20s %8s %10s %8s %10s %8s %10s %8s\n",
	      "dsl", "y", "om", "yf", "ncut", "dom", "@ ndom", 
	      "omax", "@ nomax", "omin", "@ nomin");

      fprintf(flx, "# a = %20.15Lf, b = %20.15Lf \n", a, b);
      fprintf(flx, "# meander has om = %20.15Lf after %8d iterations. \n", 
	      omega, ncutoff);
      fprintf(flx, "# %8s %20s %20s %20s %10s %8s %10s %8s %10s %8s %10s %8s\n",
	      "dsl", "y", "yf", "om", "om-omx", "ncut", "dom", "@ ndom", 
	      "omax", "@ nomax", "omin", "@ nomin");


      for (dsl=dslmin;dsl<dslmax;dsl+=dslstep) 
	{
	  y0=ymin+ystep/3000*pi; dir=1.0L; xcount=0; 
	  omx=0.0L; omold=0.0L; skip=0;
	  
	  while (y0<=ymax)
	    {
	      if (skip<speedup) {
		find_om(symmline(y0,sl)+dsl, y0, &yf, &omega, &ncutoff, &dom, 
			&ndom, &omax, &nomax, &omin, &nomin);
		if (omega<UNDEFINED+1.0) {
		  skip+=1;
		}
		else {
		  if ((fabs(omega-omx)>cutoff) && (dir*(omega-omold)<0.0)) {
		    fprintf(flx, " %10.7Lf %20.17Lf %20.13Lf %20.13Lg %20.17Lf %8d %10.3Le %8d %10.3Lg %8d %10.3Lg %8d\n", 
			    dsl, y0, omega, yf, omega-omx, ncutoff, dom, ndom, 
			    omax, nomax, omin, nomin);
		    dir *= -1.0L;
		    omx = omega;
		    xcount++; 
		  }
		  skip=0;
		  omold=omega;
		}
		
		fprintf(fl, " %10.7Lf %20.17Lf %20.13Lf %20.13Lg %8d %10.3Le %8d %10.3Lg %8d %10.3Lg %8d\n", 
			dsl, y0, omega, yf, ncutoff, dom, ndom, omax, nomax, omin, nomin);
	  }
	      else
		fprintf(fl, " %10.7Lf %20.17Lf %20.13Lf %20.13Lg %8d %10.3Le %8d %10.3Lg %8d %10.Lg %8d\n", 
			dsl, y0, UNDEFINED, 0.0, 0, UNDEFINED, 0, 0.0, 0, 0.0, 0);	  
	      
	      y0+=ystep;
	    }
	  fprintf(fl, "# %d extrema\n", xcount);
	  fprintf(flx, "# %d extrema\n", xcount);

	}
      fprintf(fl, "\n");
      fprintf(flx, "\n");

    }
  
  fclose(flx);
  fclose(fl);
} 


void makeinitials(void)      /* read initial values */
{
  FILE      *fl;
  int       i;
  char      dummy[80];

  /* read initials */
  fl=fopen(infilename, "r");

  fscanf(fl,"%ld", &n); fgets(dummy, 80, fl);        /* max. number of iterations */
  fscanf(fl,"%ld", &dn); fgets(dummy, 80, fl);       /* # of iter. used for conv. */
  fscanf(fl,"%ld", &speedup); fgets(dummy, 80, fl);  /* # of iter. used for conv. */
  fscanf(fl,"%Lg", &cutoff); fgets(dummy, 80, fl);   /* dom convergence crit. */
  fscanf(fl,"%ld", &sl); fgets(dummy, 80, fl);       /* # of symmetry-line */
  fscanf(fl,"%Lf", &dslmin); fgets(dummy, 80, fl);      /* x-offset from sl */
  fscanf(fl,"%Lf", &dslmax); fgets(dummy, 80, fl);      /* x-offset from sl */
  fscanf(fl,"%Lf", &dslstep); fgets(dummy, 80, fl);      /* x-offset from sl */
  fscanf(fl,"%Lf", &ymin); fgets(dummy, 80, fl);     /* ymin */
  fscanf(fl,"%Lf", &ymax); fgets(dummy, 80, fl);     /* ymax */
  fscanf(fl,"%Lf", &ystep); fgets(dummy, 80, fl);    /* ystep */
  fscanf(fl,"%ld", &nab); fgets(dummy, 80, fl);      /* # of (a,b) pairs */
  
  for (i=0;i<nab;i++) {                              /* read (a,b) pairs */
      fscanf(fl,"%Lf %Lf", &a0[i], &b0[i]); 
      fgets(dummy, 80, fl);
  }
  
  if (ystep<=0.0) {printf("Bad ystep!\n"); exit(2);}
  fclose(fl);
} 

void find_om(long double x0, long double y0, long double *yf, long double *om, 
	     int *ncut, long double *diffom, int *ndiffom, 
	     long double *maxom, int *nmaxom, long double *minom, int *nminom)
{
  const long double pi=4.0L*atanl(1.0L);
  int i;
  long double omo, x=x0, y=y0;

  *ncut=n+1; *diffom=(UNDEFINED); *ndiffom=n+1; 
  *maxom=-1.0e11L; *nmaxom=0;
  *minom=1.0e11L; *nminom=0;  

  y = y - b*sinl(2.0L*pi*x);
  x = x + a*(1.0L-y*y);
  omo = (x-x0);

  i=2;
  do           /* iteration of nt-map */
    {
      y = y - b*sinl(2.0L*pi*x);
      x = x + a*(1.0L-y*y);
      *om = (x-x0)/((long double)(i));
      if (*om>*maxom) {*maxom=*om; *nmaxom=i;}     /* store max. of om */
      if (*om<*minom) {*minom=*om; *nminom=i;}     /* store min. of om */
      if (*ncut<=n) {                              /* check if dom remains small */
	if (fabs(omo-*om)<cutoff) {                /* store max. dom after conv. */
	  if (fabs(omo-*om)>fabs(*diffom)) {*diffom=*om-omo; *ndiffom=i;} 	  
	}
	else {
	  *ncut=n+1;                               /* left convergence range */
	  if (fabs(omo-*om)>fabs(*diffom)) {*diffom=*om-omo; *ndiffom=i;}
	}
      }
      else if (fabs(omo-*om)<cutoff) {             /* entered convergence range */
	*ncut=i;
	*diffom=*om-omo; *ndiffom=i; 
      }
      omo=*om;
      i++;
    }
  while ((i<=n)&&((i<=*ncut+dn)||(i<=*nmaxom+dn)||(i<=*nminom+dn)));
  if ((i>n) && ((i<=*ncut+dn)||(i<=*nmaxom+dn)||(i<=*nminom+dn))) *om=(UNDEFINED);
  *yf=y;
}

long double symmline(long double y, int sln) 
{
  switch (sln) {
  case 1: return 0.0L;
  case 2: return 0.5L;
  case 3: return a/2.0L*(1.0L-y*y);
  case 4: return a/2.0L*(1.0L-y*y)+0.5;  
  default: return 0.0L;
  }
}



