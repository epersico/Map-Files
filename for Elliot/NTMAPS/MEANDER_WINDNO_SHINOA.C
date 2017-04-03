#include <stdio.h>
#include <math.h>
#define UNDEFINED -123456789.0L
#define REFINE 10

const char *infilename="meander-windno.in", 
  *outfilename="/tmp/meander-windno-shinoa.out";

long double   a, b, amin, amax, astep, bmin, bmax, bstep, cutoff;
int      n, dn, sms, speedup;
char     data_file[60];

void     makeinitials(void);
void     find_om(long double, long double, long double *, int *, long double *, int *, 
	     long double *, int *, long double *, int *);
long double   shino(long double, int);

main(int argc, char **argv)
{
  const long double pi=4.0L*atanl(1.0L);
  int i, ndom=0, nomax=0, nomin=0, ncutoff=0, skip, solutionq, solutionq0;
  int ncutoffo, ndomo, nomaxo, nomino;
  long double dom=0.0L, omega=0.0L, omax=0.0L, omin=0.0L, atmp;
  long double ao, bo, omegao, domo, omaxo, omino;
  FILE *fl, *flbc;
  char dummy[200];

  makeinitials();
  
  /* preparing output */
  fl=fopen(outfilename, "w");
  fprintf(fl, "# n     = %20d \n", n);
  fprintf(fl, "# dn    = %20d \n", dn);
  fprintf(fl, "# speedy= %20d \n", speedup);
  fprintf(fl, "# cutoff= %20.15Lf \n", cutoff);
  fprintf(fl, "# shino = %20d \n", sms);
  fprintf(fl, "# amin  = %20.15Lf \n", amin);
  fprintf(fl, "# amax  = %20.15Lf \n", amax);
  fprintf(fl, "# astep = %20.15Lf \n", astep);
  fprintf(fl, "# bmin  = %20.15Lf \n", bmin);
  fprintf(fl, "# bmax  = %20.15Lf \n", bmax);
  fprintf(fl, "# bstep = %20.15Lf \n", bstep);
  fprintf(fl, "# file  = %s \n", data_file);
  fprintf(fl, "#%15s %15s %15s %10s %10s %10s %10s %10s %10s %10s\n",
	  "a", "b", "om", "ncut", "dom", "@ ndom", 
	  "omax", "@ nomax", "omin", "@ nomin");


  if (astep<0) {
    flbc=fopen(data_file, "r");
    fgets(dummy, 200, flbc);
    fgets(dummy, 200, flbc);
    fgets(dummy, 200, flbc);
    fgets(dummy, 200, flbc);
    fgets(dummy, 200, flbc);
    fgets(dummy, 200, flbc);
    fgets(dummy, 200, flbc);
    do {
      fscanf(flbc,"%Lf %Lf %s", &a, &b, dummy);
      fgets(dummy, 200, flbc);
    } while ((b<bmin) || (a<amin));
  }
  else {
    a=amin;
    b=bmin;
  }
  
  while ((a<=amax) && (b<=bmax))
    {
      solutionq0=1;
      skip=0;
      while (a<=amax)
	{
	  if (skip<speedup) {
	    ao=a; bo=b; omegao=omega; ncutoffo=ncutoff; 
	    domo=dom; ndomo=ndom; omaxo=omax; nomaxo= nomax; omino=omin; nomino=nomin;
	    find_om(shino(a, sms), 0.0L, &omega, &ncutoff, &dom, &ndom, 
		    &omax, &nomax, &omin, &nomin);
	    if (omega<UNDEFINED+1.0) {skip+=1; solutionq=0;}
	    else {skip=0; solutionq=1;}
	    if (solutionq!=solutionq0) {
		atmp=a; solutionq=solutionq0; a=a-astep;
		for (i=1;i<=REFINE;i++) {
		    a+=astep/((long double)(REFINE));
		    find_om(shino(a, sms), 0.0L, &omega, &ncutoff, &dom, &ndom, 
			    &omax, &nomax, &omin, &nomin);
		    if (omega<UNDEFINED+1.0) solutionq=0;
		    else solutionq=1;
		    if (solutionq!=solutionq0) {
			if (solutionq==1)
			    fprintf(fl, 
			    " %15.12Lf %15.12Lf %15.8Lg %10d %10.3Le %10d %10.7Lf %10d %10.7Lf %10d\n",
			    a, b, omega, ncutoff, dom, ndom, omax, nomax, omin, nomin);
			else  
			    fprintf(fl, 
			    " %15.12Lf %15.12Lf %15.8Lg %10d %10.3Le %10d %10.7Lf %10d %10.7Lf %10d\n",
			    ao, bo, -omegao, ncutoffo, domo, ndomo, omaxo, nomaxo, omino, nomino);
			fflush(fl);
			solutionq0=solutionq;
		    }
		    ao=a; bo=b; omegao=omega; ncutoffo=ncutoff; 
		    domo=dom; ndomo=ndom; omaxo=omax; nomaxo= nomax; omino=omin; nomino=nomin;
		}
		a=atmp;
	    }
	  }
	  a+=astep;
	}
      
      b+=bstep;
      a=amin;
    }
    
  if (amax<0) fclose(flbc);
  fclose(fl);
} 


void makeinitials(void)      /* read initial values */
{
  FILE      *fl;
  char      dummy[80];

  /* read initials */
  fl=fopen(infilename, "r");

  fscanf(fl,"%ld", &n); fgets(dummy, 80, fl);        /* max. number of iterations */
  fscanf(fl,"%ld", &dn); fgets(dummy, 80, fl);       /* # of iter. used for conv. */
  fscanf(fl,"%ld", &speedup); fgets(dummy, 80, fl);  /* # of iter. used for conv. */
  fscanf(fl,"%Lg", &cutoff); fgets(dummy, 80, fl);   /* dom convergence crit. */
  fscanf(fl,"%ld", &sms); fgets(dummy, 80, fl);      /* shinohara point 1 or 2 */
  fscanf(fl,"%Lf", &amin); fgets(dummy, 80, fl);     /* amin */
  fscanf(fl,"%Lf", &amax); fgets(dummy, 80, fl);     /* amax */
  fscanf(fl,"%Lf", &astep); fgets(dummy, 80, fl);    /* astep */
  fscanf(fl,"%Lf", &bmin); fgets(dummy, 80, fl);     /* bmin */
  fscanf(fl,"%Lf", &bmax); fgets(dummy, 80, fl);     /* bmax */
  fscanf(fl,"%Lf", &bstep); fgets(dummy, 80, fl);    /* bstep */
  if (astep<0.0) astep=-astep;
  if (astep<=0.0 && bstep<=0.0) astep=1.0L;
  
  fclose(fl);
}

void find_om(long double x0, long double y0, long double *om, 
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
}


long double shino(long double aval, int shpt) 
{
  switch (shpt) {
  case 1: return aval/2.0L+0.25L;
  case 2: return aval/2.0L-0.25L;
  default: return aval/2.0L+0.25L;
  }
}
