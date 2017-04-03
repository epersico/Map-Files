#include <stdio.h>
#include <math.h>

const char *infilename="ntstmap.in";
char outfilename[160];

long double   a, b, jj0, djj, th0, dth;
int      n, skip, jjloops, thloops;

void     makeinitials(void);
void     mapstep(long double *, long double *, long double *);

main(int argc, char **argv)
{
  const long double pi=4.0L*atan(1.0L);
  int ci, cj, ck, pr, prt;
  long double jj, th, thlift, th0sh;
  FILE *fl;

  makeinitials();
  
  /* preparing output */
  fl=fopen(outfilename, "w");
  fprintf(fl, "# b  = %30.15Lf \n", a);
  fprintf(fl, "# a  = %30.15Lf \n", b);
  fprintf(fl, "# n  = %30d \n", n);
  fprintf(fl, "# dj = %30.15Lf, dth = %30.15Lf \n", djj, dth);
  fprintf(fl, "# nj = %30d, nth = %30d \n", jjloops, thloops);

  printf("----------- P r o g r e s s ------------\n");
  pr=0;
  prt=n*2*(jjloops/2)*thloops+4*n;

  /* shinohara points */
  for (ci=1; ci<=4; ci++) {
      switch (ci) {
	  case 1: {jj=-b/2.0L; th=-0.25L; break;}
	  case 2: {jj=0.0L; th=a/2.0L-0.25L; break;}
	  case 3: {jj=b/2.0L; th=0.25L; break;}
	  case 4: {jj=0.0L; th=a/2.0L+0.25L; break;}
	  default: {jj=-b/2.0L; th=-0.25L;}
      }
      
      while (th> 1.0) th=th-1.0L;
      while (th< 0.0) th=th+1.0L;
      thlift=th; th0sh=th;
      
      fprintf(fl, "\n# j0 = %10.5Lf, th0 = %10.5Lf \n", jj, th);
      fprintf(fl, "#%15s %15s %15s %15s %15s\n","step", "j", "th", "thlift", "wind");
      fprintf(fl, " %15d %15.5Lf %15.5Lf %15.5Lf %15.10Lf\n", 0, jj, th, thlift, thlift-th);
	
      for (ck=1; ck<=n; ck++)
      {
	  mapstep(&jj,&th, &thlift);
	  if ((ck/skip)*skip==ck)
	      fprintf(fl, " %15d %15.5Lf %15.5Lf %15.5Lf %15.10Lf\n", 
		      ck, jj, th, thlift, (thlift-th0sh)/((double)(ck)));
	  
	  /* progress report */
	  pr++;
	  if (floor(40.0*((double)(pr+1))/((double)(prt)))
	      >floor(40.0*((double)(pr))/((double)(prt))))
	  {
	      printf("*"); fflush(stdout);
	  }
      }
  }
  

  for (ci=0; ci<jjloops/2; ci++) {        /* several initial conditions */
    for (cj=0; cj<thloops; cj++) {        /* several initial conditions */

      /* initializations */
      jj=jj0+ci*djj;
      th=th0+cj*dth;
      /* th modulo 1 */
      while (th> 1.0) th=th-1.0L;
      while (th< 0.0) th=th+1.0L;
      /* jj modulo 1
      while (jj> 1.0) jj=jj-1.0;
      while (jj< 0.0) jj=jj+1.0; */
      thlift=th;      

      /* prepare output files */
      fprintf(fl, "\n# j0 = %10.5Lf, th0 = %10.5Lf \n", jj, th);
      fprintf(fl, "#%15s %15s %15s %15s %15s\n","step", "j", "th", "thlift", "wind");
      fprintf(fl, " %15d %15.5Lf %15.5Lf %15.5Lf %15.10Lf\n", 0, jj, th, thlift, thlift-th);
    
      /* main loop */
      for (ck=1; ck<=n; ck++)
	{
	  mapstep(&jj,&th, &thlift);
	  if ((ck/skip)*skip==ck)
	      fprintf(fl, " %15d %15.5Lf %15.5Lf %15.5Lf %15.10Lf\n", 
		      ck, jj, th, thlift, (thlift-th0-cj*dth)/((double)(ck)));

	  /* progress report */
	  pr++;
	  if (floor(40.0*((double)(pr+1))/((double)(prt)))
	      >floor(40.0*((double)(pr))/((double)(prt))))
	    {
	      printf("*"); fflush(stdout);
	    }
	}

      /* initializations */
      jj=jj0-ci*djj;
      th=th0+cj*dth;
      /* th modulo 1 */
      while (th> 1.0) th=th-1.0L;
      while (th< 0.0) th=th+1.0L;
      /* jj modulo 1
      while (jj> 1.0) jj=jj-1.0;
      while (jj< 0.0) jj=jj+1.0; */
      thlift=th;      

      /* prepare output files */
      fprintf(fl, "\n# j0 = %10.5Lf, th0 = %10.5Lf \n", jj, th);
      fprintf(fl, "#%15s %15s %15s %15s %15s\n","step", "j", "th", "thlift", "wind");
      fprintf(fl, " %15d %15.5Lf %15.5Lf %15.5Lf %15.10Lf\n", 0, jj, th, thlift, thlift-th);
    
      /* main loop */
      for (ck=1; ck<=n; ck++)
	{
	  mapstep(&jj,&th, &thlift);
	  if ((ck/skip)*skip==ck)
	      fprintf(fl, " %15d %15.5Lf %15.5Lf %15.5Lf %15.10Lf\n", 
		      ck, jj, th, thlift, (thlift-th0-cj*dth)/((double)(ck)));

	  /* progress report */
	  pr++;
	  if (floor(40.0*((double)(pr+1))/((double)(prt)))
	      >floor(40.0*((double)(pr))/((double)(prt))))
	    {
	      printf("*"); fflush(stdout);
	    }
	}
    }
  }
  printf("\n"); fflush(stdout);
  fclose(fl);
} 


void makeinitials(void)      /* read initial values */
{
  FILE      *fl;
  char      dummy[80];
  long double tmp;

  /* read initials */
  fl=fopen(infilename, "r");

  fscanf(fl,"%Lf", &a); fgets(dummy, 80, fl);        /* a */
  fscanf(fl,"%Lf", &b); fgets(dummy, 80, fl);        /* b */
  fscanf(fl,"%ld %ld", &n, &skip); fgets(dummy, 80, fl); /* # steps, skip output */
  fscanf(fl,"%Lf", &jj0); fgets(dummy, 80, fl);      /* j0 */
  fscanf(fl,"%Lf", &th0); fgets(dummy, 80, fl);      /* th0 */
  fscanf(fl,"%Lf", &djj); fgets(dummy, 80, fl);      /* dj */
  fscanf(fl,"%Lf", &dth); fgets(dummy, 80, fl);      /* dth */
  fscanf(fl,"%ld", &jjloops); fgets(dummy, 80, fl);  /* number of j loops */
  fscanf(fl,"%ld", &thloops); fgets(dummy, 80, fl);  /* number of th loops */
  fscanf(fl,"%s", outfilename);                     /* outfile */

  if (djj==0.0L) {djj=1.0L/((long double)(jjloops));}
  if (dth==0.0L) {dth=1.0L/((long double)(thloops));}
  tmp=th0;
  if (th0<-1.0L) {jj0=-b/2.0L; tmp=-0.25L;}
  if (th0<-2.0L) {jj0=0.0L; tmp=a/2.0L-0.25L;}
  if (th0<-3.0L) {jj0=b/2.0L; tmp=0.25L;}
  if (th0<-4.0L) {jj0=0.0L; tmp=a/2.0L+0.25L;}
  if (skip<1) skip=1;
  th0=tmp;
  fclose(fl);
}


void mapstep(long double *y, long double *x, long double *xx)
{
  const long double pi=4.0L*atan(1.0L);

  *y = *y - b*sin(2.0L*pi*(*x));
  *x = *x + a*(1.0L-(*y)*(*y));
  *xx = *xx + a*(1.0L-(*y)*(*y));
  
  /* x modulo 1 */
  while (*x> 1.0) *x=*x-1.0L;
  while (*x< 0.0) *x=*x+1.0L;
  /* y modulo 1
  while (*y> 1.0) *y=*y-1.0;
  while (*y< 0.0) *y=*y+1.0; */
}







