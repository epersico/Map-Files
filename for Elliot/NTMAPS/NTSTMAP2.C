#include <stdio.h>
#include <math.h>
#define NL(flnm,dumch) {do dumch=fgetc(flnm); while (dumch!='\n' && dumch!=EOF);}

const char *infilename="ntstmap2.in";
char outfilename[160];

long double   a, b, jj0[200], th0[200];
int      n, skip, xy_no;

void     makeinitials(void);
void     mapstep(long double *, long double *, long double *);
void     invmapstep(long double *, long double *, long double *);

main(int argc, char **argv)
{
  const long double pi=4.0L*atan(1.0L);
  int ci, cj, ck, pr, prt;
  long double jj, th, jj00, th00, thlift, th0sh;
  FILE *fl;

  makeinitials();
  
  /* preparing output */
  fl=fopen(outfilename, "w");
  fprintf(fl, "# b  = %30.15Lf \n", a);
  fprintf(fl, "# a  = %30.15Lf \n", b);
  fprintf(fl, "# n  = %30d \n", n);

  printf("----------- P r o g r e s s ------------\n");
  pr=0;
  prt=n*2*xy_no;

  for (ci=0; ci<xy_no; ci++) {        /* several initial conditions */

  /* forward iteration */    
      /* initializations */
      jj00=jj0[ci]; jj=jj00;
      th00=th0[ci]; th=th00;
      /* th modulo 1 */
      while (th> 0.5) th=th-1.0L;
      while (th< -0.5) th=th+1.0L;
      thlift=th;      

      /* prepare output files */
      fprintf(fl, "\n# j0 = %10.5Lf, th0 = %10.5Lf \n", jj, th);
      fprintf(fl, "#%15s %15s %15s %15s %15s\n","step", "x", "y", "xlift", "wind");
      fprintf(fl, " %15d %15.15Lf %15.15Lf %15.15Lf %15.10Lf\n", 0, th, jj, thlift, thlift-th00);
    
      /* main loop */
      for (ck=1; ck<=n; ck++)
	{
	  mapstep(&jj,&th, &thlift);
	  if ((ck/skip)*skip==ck)
	      fprintf(fl, " %15d %15.15Lf %15.15Lf %15.15Lf %15.10Lf\n", 
		      ck, th, jj, thlift, (thlift-th00)/((double)(ck)));

	  /* progress report */
	  pr++;
	  if (floor(40.0*((double)(pr+1))/((double)(prt)))
	      >floor(40.0*((double)(pr))/((double)(prt))))
	    {
	      printf("*"); fflush(stdout);
	    }
	}

  /* backward iteration */    
      /* initializations */
      jj=jj00;
      th=th00;
      /* th modulo 1 */
      while (th> 0.5) th=th-1.0L;
      while (th< -0.5) th=th+1.0L;
      thlift=th;      
    
      /* main loop */
      for (ck=1; ck<=n; ck++)
	{
	  invmapstep(&jj,&th,&thlift);
	  if ((ck/skip)*skip==ck)
	      fprintf(fl, " %15d %15.15Lf %15.15Lf %15.15Lf %15.10Lf\n", 
		      -ck, th, jj, thlift, (thlift-th00)/((double)(ck)));

	  /* progress report */
	  pr++;
	  if (floor(40.0*((double)(pr+1))/((double)(prt)))
	      >floor(40.0*((double)(pr))/((double)(prt))))
	    {
	      printf("*"); fflush(stdout);
	    }
	}

  }
  printf("\n"); fflush(stdout);
  fclose(fl);
} 


void makeinitials(void)      /* read initial values */
{
  FILE      *fl;
  char      dummy;
  long double tmp;
  int       i;

  /* read initials */
  fl=fopen(infilename, "r");

  fscanf(fl,"%Lf %Lf", &a, &b); NL(fl, dummy);       /* (a,b) */
  fscanf(fl,"%ld %ld", &n, &skip); NL(fl, dummy);    /* # steps, skip output */
  fscanf(fl,"%s", outfilename); NL(fl, dummy);       /* outfilename */
  fscanf(fl,"%d", &xy_no);  NL(fl, dummy);           /* # (x,y) values below */
  for (i=0;i<xy_no;i++) {
      fscanf(fl,"%Lf %Lf", &th0[i], &jj0[i]); NL(fl, dummy);   /* (x,y) */
      tmp=th0[i];
      if (th0[i]<-1.0L) {jj0[i]=-b/2.0L; tmp=-0.25L;}
      if (th0[i]<-2.0L) {jj0[i]=0.0L; tmp=a/2.0L-0.25L;}
      if (th0[i]<-3.0L) {jj0[i]=b/2.0L; tmp=0.25L;}
      if (th0[i]<-4.0L) {jj0[i]=0.0L; tmp=a/2.0L+0.25L;}
      if (skip<1) skip=1;
      th0[i]=tmp;
  }
  fclose(fl);
}


void mapstep(long double *y, long double *x, long double *xx)
{
  const long double pi=4.0L*atanl(1.0L);

  *y = *y - b*sinl(2.0L*pi*(*x));
  *x = *x + a*(1.0L-(*y)*(*y));
  *xx = *xx + a*(1.0L-(*y)*(*y));
  
  /* x modulo 1 */
  while (*x> 0.5) *x=*x-1.0L;
  while (*x< -0.5) *x=*x+1.0L;
}

void invmapstep(long double *y, long double *x, long double *xx)
{
  const long double pi=4.0L*atanl(1.0L);

  *x = *x - a*(1.0L-(*y)*(*y));
  *xx = *xx - a*(1.0L-(*y)*(*y));
  *y = *y + b*sinl(2.0L*pi*(*x));
  
  /* x modulo 1 */
  while (*x> 0.5) *x=*x-1.0L;
  while (*x< -0.5) *x=*x+1.0L;
}







