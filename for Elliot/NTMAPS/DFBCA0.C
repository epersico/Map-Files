#include <stdio.h>
#include <math.h>
#include <string.h>
#include <float.h>

#define NL(flnm,dumch) {do dumch=fgetc(flnm); while (dumch!='\n' && dumch!=EOF);}
#define ntmapstep(x,y,a,b) {y-=(b)*sinl(pi2*(x)); x+=(a)-(a)*(y)*(y);}
#define smoothen(x,y) {if ((x)>0.0L) y=((x)/(1.0L+(x))); else y=((x)/(1.0L-(x)));}
#define SHFT(a,b,c,d) (a)=(b);(b)=(c);(c)=(d); /* for mnbrac & brent */
#define SIGN(a,b) ((b) >= 0.0 ? fabsl(a) : -fabsl(a))
#define FMAX(a,b) ((a) > (b) ? (a) : (b))
#define REFINE 10.0L

const char *infilename="dfbca0.in";
const long double SMALL=1.0e-10L;
const int SUPPERROR=15;

void makeinitials(void);
void prepoutfiles(void);
void maincalc(void);
long double rootbyfunc(long double, int, unsigned long int, unsigned long int);
long double rootbynegfunc(long double, int, unsigned long int, unsigned long int);
long double findx0(long double, long double, int);
long double fiter2funcoe0(long double, long double, int, unsigned long int, unsigned long int);
long double fiter2funcoo0(long double, long double, int, unsigned long int, unsigned long int);
long double fiter2funceo0(long double, long double, int, unsigned long int, unsigned long int);
long double fiter2funcoe1(long double, long double, int, unsigned long int, unsigned long int);
long double fiter2funcoo1(long double, long double, int, unsigned long int, unsigned long int);
long double fiter2funceo1(long double, long double, int, unsigned long int, unsigned long int);
void zbrak(long double (*)(long double, long double, int, unsigned long int, unsigned long int), 
	   long double, int, unsigned long int, unsigned long int, 
	   long double, long double, int, long double [], long double [], int *);
long double zbrent(long double (*)(long double, long double, int, unsigned long int, unsigned long int), 
		   long double, int, unsigned long int, unsigned long int, long double, long double);
void mnbrac(long double *, long double *, long double *, 
	    long double (*)(long double, int, unsigned long int, unsigned long int), 
	    int, unsigned long int, unsigned long int);
long double brent(long double, long double, long double, 
	    long double (*)(long double, int, unsigned long int, unsigned long int), 
		  int, unsigned long int, unsigned long int, long double *);
void mkstr(char *, unsigned long int);
int qsortcompare(long double [], long double []);

long double pi2;
long double omega, a0, a00, a0res, a0min, a0max, a0step, b0, b00, ymin, ymax, ystep, y00, bxx;
long double bracdist, bracdistmin;
int sln, slnmin, slnmax, qpi, qp_no, bres, xbranch, outmode;
unsigned long int qq[500], pp[500], itstime2mod;
char iterfilename[160], iterfilenamebase[80];
long double braku[3], brakl[3];

FILE *fl0, *fl1;
long double (*fiterfunc)(long double, long double, int, unsigned long int, unsigned long int);

main(int argc, char **argv)
{
  long double b;

  /* read input */
  makeinitials();
  
  for (qpi=0; qpi<qp_no; qpi++) {
  for (sln=slnmin;sln<=slnmax;sln++) {
      
      /* prepare output files */
      prepoutfiles();
      
      /* function selection for all possible q/p-combos */
      if (2*(pp[qpi]/2)==pp[qpi]) {        /* p even (-> q odd) */
	  if (sln<=2) fiterfunc=&fiter2funcoe0;    /* I_0^i -> I_0^(1-i) after p/2, (q+-1)/2 */
	  else fiterfunc=&fiter2funcoe1;           /* I_1^i -> I_1^(i-1) after p/2, (q+-1)/2 */
      }
      else {                               /* p odd */
	  if (sln<=2) {
	      if (2*(qq[qpi]/2)==qq[qpi])      /* q even */
		  fiterfunc=&fiter2funceo0;        /* I_0^i -> I_1^i after (p+1)/2, q/2 */
	      else fiterfunc=&fiter2funcoo0;       /* I_0^i -> I_1^(1-i) after (p+1)/2, (q+-1)/2 */
	  }
	  else {
	      if (2*(qq[qpi]/2)==qq[qpi])      /* q even */
		  fiterfunc=&fiter2funceo1;        /* I_0^1 -> I_1^1 after (p+1)/2, q/2 */
	      else fiterfunc=&fiter2funcoo1;       /* I_0^1 -> I_1^0 after (p+1)/2, (q+1)/2 */
	  }
      }

      maincalc();

      if (outmode>0) fclose(fl1);      
  }

  }
  
  fclose(fl0);
  exit(0); 
} 


void makeinitials(void)      /* read and adjust initial values */
{
  FILE *fl;
  char dummy, ctmp[60];
  int i;
  long double tmp;
  
  pi2=8.0L*atanl(1.0L);

  /* read initials */
  fl=fopen(infilename, "r");
  
  fscanf(fl,"%d", &sln); NL(fl, dummy);                /* symmetry line */
  fscanf(fl,"%Lf", &a0max);  NL(fl, dummy);            /* a0-range*/
  fscanf(fl,"%Lf %Le", &a0step, &a0res);  NL(fl, dummy);          /* a0-step, a0-res */
  fscanf(fl,"%Lf %Lf", &ymin, &ymax);  NL(fl, dummy);  /* y-range */
  fscanf(fl,"%Lf", &ystep);  NL(fl, dummy);            /* ystep */
  fscanf(fl,"%Lf %Lf", &b00, &bxx);  NL(fl, dummy);    /* min b, max b */
  fscanf(fl,"%d %d", &bres, &xbranch);  NL(fl, dummy); /* # of brackets, # of extrema branches*/
  fscanf(fl,"%Le %Le", &bracdist, &bracdistmin);  NL(fl, dummy);  /* distance of first mn-brac-b */
  fscanf(fl,"%u", &itstime2mod);  NL(fl, dummy);       /* mod intervals threshold */
  fscanf(fl,"%d", &outmode);  NL(fl, dummy);           /* output mode */
  fscanf(fl,"%s", iterfilenamebase); NL(fl, dummy);    /* outfilename */
  fscanf(fl,"%d", &qp_no);  NL(fl, dummy);    /* number of p/q values below */
  for (i=0;i<qp_no;i++) {
      fscanf(fl,"%u %u", &qq[i], &pp[i]); NL(fl, dummy);                     /* q/p */
  }

  a00=0.0L;
  xbranch=1;      /* note: higher orders not implemented yet */
  if (sln==0) {slnmin=1; slnmax=4;}
  else {slnmin=sln; slnmax=sln;}
  if (bxx<b00) {tmp=b00; b00=bxx; bxx=tmp;}
  if (bres<1) bres=1;
  if (xbranch<0) xbranch=0;
  if (bracdistmin>bracdist) bracdistmin=bracdist;
  fclose(fl);

  strcpy(iterfilename,iterfilenamebase);
  strcat(iterfilename,".s");
  strcpy(ctmp,""); mkstr(&ctmp[0],sln); strcat(iterfilename,ctmp);
  
  fl0=fopen(iterfilename, "w");
  fprintf(fl0,"# Symm.l.: %2d \n", sln);
  fprintf(fl0,"# ax:      %10.7Lf \n", a0max); 
  fprintf(fl0,"# astep:   %10.7Lf , ares:    %10.3Le \n", a0step, a0res); 
  fprintf(fl0,"# ymin:    %10.7Lf , ymax:    %10.7Lf \n", ymin, ymax); 
  fprintf(fl0,"# yres:    %10.7Lf \n", ystep); 
  fprintf(fl0,"# b00:     %10.7Lf , bx:      %10.7Lf \n", b00, bxx); 
  fprintf(fl0,"# bres:    %10d , xbranch: %10d\n", bres, xbranch); 
  fprintf(fl0,"# ydist:   %10.3Le , ydistmin:%10.3Le\n", bracdist, bracdistmin); 
  fprintf(fl0,"# mod:     %10u \n", itstime2mod); 
  fprintf(fl0,"# %8s%10s%22s%22s%6s%22s\n\n", "q", "p", "a", "b", "#sol", "y");

}


void prepoutfiles(void)
{
  char ctmp[60];
    
  /* adjust initials */
  omega=((long double)(qq[qpi]))/((long double)(pp[qpi]));
  if (a00<=omega) {
    a0min=omega;
    b0=0.0L;
  }
  else {a0min=a00; b0=b00;}
  if (b0<=0.0L) b0=0.0L; else b0=b00;

  if (outmode>0) {
      strcpy(iterfilename,iterfilenamebase);
      strcpy(ctmp,""); mkstr(&ctmp[0],qq[qpi]); strcat(iterfilename,ctmp);
      strcat(iterfilename,"_");
      strcpy(ctmp,""); mkstr(&ctmp[0],pp[qpi]); strcat(iterfilename,ctmp);
      strcat(iterfilename,".s");
      strcpy(ctmp,""); mkstr(&ctmp[0],sln); strcat(iterfilename,ctmp);
      
      fl1=fopen(iterfilename, "w");
      fprintf(fl1,"# Symm.l.: %2d \n", sln);
      fprintf(fl1,"# q:       %10u , p:       %10u \n", qq[qpi], pp[qpi]);
      fprintf(fl1,"# a0:      %10.7Lf , ax:      %10.7Lf \n", a0min, a0max); 
      fprintf(fl1,"# astep:   %10.7Lf , ares:    %10.3Le \n", a0step, a0res); 
      fprintf(fl1,"# ymin:    %10.7Lf , ymax:    %10.7Lf \n", ymin, ymax); 
      fprintf(fl1,"# yres:    %10.7Lf \n", ystep); 
      fprintf(fl1,"# b0:      %10.7Lf , bx:      %10.7Lf \n", b0, bxx); 
      fprintf(fl1,"# bres:    %10d , xbranch: %10d\n", bres, xbranch); 
      fprintf(fl1,"# ydist:   %10.3Le , ydistmin:%10.3Le\n", bracdist, bracdistmin); 
      fprintf(fl1,"# mod:     %10u \n", itstime2mod); 
      fprintf(fl1,"# %23s%25s%25s%10s\n\n", "a", "b", "y", "# sol");
  }
}

void maincalc() 
{
    long double braca, bracb, bracc, sol[2*xbranch+2][2], 
	savsolb, savsoly, ssavsolb, ssavsoly, a0tmpstep, checker, bracd;
    int xcount, i, j, maxxcount;

    xcount=0; maxxcount=1;

    a0tmpstep=a0step;
    
    a0=a0min; while (a0<=a0max && maxxcount<2*xbranch) {

	if (a0<a0min+a0tmpstep) {
	    braca=sqrtl(1.0L-omega/a0);
	    if (rootbyfunc(braca, sln, qq[qpi], pp[qpi])>0.0L) braca=-braca;
	    bracd=bracdist;
	    bracb=braca-bracd; 
	}
	else {
	    if (xcount>0) {
		if (sol[xcount/2-1][0]-sol[xcount/2][0]<sol[xcount/2][0]-sol[xcount/2+1][0]) {
		    braca=sol[xcount/2-1][0];
		    bracd=(sol[xcount/2-1][0]-sol[xcount/2][0])/3.0L;
		} 
		else {
		    braca=sol[xcount/2+1][0];
		    bracd=(sol[xcount/2+1][0]-sol[xcount/2][0])/3.0L;
		}	    
		if (fabsl(bracd)<bracdistmin) bracd=bracd/fabsl(bracd)*bracdistmin;
		if (fabsl(bracd)>bracdist) bracd=bracd/fabsl(bracd)*bracdist;
		braca-=bracd/2.0L;
	    }
	    else braca=savsoly;
	    bracb=braca-bracd;
      }
	
	xcount=0;
	mnbrac(&braca, &bracb, &bracc, rootbyfunc, sln, qq[qpi], pp[qpi]);
	sol[xcount][1]=-brent(braca, bracb, bracc, rootbyfunc, sln, qq[qpi], pp[qpi], 
			      &sol[xcount][0]);
	
	(xcount)++;
	braca=sol[0][0]; bracb=braca+bracd;
	if ((checker=rootbyfunc(bracb, sln, qq[qpi], pp[qpi]))<0.0L) {
	    mnbrac(&braca, &bracb, &bracc, rootbynegfunc, sln, qq[qpi], pp[qpi]);
	    sol[xcount][1]=brent(braca, bracb, bracc, rootbynegfunc, sln, qq[qpi], pp[qpi], 
				 &sol[xcount][0]);
	    checker=rootbyfunc(sol[xcount][0]+fabs(bracd*(braca-bracb))/(bracb-braca), 
			       sln, qq[qpi], pp[qpi]);
	}
	if (checker>0.0L) (xcount)--;
	
	while (checker <0.0L && (xcount)<2*xbranch) {
	    (xcount)++;
	    braca=sol[xcount-1][0]; bracb=braca+bracd;
	    mnbrac(&braca, &bracb, &bracc, rootbyfunc, sln, qq[qpi], pp[qpi]);
	    sol[xcount][1]=-brent(braca, bracb, bracc, rootbyfunc, sln, qq[qpi], pp[qpi], 
				  &sol[xcount][0]);
	    
	    (xcount)++;
	    braca=sol[xcount-1][0]; bracb=braca+bracd;
	    if ((checker=rootbyfunc(bracb, sln, qq[qpi], pp[qpi]))<0.0L) {
		mnbrac(&braca, &bracb, &bracc, rootbynegfunc, sln, qq[qpi], pp[qpi]);
		sol[xcount][1]=brent(braca, bracb, bracc, rootbynegfunc, sln, qq[qpi], pp[qpi], 
				     &sol[xcount][0]);
		checker=rootbyfunc(sol[xcount][0]+fabs(bracd*(braca-bracb))/(bracb-braca), 
				   sln, qq[qpi], pp[qpi]);
	    }
	    if (checker>0.0L) (xcount)--;
	}    
	
	if ((xcount)<2*xbranch) {
	    (xcount)++;
	    braca=sol[0][0]; bracb=braca-bracd;
	    mnbrac(&braca, &bracb, &bracc, rootbynegfunc, sln, qq[qpi], pp[qpi]);
	    sol[xcount][1]=brent(braca, bracb, bracc, rootbynegfunc, sln, qq[qpi], pp[qpi], 
				 &sol[xcount][0]);
	    checker=rootbyfunc(sol[xcount][0]+fabs(bracd*(braca-bracb))/(bracb-braca), 
			       sln, qq[qpi], pp[qpi]);
	    if (checker>0.0L) (xcount)--;
	    
	    while (checker <0.0L && (xcount)<2*xbranch) {
		(xcount)++;
		braca=sol[xcount-1][0]; bracb=braca-bracd;
		mnbrac(&braca, &bracb, &bracc, rootbyfunc, sln, qq[qpi], pp[qpi]);
		sol[xcount][1]=-brent(braca, bracb, bracc, rootbyfunc, sln, qq[qpi], pp[qpi], 
				      &sol[xcount][0]);
		
		(xcount)++;
		braca=sol[xcount-1][0]; bracb=braca-bracd;
		if ((checker=rootbyfunc(bracb, sln, qq[qpi], pp[qpi]))<0.0L) {
		    mnbrac(&braca, &bracb, &bracc, rootbynegfunc, sln, qq[qpi], pp[qpi]);
		    sol[xcount][1]=brent(braca, bracb, bracc, rootbynegfunc, sln, qq[qpi], pp[qpi], 
					 &sol[xcount][0]);
		    checker=rootbyfunc(sol[xcount][0]+fabs(bracd*(braca-bracb))/(bracb-braca), 
				       sln, qq[qpi], pp[qpi]);
		}
		if (checker>0.0L) (xcount)--;
	    }    
	}
	
	qsort(sol, xcount+1, 2*sizeof(long double), qsortcompare);
	j=0; for (i=1; i<=xcount; i++) {
	    if (sol[j][0]-sol[i][0]>1.25L*sqrtl(LDBL_EPSILON)) {
		sol[++j][0]=sol[i][0];
		sol[j][1]=sol[i][1];
	    }	  
	} 
	if (xcount!=j) printf("Warning, found %d times same roots twice at a=%Lf.\n", xcount-j, a0); 
	xcount=j;
	
	if (outmode>0) for (i=0; i<=xcount; i++)
	    fprintf(fl1, "%25.18Lf%25.18Lf%25.18Lf%10d\n", a0, sol[i][1], sol[i][0], xcount+1);  

	if (xcount>maxxcount) {
	    if (a0tmpstep/REFINE<a0res) {
		maxxcount=xcount;
		fprintf(fl0, "%10u%10u%22.18Lf%22.18Lf%6d%22.18Lf\n", 
			qq[qpi], pp[qpi], a0-a0tmpstep, savsolb, xcount+1, savsoly);
		fflush(fl0);
	    } 
	    else { a0-=2.0L*a0tmpstep; a0tmpstep/=REFINE; savsolb=ssavsolb; savsoly=ssavsoly;}
	}
	else { ssavsolb=savsolb; ssavsoly=savsoly; savsolb=sol[0][1]; savsoly=sol[0][0]; }
	a0+=a0tmpstep; 
    }
}


long double rootbyfunc(long double y, int sss, unsigned long int qqq, unsigned long int ppp)
{
    int brakno=2;
    zbrak(*fiterfunc, y, sss, qqq, ppp, b0, bxx, bres, &brakl[0], &braku[0], &brakno);
    if (brakno>0) return -zbrent(*fiterfunc, y, sss, qqq, ppp, brakl[brakno], braku[brakno]);
    else return fabsl(sqrtl(1.0L-omega/a0)-fabsl(y));
}

long double rootbynegfunc(long double y, int sss, unsigned long int qqq, unsigned long int ppp)
{
    int brakno=2;
    zbrak(*fiterfunc, y, sss, qqq, ppp, b0, bxx, bres, &brakl[0], &braku[0], &brakno);
    if (brakno>0) return zbrent(*fiterfunc, y, sss, qqq, ppp, brakl[brakno], braku[brakno]);
    else return 1.0L+fabsl(sqrtl(1.0L-omega/a0)-fabsl(y));
}


/* find x0(y0) so that point lies on symmetry line */
inline long double findx0(long double aa0, long double yy0, int symmln)
{
  long double tmp;

  switch (symmln) {
    case 1: {tmp=0.0L; break;}
    case 2: {tmp=0.5L; break;}
    case 3: {tmp=aa0/2.0L*(1.0L-yy0*yy0); break;}
    case 4: {tmp=aa0/2.0L*(1.0L-yy0*yy0)+0.5L; break;}
    default: {tmp=0.0L; printf("Uhoh - sln...\n");}
  }

  /* x modulo 1
  while (tmp> 0.5L-SMALL) tmp=tmp-1.0L;
  while (tmp< -0.5-SMALL) tmp=tmp+1.0L; */

  return tmp;
}

long double fiter2funcoe0(long double b, long double y, 
			  int sss, unsigned long int qqq, unsigned long int ppp) 
{
  register unsigned long int i, j, jj, qf, pf;
  register long double xxtmp, yytmp;

  pf=ppp; if (sss==1) qf=qqq-1; else qf=qqq+1;
  xxtmp=findx0(a0, y, sss); yytmp=y; jj=itstime2mod; j=jj*pf/qf;
  for (i=1;i<=ppp/2;i++) {
    ntmapstep(xxtmp, yytmp, a0, b); 
    if (i>=j){ xxtmp-=itstime2mod; jj+=itstime2mod; j=jj*pf/qf; } /* mod enat to have subtracted (q+-1)/2 in the end */ 
  }
  xxtmp-=((long double)(qf/2-(jj-itstime2mod)));
  if (sss==1) {smoothen(xxtmp-0.5L,xxtmp);}
  else {smoothen(xxtmp,xxtmp);}
  return xxtmp;
}

long double fiter2funcoe1(long double b, long double y, 
			  int sss, unsigned long int qqq, unsigned long int ppp) 
{
  register unsigned long int i, j, jj, qf, pf;
  register long double xxtmp, yytmp;

  pf=ppp; if (sss==3) qf=qqq-1; else qf=qqq+1;
  xxtmp=findx0(a0, y, sss); yytmp=y; jj=itstime2mod; j=jj*pf/qf;
  for (i=1;i<=ppp/2;i++) {
    ntmapstep(xxtmp, yytmp, a0, b); 
    if (i>=j){ xxtmp-=itstime2mod; jj+=itstime2mod; j=jj*pf/qf; } /* mod enat to have subtracted (q+-1)/2 in the end */ 
  }
  xxtmp-=((long double)(qf/2-(jj-itstime2mod)));
  if (sss==3) {smoothen(xxtmp-a0/2.0L*(1.0L-yytmp*yytmp)-0.5L,xxtmp);}
  else {smoothen(xxtmp-a0/2.0L*(1.0L-yytmp*yytmp),xxtmp);}
  return xxtmp;
}

long double fiter2funcoo0(long double b, long double y, 
			  int sss, unsigned long int qqq, unsigned long int ppp) 
{
  register unsigned long int i, j, jj, qf, pf;
  register long double xxtmp, yytmp;

  pf=ppp+1; if (sss==1) qf=qqq-1; else qf=qqq+1;
  xxtmp=findx0(a0, y, sss); yytmp=y; jj=itstime2mod; j=jj*pf/qf;
  for (i=1;i<=(ppp+1)/2;i++) {
    ntmapstep(xxtmp, yytmp, a0, b); 
    if (i>=j){ xxtmp-=itstime2mod; jj+=itstime2mod; j=jj*pf/qf; } /* mod enat to have subtracted (q+-1)/2 in the end */ 
  }
  xxtmp-=((long double)(qf/2-(jj-itstime2mod)));
  if (sss==1) {smoothen(xxtmp-a0/2.0L*(1.0L-yytmp*yytmp)-0.5L,xxtmp);}
  else {smoothen(xxtmp-a0/2.0L*(1.0L-yytmp*yytmp),xxtmp);}
  return xxtmp;
}

long double fiter2funcoo1(long double b, long double y, 
			  int sss, unsigned long int qqq, unsigned long int ppp) 
{
  register unsigned long int i, j, jj, qf, pf;
  register long double xxtmp, yytmp;

  pf=ppp-1; if (sss==3) qf=qqq-1; else qf=qqq+1;
  xxtmp=findx0(a0, y, sss); yytmp=y; jj=itstime2mod; j=jj*pf/qf;
  for (i=1;i<=(ppp-1)/2;i++) {
    ntmapstep(xxtmp, yytmp, a0, b); 
    if (i>=j){ xxtmp-=itstime2mod; jj+=itstime2mod; j=jj*pf/qf; } /* mod enat to have subtracted (q+-1)/2 in the end */ 
  }
  xxtmp-=((long double)(qf/2-(jj-itstime2mod)));
  if (sss==3) {smoothen(xxtmp-0.5L,xxtmp);}
  else {smoothen(xxtmp, xxtmp);}
  return xxtmp;
}

long double fiter2funceo0(long double b, long double y, 
			  int sss, unsigned long int qqq, unsigned long int ppp) 
{
  register unsigned long int i, j, jj, qf, pf;
  register long double xxtmp, yytmp;

  pf=ppp+1; qf=qqq;
  xxtmp=findx0(a0, y, sss); yytmp=y; jj=itstime2mod; j=jj*pf/qf;
  for (i=1;i<=(ppp+1)/2;i++) {
    ntmapstep(xxtmp, yytmp, a0, b); 
    if (i>=j){ xxtmp-=itstime2mod; jj+=itstime2mod; j=jj*pf/qf; } /* mod enat to have subtracted q/2 in the end */ 
  }
  xxtmp-=((long double)(qf/2-(jj-itstime2mod)));
  if (sss==1) {smoothen(xxtmp-a0/2.0L*(1.0L-yytmp*yytmp), xxtmp);}
  else {smoothen(xxtmp-a0/2.0L*(1.0L-yytmp*yytmp)-0.5L, xxtmp);}
  return xxtmp;
}

long double fiter2funceo1(long double b, long double y, 
			  int sss, unsigned long int qqq, unsigned long int ppp) 
{
  register unsigned long int i, j, jj, qf, pf;
  register long double xxtmp, yytmp;

  pf=ppp-1; qf=qqq;
  xxtmp=findx0(a0, y, sss); yytmp=y; jj=itstime2mod; j=jj*pf/qf;
  for (i=1;i<=(ppp-1)/2;i++) {
    ntmapstep(xxtmp, yytmp, a0, b); 
    if (i>=j){ xxtmp-=itstime2mod; jj+=itstime2mod; j=jj*pf/qf; } /* mod enat to have subtracted q/2 in the end */ 
  }
  xxtmp-=((long double)(qf/2-(jj-itstime2mod)));
  if (sss==3) {smoothen(xxtmp,xxtmp);}
  else {smoothen(xxtmp-0.5L, xxtmp);}
  return xxtmp;
}

void mkstr(char *cnummer, unsigned long int nummer)
{
    unsigned long int tmp=nummer;
    char tmpstr[160], nostr[160];

    strcpy(tmpstr,"");
    while (tmp>0) {
	switch (tmp-10*(tmp/10)) {
	    case 1: {strcpy(nostr,"1"); break;}
	    case 2: {strcpy(nostr,"2"); break;}
	    case 3: {strcpy(nostr,"3"); break;}
	    case 4: {strcpy(nostr,"4"); break;}
	    case 5: {strcpy(nostr,"5"); break;}
	    case 6: {strcpy(nostr,"6"); break;}
	    case 7: {strcpy(nostr,"7"); break;}
	    case 8: {strcpy(nostr,"8"); break;}
	    case 9: {strcpy(nostr,"9"); break;}
	    case 0: {strcpy(nostr,"0"); break;}
	}
	strcat(nostr,tmpstr);
	strcpy(tmpstr,nostr);
      	tmp/=10;
    }
    strcpy(cnummer,tmpstr);
}

void zbrak(long double (*fx)(long double, long double, int, unsigned long int, unsigned long int), 
	   long double par0, int par1, unsigned long int par2, unsigned long int par3, 
	   long double x1, long double x2, int nn, 
	   long double xb1[], long double xb2[], int *nb)
{
  int nbb, i;
  long double x, fp, fc, dx;
  
  nbb=0;
  dx=(x2-x1)/nn;
  fp=(*fx)(x=x1, par0, par1, par2, par3);
  for (i=1;i<=nn;i++) {
    fc=(*fx)(x += dx, par0, par1, par2, par3);
    if (fc*fp <= 0.0) {
      xb1[++nbb]=x-dx;
      xb2[nbb]=x;
      if (*nb==nbb) return;
    }
    fp=fc;
  }
  *nb = nbb;
}

long double zbrent(long double (*func)(long double, long double, 
				       int, unsigned long int, unsigned long int), 
		   long double par0, int par1, unsigned long int par2, unsigned long int par3, 
		   long double x1, long double x2)
{
    long double ERROR=123456789.0L, EPS=LDBL_EPSILON, ztol=LDBL_EPSILON;
    int ITMAX=100;

    int i;
    long double a=x1,b=x2,c=x2,d,e,min1,min2;
    long double fa, fb, fc,p,q,r,s,tol1,xm;

    fa = (*func)(a, par0, par1, par2, par3); fb = (*func)(b, par0, par1, par2, par3);
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
	fb=(*func)(b, par0, par1, par2, par3);
    }
    return ERROR;
}


void mnbrac(long double *ax, long double *bx, long double *cx, 
	    long double (*func)(long double, int, unsigned long int, unsigned long int), 
		   int par1, unsigned long int par2, unsigned long int par3)
{
  long double GLIMIT=10.0L, EPS=LDBL_EPSILON, TINY=EPS;
  long double GOLD=1.6180339887498948482045868343656381177203091798058L;
  
  long double ulim, u, r, q, fu, dum;
  long double fa, fb, fc;

  fa = (*func)(*ax, par1, par2, par3); fb = (*func)(*bx, par1, par2, par3);
  if (fb > fa) {
    SHFT(dum, *ax, *bx, dum);
    SHFT(dum, fb, fa, dum);
  }
  *cx = (*bx) + GOLD*(*bx-*ax);
  fc = (*func)(*cx, par1, par2, par3);
  while (fb > fc) {
    r = (*bx - *ax)*(fb-fc);
    q = (*bx - *cx)*(fb-fa);
    u = (*bx) - 
      ((*bx - *cx)*q - (*bx - *ax)*r)/(2.0*SIGN(FMAX(fabsl(q-r), TINY),q-r));
    ulim = (*bx) + GLIMIT*(*cx-*bx);
    if ((*bx-u)*(u-*cx) > 0.0) {
      fu = (*func)(u, par1, par2, par3);
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
      fu = (*func)(u, par1, par2, par3);
    } else if ((*cx-u)*(u-ulim) > 0.0) {
      fu = (*func)(u, par1, par2, par3);
      if (fu < fc) {
	SHFT(*bx, *cx, u, *cx+GOLD*(*cx-*bx));
	SHFT(fb, fc, fu, (*func)(u, par1, par2, par3));
      }
    } else if ((u-ulim)*(ulim-*cx) >= 0.0) {
      u = ulim;
      fu = (*func)(u, par1, par2, par3);
    } else {
      u = (*cx) + GOLD*(*cx-*bx);
      fu = (*func)(u, par1, par2, par3);
    }
    SHFT(*ax, *bx, *cx, u);
    SHFT(fa, fb, fc, fu);
  }
}

long double brent(long double ax, long double bx, long double cx, 
	    long double (*f)(long double, int, unsigned long int, unsigned long int), 
		   int par1, unsigned long int par2, unsigned long int par3, long double *xmin)
{
  long double ITMAX=100, ZEPS=LDBL_EPSILON, tol=sqrtl(ZEPS), ERROR=123456789.0L;
  long double CGOLD=0.3819660112501051517954131656343618822796908201942L;
  int iter;
  long double a,b,d,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm,e=0.0;

  a=(ax < cx ? ax : cx);
  b=(ax > cx ? ax : cx);
  x=w=v=bx;
  fw=fv=fx=(*f)(x, par1, par2, par3); 
  for (iter=1; iter<=ITMAX; iter++) {
    xm=0.5*(a+b);
    tol2= 2.0*(tol1=tol*fabsl(x)+ZEPS);
    if (fabsl(x-xm) <= (tol2-0.5*(b-a))) {
      *xmin = x;
      return fx;
    }
    if (fabsl(e) > tol1) {
      r = (x-w)*(fx-fv);
      q = (x-v)*(fx-fw);
      p = (x-v)*q - (x-w)*r;
      q = 2.0*(q-r);
      if (q > 0.0) p = -p;
      q = fabsl(q);
      etemp=e;
      e=d;
      if (fabsl(p) >= fabsl(0.5*q*etemp) || p <= q*(a-x) || p >= q*(b-x))
	d = CGOLD*(e=(x >= xm ? a-x : b-x));
      else {
	d = p/q;
	u = x+d;
	if (u-a < tol2 || b-u < tol2)
	  d = SIGN(tol1,xm-x);
      }
    } else {
      d = CGOLD*(e=(x >= xm ? a-x : b-x));
    }
    u = (fabsl(d) >= tol1 ? x+d : x+SIGN(tol1,d));
    fu = (*f)(u, par1, par2, par3); 
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
  return ERROR;
}

int qsortcompare(long double a[2], long double b[2])
{
  return (a[0]<b[0]? 1 : (a[0]>b[0] ? -1 : 0));
}
