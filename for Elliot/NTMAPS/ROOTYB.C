#include <stdio.h>
#include <math.h>
#include <string.h>
#include <float.h>

#define NL(flnm,dumch) {do dumch=fgetc(flnm); while (dumch!='\n' && dumch!=EOF);}
#define ntmapstep(x,y,a,b) {y-=(b)*sinl(pi2*(x)); x+=(a)-(a)*(y)*(y);}
#define smoothen(x,y) {if (x>0.0L) y=((x)/(1.0L+(x))); else y=((x)/(1.0L-(x)));}

const char *infilename="rootyb.in";
const long double SMALL=1.0e-10L;
const int SUPPERROR=15;

void makeinitials(void);
void prepoutfiles(void);
void maincalc(void);
long double findx0(long double, long double, int);
long double fiter2funcoe0(long double, int, unsigned long int, unsigned long int);
long double fiter2funcoo0(long double, int, unsigned long int, unsigned long int);
long double fiter2funceo0(long double, int, unsigned long int, unsigned long int);
long double fiter2funcoe1(long double, int, unsigned long int, unsigned long int);
long double fiter2funcoo1(long double, int, unsigned long int, unsigned long int);
long double fiter2funceo1(long double, int, unsigned long int, unsigned long int);
void zbrak(long double (*)(long double, int, unsigned long int, unsigned long int), 
	   int, unsigned long int, unsigned long int, 
	   long double, long double, int, long double [], long double [], int *);
long double zbrent(long double (*)(long double, int, unsigned long int, unsigned long int), 
		   int, unsigned long int, unsigned long int, long double, long double);
void mkstr(char *, unsigned long int);

long double pi2;
long double omega, a0, a00, b0, b00, bmin, bmax, bstep, y00, yxx;
int sln, slnmin, slnmax, qpi, qp_no, yres, yresx;
unsigned long int qq[500], pp[500], itstime2mod;
char iterfilename[160], iterfilenamebase[80];

FILE *fl1;
long double (*fiterfunc)(long double, int, unsigned long int, unsigned long int);

main(int argc, char **argv)
{
    
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
      
      for (b0=bmin; b0<=bmax; b0+=bstep) maincalc();
  }
  }
  
  fclose(fl1);
  exit(0); 
} 


void makeinitials(void)      /* read and adjust initial values */
{
  FILE *fl;
  char dummy;
  int i;
  long double tmp;
  
  pi2=8.0L*atanl(1.0L);

  /* read initials */
  fl=fopen(infilename, "r");
  
  fscanf(fl,"%d", &sln); NL(fl, dummy);                /* shinohara point */
  fscanf(fl,"%Lf", &a00);  NL(fl, dummy);              /* a0 */
  fscanf(fl,"%Lf %Lf", &b00, &bmax);  NL(fl, dummy);  /* b-range */
  fscanf(fl,"%Lf", &bstep);  NL(fl, dummy);            /* bstep */
  fscanf(fl,"%Lf", &y00);  NL(fl, dummy);              /* min y */
  fscanf(fl,"%Lf", &yxx);  NL(fl, dummy);              /* max y */
  fscanf(fl,"%d %d", &yres, &yresx);  NL(fl, dummy);   /* # of brackets, # used */
  fscanf(fl,"%u", &itstime2mod);  NL(fl, dummy);       /* mod intervals threshold */
  fscanf(fl,"%s", iterfilenamebase); NL(fl, dummy);    /* outfilename */
  fscanf(fl,"%d", &qp_no);  NL(fl, dummy);    /* number of p/q values below */
  for (i=0;i<qp_no;i++) {
      fscanf(fl,"%u %u", &qq[i], &pp[i]); NL(fl, dummy);                     /* q/p */
  }

  if (sln==0) {slnmin=1; slnmax=4;}
  else {slnmin=sln; slnmax=sln;}
  if (yxx<y00) {tmp=y00; y00=yxx; yxx=tmp;}
  if (yres<1) yres=1;
  fclose(fl);
}


void prepoutfiles(void)
{
  char ctmp[60];
    
  /* adjust initials */
  omega=((long double)(qq[qpi]))/((long double)(pp[qpi]));
  if (a00<=omega) {
    a0=omega;
    bmin=0.0L;
  }
  else {a0=a00; bmin=b00;}
  if (bmin<=0.0L) bmin=0.0L; else bmin=b00;

  strcpy(iterfilename,iterfilenamebase);
  strcpy(ctmp,""); mkstr(&ctmp[0],qq[qpi]); strcat(iterfilename,ctmp);
  strcat(iterfilename,"_");
  strcpy(ctmp,""); mkstr(&ctmp[0],pp[qpi]); strcat(iterfilename,ctmp);
  strcat(iterfilename,".s");
  strcpy(ctmp,""); mkstr(&ctmp[0],sln); strcat(iterfilename,ctmp);

  fl1=fopen(iterfilename, "w");
  fprintf(fl1,"# Symm.l.: %2d \n", sln);
  fprintf(fl1,"# q:       %10u , p:       %10u \n", qq[qpi], pp[qpi]);
  fprintf(fl1,"# a0:      %10.7Lf , bstep:   %10.7Lf \n", a0, bstep); 
  fprintf(fl1,"# bmin:    %10.7Lf , bmax:    %10.7Lf \n", bmin, bmax); 
  fprintf(fl1,"# y0:      %10.7Lf , yx:      %10.7Lf \n", y00, yxx); 
  fprintf(fl1,"# yres:    %10d , yresx:   %10u \n", yres, yresx); 
  fprintf(fl1,"# mod:     %10u \n", itstime2mod); 
  fprintf(fl1,"# %23s%25s%15d\n\n", "b", "y", "# sol");
}

void maincalc(void) 
{
    long double braku[yres+1], brakl[yres+1];
    int brakno=(yresx==0 ? yres : yresx), i;

    zbrak(*fiterfunc, sln, qq[qpi], pp[qpi], y00, yxx, yres, &brakl[0], &braku[0], &brakno);
    if (brakno>0) fprintf(fl1, "%25.18Lf%25.18Lf%15d\n", b0,
			  zbrent(*fiterfunc, sln, qq[qpi], pp[qpi], brakl[brakno], braku[brakno]),
			  brakno);
    if (yresx==0) for (i=1; i<=brakno-1; i++) 
	fprintf(fl1, "%25.18Lf%25.18Lf\n", b0, 
		zbrent(*fiterfunc, sln, qq[qpi], pp[qpi], brakl[i], braku[i]));
    
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


long double fiter2funcoe0(long double yy0, int sss, unsigned long int qqq, unsigned long int ppp) 
{
  register unsigned long int i, j, jj, qf, pf;
  register long double xxtmp, yytmp;

  pf=ppp; if (sss==1) qf=qqq-1; else qf=qqq+1;
  xxtmp=findx0(a0, yy0, sss); yytmp=yy0; jj=itstime2mod; j=(unsigned long int)((double)(jj)*(((double)(pf))/((double)(qf))));
  for (i=1;i<=ppp/2;i++) {
    ntmapstep(xxtmp, yytmp, a0, b0); 
    if (i>=j){ xxtmp-=itstime2mod; jj+=itstime2mod; j=(unsigned long int)((double)(jj)*(((double)(pf))/((double)(qf)))); } /* mod enat to have subtracted (q+-1)/2 in the end */ 
  }
  xxtmp-=((long double)(qf/2-(jj-itstime2mod)));
  if (sss==1) {smoothen(xxtmp-0.5L,xxtmp);}
  else {smoothen(xxtmp,xxtmp);}
  return xxtmp;
}

long double fiter2funcoe1(long double yy0, int sss, unsigned long int qqq, unsigned long int ppp) 
{
  register unsigned long int i, j, jj, qf, pf;
  register long double xxtmp, yytmp;

  pf=ppp; if (sss==3) qf=qqq-1; else qf=qqq+1;
  xxtmp=findx0(a0, yy0, sss); yytmp=yy0; jj=itstime2mod; j=(unsigned long int)((double)(jj)*(((double)(pf))/((double)(qf))));
  for (i=1;i<=ppp/2;i++) {
    ntmapstep(xxtmp, yytmp, a0, b0); 
    if (i>=j){ xxtmp-=itstime2mod; jj+=itstime2mod; j=(unsigned long int)((double)(jj)*(((double)(pf))/((double)(qf)))); } /* mod enat to have subtracted (q+-1)/2 in the end */ 
  }
  xxtmp-=((long double)(qf/2-(jj-itstime2mod)));
  if (sss==3) {smoothen(xxtmp-a0/2.0L*(1.0L-yytmp*yytmp)-0.5L,xxtmp);}
  else {smoothen(xxtmp-a0/2.0L*(1.0L-yytmp*yytmp),xxtmp);}
  return xxtmp;
}

long double fiter2funcoo0(long double yy0, int sss, unsigned long int qqq, unsigned long int ppp) 
{
  register unsigned long int i, j, jj, qf, pf;
  register long double xxtmp, yytmp;

  pf=ppp+1; if (sss==1) qf=qqq-1; else qf=qqq+1;
  xxtmp=findx0(a0, yy0, sss); yytmp=yy0; jj=itstime2mod; j=(unsigned long int)((double)(jj)*(((double)(pf))/((double)(qf))));
  for (i=1;i<=(ppp+1)/2;i++) {
    ntmapstep(xxtmp, yytmp, a0, b0); 
    if (i>=j){ xxtmp-=itstime2mod; jj+=itstime2mod; j=(unsigned long int)((double)(jj)*(((double)(pf))/((double)(qf)))); } /* mod enat to have subtracted (q+-1)/2 in the end */ 
  }
  xxtmp-=((long double)(qf/2-(jj-itstime2mod)));
  if (sss==1) {smoothen(xxtmp-a0/2.0L*(1.0L-yytmp*yytmp)-0.5L,xxtmp);}
  else {smoothen(xxtmp-a0/2.0L*(1.0L-yytmp*yytmp),xxtmp);}
  return xxtmp;
}

long double fiter2funcoo1(long double yy0, int sss, unsigned long int qqq, unsigned long int ppp) 
{
  register unsigned long int i, j, jj, qf, pf;
  register long double xxtmp, yytmp;

  pf=ppp-1; if (sss==3) qf=qqq-1; else qf=qqq+1;
  xxtmp=findx0(a0, yy0, sss); yytmp=yy0; jj=itstime2mod; j=(unsigned long int)((double)(jj)*(((double)(pf))/((double)(qf))));
  for (i=1;i<=(ppp-1)/2;i++) {
    ntmapstep(xxtmp, yytmp, a0, b0); 
    if (i>=j){ xxtmp-=itstime2mod; jj+=itstime2mod; j=(unsigned long int)((double)(jj)*(((double)(pf))/((double)(qf)))); } /* mod enat to have subtracted (q+-1)/2 in the end */ 
  }
  xxtmp-=((long double)(qf/2-(jj-itstime2mod)));
  if (sss==3) {smoothen(xxtmp-0.5L,xxtmp);}
  else {smoothen(xxtmp, xxtmp);}
  return xxtmp;
}

long double fiter2funceo0(long double yy0, int sss, unsigned long int qqq, unsigned long int ppp) 
{
  register unsigned long int i, j, jj, qf, pf;
  register long double xxtmp, yytmp;

  pf=ppp+1; qf=qqq;
  xxtmp=findx0(a0, yy0, sss); yytmp=yy0; jj=itstime2mod; j=(unsigned long int)((double)(jj)*(((double)(pf))/((double)(qf))));
  for (i=1;i<=(ppp+1)/2;i++) {
    ntmapstep(xxtmp, yytmp, a0, b0); 
    if (i>=j){ xxtmp-=itstime2mod; jj+=itstime2mod; j=(unsigned long int)((double)(jj)*(((double)(pf))/((double)(qf)))); } /* mod enat to have subtracted q/2 in the end */ 
  }
  xxtmp-=((long double)(qf/2-(jj-itstime2mod)));
  if (sss==1) {smoothen(xxtmp-a0/2.0L*(1.0L-yytmp*yytmp), xxtmp);}
  else {smoothen(xxtmp-a0/2.0L*(1.0L-yytmp*yytmp)-0.5L, xxtmp);}
  return xxtmp;
}

long double fiter2funceo1(long double yy0, int sss, unsigned long int qqq, unsigned long int ppp) 
{
  register unsigned long int i, j, jj, qf, pf;
  register long double xxtmp, yytmp;

  pf=ppp-1; qf=qqq;
  xxtmp=findx0(a0, yy0, sss); yytmp=yy0; jj=itstime2mod; j=(unsigned long int)((double)(jj)*(((double)(pf))/((double)(qf))));
  for (i=1;i<=(ppp-1)/2;i++) {
    ntmapstep(xxtmp, yytmp, a0, b0); 
    if (i>=j){ xxtmp-=itstime2mod; jj+=itstime2mod; j=(unsigned long int)((double)(jj)*(((double)(pf))/((double)(qf)))); } /* mod enat to have subtracted q/2 in the end */ 
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
void zbrak(long double (*fx)(long double, int, unsigned long int, unsigned long int), 
	   int par1, unsigned long int par2, unsigned long int par3, 
	   long double x1, long double x2, int nn, 
	   long double xb1[], long double xb2[], int *nb)
{
  int nbb, i;
  long double x, fp, fc, dx;
  
  nbb=0;
  dx=(x2-x1)/nn;
  fp=(*fx)(x=x1, par1, par2, par3);
  for (i=1;i<=nn;i++) {
    fc=(*fx)(x += dx, par1, par2, par3);
    if (fc*fp <= 0.0) {
      xb1[++nbb]=x-dx;
      xb2[nbb]=x;
      if (*nb==nbb) return;
    }
    fp=fc;
  }
  *nb = nbb;
}

long double zbrent(long double (*func)(long double, int, unsigned long int, unsigned long int), 
		   int par1, unsigned long int par2, unsigned long int par3, 
		   long double x1, long double x2)
{
    long double ERROR=123456789.0L, EPS=LDBL_EPSILON, ztol=LDBL_EPSILON;
    int ITMAX=100;

    int i;
    long double a=x1,b=x2,c=x2,d,e,min1,min2;
    long double fa, fb, fc,p,q,r,s,tol1,xm;

    fa = (*func)(a, par1, par2, par3); fb = (*func)(b, par1, par2, par3);
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
	fb=(*func)(b, par1, par2, par3);
    }
    return ERROR;
}

