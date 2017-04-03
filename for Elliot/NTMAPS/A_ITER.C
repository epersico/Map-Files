#include <stdio.h>
#include <math.h>
#include <string.h>

#define NL(flnm,dumch) {do dumch=fgetc(flnm); while (dumch!='\n' && dumch!=EOF);}
#define ntmapstep(x,y,a,b) {y-=(b)*sinl(pi2*(x)); x+=(a)-(a)*(y)*(y);}
#define smoothen(x,y) {if (x>0.0L) y=((x)/(1.0L+(x))); else y=((x)/(1.0L-(x)));}

const char *infilename="a_iter.in";
const long double SMALL=1.0e-10L;
const int SUPPERROR=15;

void makeinitials(void);
void prepoutfiles(void);
long double findx0(long double, long double, int);
long double fiter1funceven0(long double, int, unsigned long int, unsigned long int);
long double fiter1funceven1(long double, int, unsigned long int, unsigned long int);
long double fiter1funcodd0(long double, int, unsigned long int, unsigned long int);
long double fiter1funcodd1(long double, int, unsigned long int, unsigned long int);
long double fiter2funcoe0(long double, int, unsigned long int, unsigned long int);
long double fiter2funcoo0(long double, int, unsigned long int, unsigned long int);
long double fiter2funceo0(long double, int, unsigned long int, unsigned long int);
long double fiter2funcoe1(long double, int, unsigned long int, unsigned long int);
long double fiter2funcoo1(long double, int, unsigned long int, unsigned long int);
long double fiter2funceo1(long double, int, unsigned long int, unsigned long int);
void mkstr(char *, unsigned long int);

long double pi2;
long double omega, y00, b0, b00, a00, a0, axx, a, ares;
int sln, slnmin, slnmax, qpi, qp_no;
unsigned long int qq[500], pp[500], itstime2mod;
char iterfilename[160], iterfilenamebase[80];

FILE *fl1;

main(int argc, char **argv)
{
    
  long double (*fiter1func)(long double, int, unsigned long int, unsigned long int), 
    (*fiter2func)(long double, int, unsigned long int, unsigned long int);
  long double tmp, tmp2;

  /* read input */
  makeinitials();
  
  for (qpi=0; qpi<qp_no; qpi++) {
  for (sln=slnmin;sln<=slnmax;sln++) {
      
      /* prepare output files */
      prepoutfiles();
      
      /* function selection for all possible q/p-combos */
      if (2*(pp[qpi]/2)==pp[qpi]) {        /* p even (-> q odd) */
	  if (sln<=2) {
	      fiter1func=&fiter1funceven0;          /* I0 -> I0 after p/2 */
	      fiter2func=&fiter2funcoe0;            /* I_0^i -> I_0^(1-i) after p/2, (q+-1)/2 */
	  }      
	  else {
	      fiter1func=&fiter1funceven1;          /* I1 -> I1 after p/2 */
	      fiter2func=&fiter2funcoe1;            /* I_1^i -> I_1^(i-1) after p/2, (q+-1)/2 */
	  }
      }
      else {                               /* p odd */
	  if (sln<=2) {
	      fiter1func=&fiter1funcodd0;           /* I0 -> I1 after (p+1)/2 */
	      if (2*(qq[qpi]/2)==qq[qpi])      /* q even */
		  fiter2func=&fiter2funceo0;        /* I_0^i -> I_1^i after (p+1)/2, q/2 */
	      else fiter2func=&fiter2funcoo0;       /* I_0^i -> I_1^(1-i) after (p+1)/2, (q+-1)/2 */
	  }
	  else {
	      fiter1func=&fiter1funcodd1;           /* I1 -> I0 after (p-1)/2 */
	      if (2*(qq[qpi]/2)==qq[qpi])      /* q even */
		  fiter2func=&fiter2funceo1;        /* I_0^1 -> I_1^1 after (p+1)/2, q/2 */
	      else fiter2func=&fiter2funcoo1;       /* I_0^1 -> I_1^0 after (p+1)/2, (q+1)/2 */
	  }
      }
      
      
      for (a=a0+ares; a<=axx; a+=ares) {
	  tmp=fiter2func(y00,sln,qq[qpi],pp[qpi]); smoothen(tmp,tmp2);
	  fprintf(fl1,"%20.13Lf%25.15Lg%25.15Lg%25.15Lg\n", 
	     a, fiter1func(y00,sln,qq[qpi],pp[qpi]), tanhl(tmp), tmp2);
      }
      
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
  fscanf(fl,"%Lf", &y00);  NL(fl, dummy);               /* y00 */
  fscanf(fl,"%Lf", &b00);  NL(fl, dummy);              /* b0 */
  fscanf(fl,"%Lf", &a00);  NL(fl, dummy);              /* start a */
  fscanf(fl,"%Lf", &axx);  NL(fl, dummy);              /* max a */
  fscanf(fl,"%Lg", &ares);  NL(fl, dummy);             /* a resolution */
  fscanf(fl,"%u", &itstime2mod);  NL(fl, dummy);       /* mod intervals threshold */
  fscanf(fl,"%s", iterfilenamebase); NL(fl, dummy);    /* outfilename */
  fscanf(fl,"%d", &qp_no);  NL(fl, dummy);    /* number of p/q values below */
  for (i=0;i<qp_no;i++) {
      fscanf(fl,"%u %u", &qq[i], &pp[i]); NL(fl, dummy);                     /* q/p */
  }

  if (sln==0) {slnmin=1; slnmax=4;}
  else {slnmin=sln; slnmax=sln;}
  if (axx<a00) {tmp=a00; a00=axx; axx=tmp;}
  if (ares<SMALL) ares=axx-a00+SMALL;
  fclose(fl);
}


void prepoutfiles(void)
{
  char ctmp[60];
    
  /* adjust initials */
  omega=((long double)(qq[qpi]))/((long double)(pp[qpi]));
  if (a00<=omega) {
    a0=omega;
    b0=0.0L;
  }
  else {a0=a00; b0=b00;}
  if (b0<=0.0L) b0=0.0L; else b0=b00;

  strcpy(iterfilename,iterfilenamebase);
  strcpy(ctmp,""); mkstr(&ctmp[0],qq[qpi]); strcat(iterfilename,ctmp);
  strcat(iterfilename,"_");
  strcpy(ctmp,""); mkstr(&ctmp[0],pp[qpi]); strcat(iterfilename,ctmp);
  strcat(iterfilename,".s");
  strcpy(ctmp,""); mkstr(&ctmp[0],sln); strcat(iterfilename,ctmp);

  fl1=fopen(iterfilename, "w");
  fprintf(fl1,"# Symm.l.: %2d \n", sln);
  fprintf(fl1,"# q:       %10u , p:       %10u \n", qq[qpi], pp[qpi]);
  fprintf(fl1,"# y00:      %10.7Lf , b0:      %10.7Lf \n", y00, b0); 
  fprintf(fl1,"# a0:      %10.7Lf , ax:      %10.7Lf \n", a0, axx); 
  fprintf(fl1,"# ares:    %10.7Lf , mod:     %10u \n", ares, itstime2mod); 
  fprintf(fl1,"# %18s%25s%25s%25s\n\n", "y", "f_iter1", "f_iter2", "f_iter2b");

}

/* find x0(y00) so that point lies on symmetry line */
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

long double fiter1funceven0(long double yy0, int sss, unsigned long int qqq, unsigned long int ppp) 
{
  int i;
  register long double xxtmp, yytmp;
  xxtmp=findx0(a, yy0, sss); yytmp=yy0;
  for (i=1;i<=ppp/2;i++) {
    ntmapstep(xxtmp, yytmp, a, b0);
    if (fabsl(xxtmp)>itstime2mod) 
	xxtmp=dreml(xxtmp,1.0L); /* mod every now and then to keep x handy) */ 
  }  
  return sinl(pi2*xxtmp);
}

long double fiter1funceven1(long double yy0, int sss, unsigned long int qqq, unsigned long int ppp) 
{
  int i;
  register long double xxtmp, yytmp;
  xxtmp=findx0(a, yy0, sss); yytmp=yy0;
  for (i=1;i<=ppp/2;i++) {
    ntmapstep(xxtmp, yytmp, a, b0);
    if (fabsl(xxtmp)>itstime2mod) 
	xxtmp=dreml(xxtmp,1.0L); /* mod every now and then to keep x handy) */ 
  }  
  return sinl(pi2*(xxtmp-a/2.0L*(1.0L-yytmp*yytmp)));
}

long double fiter1funcodd0(long double yy0, int sss, unsigned long int qqq, unsigned long int ppp) 
{
  int i;
  register long double xxtmp, yytmp;
  xxtmp=findx0(a, yy0, sss); yytmp=yy0;
  for (i=1;i<=ppp/2+1;i++) {
    ntmapstep(xxtmp, yytmp, a, b0);
    if (fabsl(xxtmp)>itstime2mod) 
	xxtmp=dreml(xxtmp,1.0L); /* mod every now and then to keep x handy) */ 
  }  
  return sinl(pi2*(xxtmp-a/2.0L*(1.0L-yytmp*yytmp)));
}

long double fiter1funcodd1(long double yy0, int sss, unsigned long int qqq, unsigned long int ppp) 
{
  int i;
  register long double xxtmp, yytmp;
  xxtmp=findx0(a, yy0, sss); yytmp=yy0;
  for (i=1;i<=ppp/2;i++) {
    ntmapstep(xxtmp, yytmp, a, b0);
    if (fabsl(xxtmp)>itstime2mod) 
	xxtmp=dreml(xxtmp,1.0L); /* mod every now and then to keep x handy) */ 
  }  
  return sinl(pi2*xxtmp);
}

long double fiter2funcoe0(long double yy0, int sss, unsigned long int qqq, unsigned long int ppp) 
{
  register unsigned long int i, j, jj, qf, pf;
  register long double xxtmp, yytmp;

  pf=ppp; if (sss==1) qf=qqq-1; else qf=qqq+1;
  xxtmp=findx0(a, yy0, sss); yytmp=yy0; jj=itstime2mod; j=jj*pf/qf;
  for (i=1;i<=ppp/2;i++) {
    ntmapstep(xxtmp, yytmp, a, b0); 
    if (i>=j){ xxtmp-=itstime2mod; jj+=itstime2mod; j=jj*pf/qf; } /* mod enat to have subtracted (q+-1)/2 in the end */ 
  }
  xxtmp-=((long double)(qf/2-(jj-itstime2mod)));
  if (sss==1) return xxtmp-0.5L;
  else return xxtmp;
}

long double fiter2funcoe1(long double yy0, int sss, unsigned long int qqq, unsigned long int ppp) 
{
  register unsigned long int i, j, jj, qf, pf;
  register long double xxtmp, yytmp;

  pf=ppp; if (sss==3) qf=qqq-1; else qf=qqq+1;
  xxtmp=findx0(a, yy0, sss); yytmp=yy0; jj=itstime2mod; j=jj*pf/qf;
  for (i=1;i<=ppp/2;i++) {
    ntmapstep(xxtmp, yytmp, a, b0); 
    if (i>=j){ xxtmp-=itstime2mod; jj+=itstime2mod; j=jj*pf/qf; } /* mod enat to have subtracted (q+-1)/2 in the end */ 
  }
  xxtmp-=((long double)(qf/2-(jj-itstime2mod)));
  if (sss==3) return xxtmp-a/2.0L*(1.0L-yytmp*yytmp)-0.5L;
  else return xxtmp-a/2.0L*(1.0L-yytmp*yytmp);
}

long double fiter2funcoo0(long double yy0, int sss, unsigned long int qqq, unsigned long int ppp) 
{
  register unsigned long int i, j, jj, qf, pf;
  register long double xxtmp, yytmp;

  pf=ppp+1; if (sss==1) qf=qqq-1; else qf=qqq+1;
  xxtmp=findx0(a, yy0, sss); yytmp=yy0; jj=itstime2mod; j=jj*pf/qf;
  for (i=1;i<=(ppp+1)/2;i++) {
    ntmapstep(xxtmp, yytmp, a, b0); 
    if (i>=j){ xxtmp-=itstime2mod; jj+=itstime2mod; j=jj*pf/qf; } /* mod enat to have subtracted (q+-1)/2 in the end */ 
  }
  xxtmp-=((long double)(qf/2-(jj-itstime2mod)));
  if (sss==1) return xxtmp-a/2.0L*(1.0L-yytmp*yytmp)-0.5L;
  else return xxtmp-a/2.0L*(1.0L-yytmp*yytmp);
}

long double fiter2funcoo1(long double yy0, int sss, unsigned long int qqq, unsigned long int ppp) 
{
  register unsigned long int i, j, jj, qf, pf;
  register long double xxtmp, yytmp;

  pf=ppp-1; if (sss==3) qf=qqq-1; else qf=qqq+1;
  xxtmp=findx0(a, yy0, sss); yytmp=yy0; jj=itstime2mod; j=jj*pf/qf;
  for (i=1;i<=(ppp-1)/2;i++) {
    ntmapstep(xxtmp, yytmp, a, b0); 
    if (i>=j){ xxtmp-=itstime2mod; jj+=itstime2mod; j=jj*pf/qf; } /* mod enat to have subtracted (q+-1)/2 in the end */ 
  }
  xxtmp-=((long double)(qf/2-(jj-itstime2mod)));
  if (sss==3) return xxtmp-0.5L;
  else return xxtmp;
}

long double fiter2funceo0(long double yy0, int sss, unsigned long int qqq, unsigned long int ppp) 
{
  register unsigned long int i, j, jj, qf, pf;
  register long double xxtmp, yytmp;

  pf=ppp+1; qf=qqq;
  xxtmp=findx0(a, yy0, sss); yytmp=yy0; jj=itstime2mod; j=jj*pf/qf;
  for (i=1;i<=(ppp+1)/2;i++) {
    ntmapstep(xxtmp, yytmp, a, b0); 
    if (i>=j){ xxtmp-=itstime2mod; jj+=itstime2mod; j=jj*pf/qf; } /* mod enat to have subtracted q/2 in the end */ 
  }
  xxtmp-=((long double)(qf/2-(jj-itstime2mod)));
  if (sss==1) return xxtmp-a/2.0L*(1.0L-yytmp*yytmp);
  else return xxtmp-a/2.0L*(1.0L-yytmp*yytmp)-0.5L;
}

long double fiter2funceo1(long double yy0, int sss, unsigned long int qqq, unsigned long int ppp) 
{
  register unsigned long int i, j, jj, qf, pf;
  register long double xxtmp, yytmp;

  pf=ppp-1; qf=qqq;
  xxtmp=findx0(a, yy0, sss); yytmp=yy0; jj=itstime2mod; j=jj*pf/qf;
  for (i=1;i<=(ppp-1)/2;i++) {
    ntmapstep(xxtmp, yytmp, a, b0); 
    if (i>=j){ xxtmp-=itstime2mod; jj+=itstime2mod; j=jj*pf/qf; } /* mod enat to have subtracted q/2 in the end */ 
  }
  xxtmp-=((long double)(qf/2-(jj-itstime2mod)));
  if (sss==3) return xxtmp;
  else return xxtmp-0.5L;
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

