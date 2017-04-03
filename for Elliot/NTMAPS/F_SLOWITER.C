#include <stdio.h>
#include <math.h>
#include <string.h>
#define NL(flnm,dumch) {do dumch=fgetc(flnm); while (dumch!='\n' && dumch!=EOF);}
#define itstime2mod 1000.0L

const char *infilename="f_iter.in";
const long double SMALL=1.0e-10L;
const int SUPPERROR=15;

void makeinitials(void);
void prepoutfiles(void);
void ntmapstep(long double *, long double *, long double);
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

long double pi;
long double omega, a0, a00, b0, b00, y00, yxx, yres, smooththresh;
int sln, slnmin, slnmax, qpi, qp_no;
unsigned long int qq[500], pp[500];
char iterfilename[160], iterfilenamebase[80];

FILE *fl1;

main(int argc, char **argv)
{
    
  long double (*fiter1func)(long double, int, unsigned long int, unsigned long int), 
    (*fiter2func)(long double, int, unsigned long int, unsigned long int);
  long double smoothtmp, tmp2, y;

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
      
      
      tmp2=fiter2func(y00,sln,qq[qpi],pp[qpi]);
      if (fabsl(tmp2)<smooththresh) smoothtmp=tmp2;
      else smoothtmp=(tmp2/fabsl(tmp2))*smooththresh;
      for (y=y00+yres; y<=yxx; y+=yres) {
	  tmp2=fiter2func(y,sln,qq[qpi],pp[qpi]);
	  if (fabsl(tmp2)<smooththresh) smoothtmp=tmp2;
	  fprintf(fl1,"%20.13Lf%25.15Lg%25.15Lg%25.15Lg\n", 
		  y, fiter1func(y,sln,qq[qpi],pp[qpi]), tmp2, smoothtmp);
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
  
  pi=4.0L*atanl(1.0L);

  /* read initials */
  fl=fopen(infilename, "r");
  
  fscanf(fl,"%Ld", &sln); NL(fl, dummy);               /* shinohara point */
  fscanf(fl,"%Lf", &a00);  NL(fl, dummy);              /* a0 */
  fscanf(fl,"%Lf", &b00);  NL(fl, dummy);              /* b0 */
  fscanf(fl,"%Lf", &y00);  NL(fl, dummy);              /* start y */
  fscanf(fl,"%Lf", &yxx);  NL(fl, dummy);              /* max y */
  fscanf(fl,"%Lg", &yres);  NL(fl, dummy);             /* y resolution */
  fscanf(fl,"%Lf", &smooththresh);  NL(fl, dummy);     /* smoothening threshold */
  fscanf(fl,"%s", iterfilenamebase); NL(fl, dummy);   /* outfilename */
  fscanf(fl,"%d", &qp_no);  NL(fl, dummy);    /* number of p/q values below */
  for (i=0;i<qp_no;i++) {
      fscanf(fl,"%u %u", &qq[i], &pp[i]); NL(fl, dummy);                     /* q/p */
  }

  if (sln==0) {slnmin=1; slnmax=4;}
  else {slnmin=sln; slnmax=sln;}
  if (yxx<y00) {tmp=y00; y00=yxx; yxx=tmp;}
  if (yres<SMALL) yres=yxx-y00+SMALL;
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
  fprintf(fl1,"# a0:      %10.7Lf , b0:      %10.7Lf \n", a0, b0); 
  fprintf(fl1,"# y0:      %10.7Lf , yx:      %10.7Lf \n", y00, yxx); 
  fprintf(fl1,"# yres:    %10.7Lf , smooth:  %20.13Lg \n", yres, smooththresh); 
  fprintf(fl1,"# %18s%25s%25s%25s\n\n", "y", "f_iter1", "f_iter2", "f_iter2s");

}

/* single iteration of NT-map without mod */
void ntmapstep(long double *x, long double *y, long double bb)
{

  *y = *y - bb*sinl(2.0L*pi*(*x));
  *x = *x + a0*(1.0L-(*y)*(*y));
  
  /* x modulo 1 
  while (*x> 0.5001L) {*x=*x-1.0L; (*lft)++;}
  while (*x< -0.5L) {*x=*x+1.0L; (*lft)--;} */
}

/* find x0(y0) so that point lies on symmetry line */
long double findx0(long double aa0, long double yy0, int symmln)
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
  long double xxtmp, yytmp;
  xxtmp=findx0(a0, yy0, sss); yytmp=yy0;
  for (i=1;i<=ppp/2;i++) {
    ntmapstep(&xxtmp, &yytmp, b0);
    if (fabsl(xxtmp)>itstime2mod) 
	xxtmp=dreml(xxtmp,1.0L); /* mod every now and then to keep x handy) */ 
  }  
  return sinl(2.0L*pi*xxtmp);
}

long double fiter1funceven1(long double yy0, int sss, unsigned long int qqq, unsigned long int ppp) 
{
  int i;
  long double xxtmp, yytmp;
  xxtmp=findx0(a0, yy0, sss); yytmp=yy0;
  for (i=1;i<=ppp/2;i++) {
    ntmapstep(&xxtmp, &yytmp, b0);
    if (fabsl(xxtmp)>itstime2mod) 
	xxtmp=dreml(xxtmp,1.0L); /* mod every now and then to keep x handy) */ 
  }  
  return sinl(2.0L*pi*(xxtmp-a0/2.0L*(1.0L-yytmp*yytmp)));
}

long double fiter1funcodd0(long double yy0, int sss, unsigned long int qqq, unsigned long int ppp) 
{
  int i;
  long double xxtmp, yytmp;
  xxtmp=findx0(a0, yy0, sss); yytmp=yy0;
  for (i=1;i<=ppp/2+1;i++) {
    ntmapstep(&xxtmp, &yytmp, b0);
    if (fabsl(xxtmp)>itstime2mod) 
	xxtmp=dreml(xxtmp,1.0L); /* mod every now and then to keep x handy) */ 
  }  
  return sinl(2.0L*pi*(xxtmp-a0/2.0L*(1.0L-yytmp*yytmp)));
}

long double fiter1funcodd1(long double yy0, int sss, unsigned long int qqq, unsigned long int ppp) 
{
  int i;
  long double xxtmp, yytmp;
  xxtmp=findx0(a0, yy0, sss); yytmp=yy0;
  for (i=1;i<=ppp/2;i++) {
    ntmapstep(&xxtmp, &yytmp, b0);
    if (fabsl(xxtmp)>itstime2mod) 
	xxtmp=dreml(xxtmp,1.0L); /* mod every now and then to keep x handy) */ 
  }  
  return sinl(2.0L*pi*xxtmp);
}

long double fiter2funcoe0(long double yy0, int sss, unsigned long int qqq, unsigned long int ppp) 
{
  int i, j=1.0L;
  long double xxtmp, yytmp, omh;

  if (sss==1) omh=((long double)(qqq-1))/((long double)(ppp));
  else omh=((long double)(qqq+1))/((long double)(ppp));
  xxtmp=findx0(a0, yy0, sss); yytmp=yy0;
  for (i=1;i<=ppp/2;i++) {
    ntmapstep(&xxtmp, &yytmp, b0); 
    if (i>=j/omh){ xxtmp-=1.0L; j++; } /* mod enat to have subtracted (q+-1)/2 in the end */ 
  }
  if (sss==1) return tanhl(xxtmp-0.5L);
  else return tanhl(xxtmp);
}

long double fiter2funcoe1(long double yy0, int sss, unsigned long int qqq, unsigned long int ppp) 
{
  int i, j=1.0L;
  long double xxtmp, yytmp, omh;

  if (sss==3) omh=((long double)(qqq-1))/((long double)(ppp));
  else omh=((long double)(qqq+1))/((long double)(ppp));
  xxtmp=findx0(a0, yy0, sss); yytmp=yy0;
  for (i=1;i<=ppp/2;i++) {
    ntmapstep(&xxtmp, &yytmp, b0); 
    if (i>=j/omh){ xxtmp-=1.0L; j++; } /* mod enat to have subtracted (q+-1)/2 in the end */ 
  }
  if (sss==3) return tanhl(xxtmp-a0/2.0L*(1.0L-yytmp*yytmp)-0.5L);
  else return tanhl(xxtmp-a0/2.0L*(1.0L-yytmp*yytmp));
}

long double fiter2funcoo0(long double yy0, int sss, unsigned long int qqq, unsigned long int ppp) 
{
  int i, j=1.0L;
  long double xxtmp, yytmp, omh;

  if (sss==1) omh=((long double)(qqq-1))/((long double)(ppp+1));
  else omh=((long double)(qqq+1))/((long double)(ppp+1));
  xxtmp=findx0(a0, yy0, sss); yytmp=yy0;
  for (i=1;i<=(ppp+1)/2;i++) {
    ntmapstep(&xxtmp, &yytmp, b0); 
    if (i>=j/omh){ xxtmp-=1.0L; j++; } /* mod enat to have subtracted (q+-1)/2 in the end */ 
  }
  if (sss==1) return tanhl(xxtmp-a0/2.0L*(1.0L-yytmp*yytmp)-0.5L);
  else return tanhl(xxtmp-a0/2.0L*(1.0L-yytmp*yytmp));
}

long double fiter2funcoo1(long double yy0, int sss, unsigned long int qqq, unsigned long int ppp) 
{
  int i, j=1.0L;
  long double xxtmp, yytmp, omh;

  if (sss==3) omh=((long double)(qqq-1))/((long double)(ppp-1));
  else omh=((long double)(qqq+1))/((long double)(ppp-1));
  xxtmp=findx0(a0, yy0, sss); yytmp=yy0;
  for (i=1;i<=(ppp-1)/2;i++) {
    ntmapstep(&xxtmp, &yytmp, b0); 
    if (i>=j/omh){ xxtmp-=1.0L; j++; } /* mod enat to have subtracted (q+-1)/2 in the end */ 
  }
  if (sss==3) return tanhl(xxtmp-0.5L);
  else return tanhl(xxtmp);
}

long double fiter2funceo0(long double yy0, int sss, unsigned long int qqq, unsigned long int ppp) 
{
  int i, j=1.0L;
  long double xxtmp, yytmp, omh;

  omh=((long double)(qqq))/((long double)(ppp+1));
  xxtmp=findx0(a0, yy0, sss); yytmp=yy0;
  for (i=1;i<=(ppp+1)/2;i++) {
    ntmapstep(&xxtmp, &yytmp, b0); 
    if (i>=j/omh){ xxtmp-=1.0L; j++; } /* mod enat to have subtracted q/2 in the end */ 
  }
  if (sss==1) return tanhl(xxtmp-a0/2.0L*(1.0L-yytmp*yytmp));
  else return tanhl(xxtmp-a0/2.0L*(1.0L-yytmp*yytmp)-0.5L);
}

long double fiter2funceo1(long double yy0, int sss, unsigned long int qqq, unsigned long int ppp) 
{
  int i, j=1.0L;
  long double xxtmp, yytmp, omh;

  omh=((long double)(qqq))/((long double)(ppp-1));
  xxtmp=findx0(a0, yy0, sss); yytmp=yy0;
  for (i=1;i<=(ppp-1)/2;i++) {
    ntmapstep(&xxtmp, &yytmp, b0); 
    if (i>=j/omh){ xxtmp-=1.0L; j++; } /* mod enat to have subtracted q/2 in the end */ 
  }
  if (sss==3) return tanhl(xxtmp);
  else return tanhl(xxtmp-0.5L);
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

