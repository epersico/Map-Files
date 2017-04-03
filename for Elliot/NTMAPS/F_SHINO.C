#include <stdio.h>
#include <math.h>
#include <string.h>
#define NL(flnm,dumch) {do dumch=fgetc(flnm); while (dumch!='\n' && dumch!=EOF);}

const char *infilename="f_shino.in";
const long double SMALL=1.0e-10L;
const int SUPPERROR=15;

void makeinitials(void);
void prepoutfiles(void);
void assignshino(long double, int);
void ntmapstep(long double *, long double *, long double);
long double xrootfunc(long double, long double *, long double *);
long double yrootfunc(long double);
 
int rootbybrak(long double (*)(long double, long double *, long double *), 
	       long double *, long double, long double, int, int);

void zbrakwrite(long double (*)(long double, long double *, long double *),  
	   long double, long double, int, long double[], long double[], int *);
long double rtbis(long double (*)(long double, long double *, long double *),   
		  long double, long double, long double);
void mkstr(char *, unsigned long int);

long double pi;
long double omega, a0, b0, bx, a00, b00, bx0, bxx, xacc, yacc;
long double b;
long double shinox, shinoyf, shinopx, shinopyf;
int shp, root_no, brakno, qp_no, qpi;
unsigned long int qq[500], pp[500];
char shinofilename[160], shinofilenamebase[80], shinofilenameend[20];

FILE *fl1, *fl2;

main(int argc, char **argv)
{
    int rootexit=1;
    
    /* read input */
    makeinitials();

for (qpi=0; qpi<qp_no; qpi++) {
	
    /* prepare output files */
    prepoutfiles();

    if (a0<=omega+SMALL) {b=0; rootexit=1;}
    else {
	assignshino(a0,shp);
	rootexit=rootbybrak(&xrootfunc, &b, b0, bx, brakno, root_no);
    }
    fprintf(fl2,"%15u%15u%20.15Lf%10d%20.10Lg\n", qq[qpi], pp[qpi], b, rootexit, yrootfunc(b));
    if (rootexit<=0)
	printf("Warning, no yroot found for %10u/%-10u among %d xroot candidates!\n",
	       qq[qpi],pp[qpi], -rootexit);

}
 fclose(fl2);
 fclose(fl1);
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
  
  fscanf(fl,"%Ld", &shp);  NL(fl, dummy);     /* shinohara point */
  fscanf(fl,"%Lf", &a00);  NL(fl, dummy);      /* a0 */
  fscanf(fl,"%Lf", &b00);  NL(fl, dummy);      /* b0 */
  fscanf(fl,"%Lf", &bx0);  NL(fl, dummy);      /* bmax */
  fscanf(fl,"%Lf", &bxx);  NL(fl, dummy);     /* bmaxx */
  fscanf(fl,"%d", &brakno);  NL(fl, dummy);   /* # of zbrak intervals */
  fscanf(fl,"%Le", &xacc);  NL(fl, dummy);    /* absolute root xacc */
  fscanf(fl,"%Le", &yacc);  NL(fl, dummy);    /* absolute root yacc */
  fscanf(fl,"%d", &root_no); NL(fl, dummy);   /* max "brak" intervals stored */
  fscanf(fl,"%s %s", shinofilenamebase, shinofilenameend); NL(fl, dummy);   /* outfilename */
  fscanf(fl,"%d", &qp_no);  NL(fl, dummy);    /* number of p/q values below */
  for (i=0;i<qp_no;i++) {
      fscanf(fl,"%u %u", &qq[i], &pp[i]); NL(fl, dummy);                     /* q/p */
  }
  
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
  if (bx0<b0) bx=b0; else bx=bx0;
  if (brakno<=1) brakno=1;
  if (yacc<xacc) yacc=10.0L*xacc;
  if (bx>bxx) bx=bxx;

  strcpy(shinofilename,shinofilenamebase);
  strcpy(ctmp,""); mkstr(&ctmp[0],qq[qpi]); strcat(shinofilename,ctmp);
  strcat(shinofilename,"_");
  strcpy(ctmp,""); mkstr(&ctmp[0],pp[qpi]); strcat(shinofilename,ctmp);
  strcat(shinofilename,shinofilenameend);

  fl1=fopen(shinofilename, "w");
  fprintf(fl1,"# Sh. pt.: %2d \n", shp);
  fprintf(fl1,"# q:       %10u , p:       %10u \n", qq[qpi], pp[qpi]);
  fprintf(fl1,"# a0:      %10.7Lf , b0:      %10.7Lf \n", a0, b0); 
  fprintf(fl1,"# bx:      %10.7Lf , bxx:     %10.7Lf \n", bx, bxx); 
  fprintf(fl1,"# root#:   %10d \n", root_no); 
  fprintf(fl1,"# %18s%25s%25s%25s%25s\n\n", 
	  "b", "xrootfunc", "smooth xrf", "yrootfct", "xyrootfunc");

  if (qpi==0) {
      strcpy(shinofilename,shinofilenamebase);
      strcat(shinofilename,shinofilenameend);
      fl2=fopen(shinofilename, "w");
      fprintf(fl2,"# Sh. pt.: %2d \n", shp);
      fprintf(fl2,"# a0:      %10.7Lf , b0:      %10.7Lf \n", a00, b00); 
      fprintf(fl2,"# bx:      %10.7Lf , bxx:     %10.7Lf \n", bx0, bxx); 
      fprintf(fl2,"# brakno:  %10d \n", brakno); 
      fprintf(fl2,"# xacc:    %10.3Le , yacc:    %10.3Le \n", xacc, yacc); 
      fprintf(fl2,"# root#:   %10d \n", root_no); 
      fprintf(fl2,"# %13s%15s%20s%10s%20s\n\n", "q", "p", "b", "exitcode", "yrootfct");
  }
}

/* assign shinohara points or prefactors depending on symmetry line and a */
void assignshino(long double aa, int ss) 
{
	switch (ss) {
	    case 1: {
		shinox=-0.25L; shinoyf=-0.5L; 
		shinopx=-0.25L+0.5L; shinopyf=((2*(qq[qpi]/2)==qq[qpi]) ? -0.5L : 0.5L); 
		break;}
	    case 2: {
		shinox=aa/2.0L-0.25L; shinoyf=0.0L; 
		shinopx=aa/2.0L-0.25L+0.5L; shinopyf=0.0L; 
		break;}
	    case 3: {
		shinox=0.25L; shinoyf=0.5L; 
		shinopx=0.25L+0.5L; shinopyf=((2*(qq[qpi]/2)==qq[qpi]) ? -0.5L : 0.5L);
		break;}
	    case 4: {
		shinox=aa/2.0L+0.25L; shinoyf=0.0L; 
		shinopx=aa/2.0L+0.25L+0.5L; shinopyf=0.0L; 
		break;}
	    default: {
		shinox=-0.25L; shinoyf=-0.5L;
		shinopx=aa/2.0L+0.25L+0.5L; shinopyf=0.0L;}
	}
}



/* single iteration of NT-map */
void ntmapstep(long double *x, long double *y, long double bb)
{

  *y = *y - bb*sinl(2.0L*pi*(*x));
  *x = *x + a0*(1.0L-(*y)*(*y));
  
  /* x modulo 1 
  while (*x> 0.5001L) {*x=*x-1.0L; (*lft)++;}
  while (*x< -0.5L) {*x=*x+1.0L; (*lft)--;} */
}


/*  p/2 map iterations minus shino(m+q) for root search */
long double xrootfunc(long double bb, long double *yyrf, long double *xxyyrf)
{
  long double i, j=1.0L;
  long double xx, yy;

  xx=shinox; yy=shinoyf*bb;
  for (i=1.0L;i<=(long double) (pp[qpi]/2);i++) {
      ntmapstep(&xx, &yy, bb);
      if (i>=j/omega){
	  xx-=1.0L;
	  j++;
      }
  }
  *yyrf=(yy-shinopyf*bb); 
  *xxyyrf=(xx-shinopx)*(xx-shinopx)+(yy-shinopyf*bb)*(yy-shinopyf*bb);
  return (xx-shinopx); 
}

/*  p/2 map iterations minus shino(m+q) for root search */
long double yrootfunc(long double bb)
{
  long double i, j=1.0L;
  long double xx, yy;

  xx=shinox; yy=shinoyf*bb;
  for (i=1.0L;i<=(long double) (pp[qpi]/2);i++) {
      ntmapstep(&xx, &yy, bb);
      if (i>=j/omega){
	  xx-=1.0L;
	  j++;
      }
  }
  return (yy-shinopyf*bb); 
}

/* zbrak bracketing and root search and test */
int rootbybrak(long double (*func)(long double, long double *, long double *), long double *bb,
	       long double x1, long double x2, int intervals, int maxroots)
{
  long double brackl[maxroots], bracku[maxroots];
  int i, maxbrack=maxroots;

  zbrakwrite(func, x1, x2, intervals, &brackl[0], &bracku[0], &maxbrack);
  i=1; 
  do {
      *bb=rtbis(func, brackl[i], bracku[i], xacc);
      i++;
  }
  while (i<=maxbrack && fabsl(yrootfunc(*bb))>yacc);
  if (fabsl(yrootfunc(*bb))>yacc) return -i;
  else return i;
}

/* zbrak from Num. Rec., modified to write output at each step */
void zbrakwrite(long double (*fx)(long double, long double *, long double *), 
		long double x1, long double x2, int nn, long double xb1[], long double xb2[], int *nb)
{
  int nbb, i;
  long double x, fp, fc, dx, xsrf, yrf, xyrf, tmp=10.0L;
  
  nbb=0;
  dx=(x2-x1)/nn;
  fp=(*fx)(x=x1, &yrf, &xyrf);
  for (i=1;i<=nn;i++) {
      fc=(*fx)(x += dx, &yrf, &xyrf); 
      xsrf=fc; if (fabsl(xsrf)>10.0L) xsrf=tmp; else tmp=xsrf;
    fprintf(fl1,"%20.15Lf%25.15Lg%25.15Lg%25.15Lg%25.15Lg\n", x, fc, xsrf, yrf, xyrf);
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
long double rtbis(long double (*func)(long double, long double *, long double *), 
		  long double x1, long double x2, long double xxacc)
{
  const int JMAX=40;
  int nrerror;
  int j;
  long double dx, f, fmid, xmid, rtb, dummy1, dummy2;

  f=(*func)(x1, &dummy1, &dummy2);
  fmid=(*func)(x2, &dummy1, &dummy2);
  if (f*fmid >= 0.0) nrerror=1; /* "Root must be bracketed for bisection in rtbis" */
  rtb = f < 0.0 ? (dx=x2-x1, x1) : (dx=x1-x2, x2); 
  for (j=1;j<=JMAX;j++) {
    fmid=(*func)(xmid=rtb+(dx *= 0.5L), &dummy1, &dummy2); 
    if (fmid <= 0.0) rtb=xmid;
    if (fabsl(dx) < xxacc || fmid == 0.0) return rtb;
  }
  nrerror=2; /* "Too many bisections in rtbis" */
  return 0.0L;
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

