#include <stdio.h>
#include <math.h>
#include <string.h>
#define NL(flnm,dumch) {do dumch=fgetc(flnm); while (dumch!='\n' && dumch!=EOF);}
#define ROOT_NO 2000
#define SMALL (1.0e-15L)
#define TOOLONG 10000
#define EVEN(intval) (2*(intval/2)==intval)

const char *infilename="findres-es.in";
char outfilename[65];

void makeinitials(void);
void prepoutfiles(void);
void ntmapstep(long double *, long double *, long double *);
void ntdm(long double, long double, long double [2][2]);
long double findx0(long double);
long double rootfunceven1(long double);
long double rootfunceven2(long double);
long double rootfuncodd1(long double);
long double rootfuncodd2(long double);
int findshinofp(int, long double *, long double *, long double *, long double *);
long double wind(long double, long double);
long double residue(int, int, long double *);
int routing(int, unsigned long int, unsigned long int);
  
int zbrac(long double (*)(long double), long double *, long double, 
	  long double *, long double);
void zbrak(long double (*)(long double), long double, long double, int, long double[], long double[], int *);
long double rtbis(long double (*)(long double), long double, long double, long double);
void mkstr(char *, unsigned long int);


long double a0, b0, acc, lacc, brakl0, braku0, bracs, bracx, pi, lift=0.0L;
long double xorbs[5][ROOT_NO], yorbs[5][ROOT_NO];
int  ii, sln, brakn, omno; 
unsigned long int qq[100], pp[100], evalcount=0; 
void (*map)(long double *, long double *, long double *);
FILE *fl1;

main(int argc, char **argv)
{
  long double pery, phalf[2], qld;
  long double brackl[ROOT_NO], bracku[ROOT_NO];
  int maxbrack=ROOT_NO;
  int i, j, dir, steps, zbexit[ROOT_NO], meql[ROOT_NO];
  long double (*rootfunc)(long double);
 
  map=&ntmapstep;

  /* read input */
  makeinitials();

  /* prepare output files */
  fl1=fopen(outfilename, "w");
  prepoutfiles();

  /* start y-values on two sln's from lowest order qq/pp via zbrak-bracketing */
  qld=(long double) (qq[0]);  
  for (j=0;j<2;j++) {
      if (j==0) {
	  if (EVEN(pp[0])) {sln=1; rootfunc=&rootfunceven1;} 
	  else {sln=1; rootfunc=&rootfuncodd1;}
      }
      else {
	  if (EVEN(pp[0])) {sln=3; rootfunc=&rootfunceven2;} 
	  else {sln=2; rootfunc=&rootfuncodd1;}
      }

      maxbrack=ROOT_NO-1; 
      zbrak(rootfunc, brakl0, braku0, brakn, &brackl[0], &bracku[0], &maxbrack);
      meql[j]=0; for (i=1;i<=maxbrack;i++) {
	  pery=rtbis(rootfunc, brackl[i], bracku[i], acc);
	  lift=wind(findx0(pery), pery); 
	  if (fabsl(lift-qld)<0.1L) {
	      yorbs[sln][++meql[j]]=pery; 
	      xorbs[sln][meql[j]]=findx0(pery);
	  }
      }
      
      fprintf(fl1,"%15u %15u %4d %4d\n", qq[0], pp[0], meql[j], sln);
      for (i=1;i<=meql[j];i++) {
	  fprintf(fl1,"%20.15Lf %20.15Lf %20.12Lg\n", 
		  xorbs[sln][i], yorbs[sln][i], residue(i, sln, &phalf[i]));
      }
      fprintf(fl1,"\n"); fflush(fl1);

      sln=routing(sln, qq[0],pp[0]);
      fprintf(fl1,"%15u %15u %4d %4d\n", qq[0], pp[0], meql[j], sln);
      for (i=1;i<=meql[j];i++) {
	  yorbs[sln][i]=phalf[i];
	  xorbs[sln][i]=findx0(yorbs[sln][i]);
	  fprintf(fl1,"%20.15Lf %20.15Lf %20.12Lg\n", 
		  xorbs[sln][i], yorbs[sln][i], residue(i, sln, &phalf[i]));
      }
      fprintf(fl1,"\n"); fflush(fl1);
  }

  
  /* further y-values for higher order qq/pp from previous ones */
  for (ii=2;ii<omno;ii+=2) {
      qld=(long double) (qq[ii]);  

      /* find root number and brackets */
      maxbrack=0; 
      for (j=0;j<2;j++) {
	  if (j==0) {
	      if (EVEN(pp[ii])) {sln=1; rootfunc=&rootfunceven1;} 
	      else {sln=1; rootfunc=&rootfuncodd1;}
	  }
	  else {
	      if (EVEN(pp[ii])) {sln=3; rootfunc=&rootfunceven2;} 
	      else {sln=2; rootfunc=&rootfuncodd1;}
	  }
	  
	  for (i=1;i<=meql[j];i++) {
	      if (i==1) brackl[2]=brakl0; else brackl[2]=yorbs[sln][i-1];
	      if (i==meql[j]) bracku[2]=braku0; else bracku[2]=yorbs[sln][i+1];
	      brackl[3]=brackl[2]; bracku[3]=bracku[2];
	      brackl[0]=yorbs[sln][i]; bracku[0]=brackl[0]+bracs;
	      brackl[1]=yorbs[sln][i]-bracs; bracku[1]=brackl[1]+bracs;
	      dir=0; steps=0;
	      do {
		  zbexit[i]=zbrac(rootfunc, &brackl[dir], brackl[2], &bracku[dir], bracku[2]);
		  if (zbexit[i]==1) {
		      pery=rtbis(rootfunc, brackl[dir], bracku[dir], acc);
		      lift=wind(findx0(pery), pery);
		      if (fabsl(lift-qld)>0.1L) {
			  if (dir==1) {
			      if (bracku[0]+bracs<bracku[3]) 
			      {dir=0; brackl[2]=bracku[0]; brackl[0]=bracku[0]; bracku[0]+=bracs;}
			      else if (brackl[1]-bracs>brackl[3]) 
			      {bracku[2]=brackl[1]; bracku[1]=brackl[1]; brackl[1]-=bracs;}
			      else zbexit[i]=0;
			  }
			  else {
			      if (brackl[1]-bracs>brackl[3]) 
			      {dir=1; bracku[2]=brackl[1]; bracku[1]=brackl[1]; brackl[1]-=bracs;}
			      else if (bracku[0]+bracs<bracku[3]) 
			      {brackl[2]=bracku[0]; brackl[0]=bracku[0]; bracku[0]+=bracs;}
			      else zbexit[i]=0;
			  }
		      }
		  }
/*		  printf("%Lf %u %.0Lf %d\n", pery, dir, lift-qld, zbexit[i]); */
		  
	      }
	      while (zbexit[i]==1 && fabsl(lift-qld)>0.1L && steps++<TOOLONG);	      
	      if (fabsl(lift-qld)<0.1L) {
		  yorbs[sln][i]=pery;
		  xorbs[sln][i]=findx0(pery);
/*		  printf("%Lf %u %.0Lf %d %d %d ok\n",pery, dir, lift-qld, zbexit[i], sln, ii); */
	      }
	      else zbexit[i]=-fabsl(lift-qld);
/*	      else printf("%Lf %u %.0Lf %d %d %d !!!!!\n",pery, dir, lift-qld, zbexit[i], sln, ii); */
	  }

	  fprintf(fl1,"%15u %15u %4d %4d\n", qq[ii], pp[ii], meql[j], sln);
	  for (i=1;i<=meql[j];i++) {
	      if (zbexit[i]==1) {
		  fprintf(fl1,"%20.15Lf %20.15Lf %20.12Lg %4d %10u\n", 
			  xorbs[sln][i], yorbs[sln][i], 
			  residue(i, sln, &phalf[i]), zbexit[i], steps);
	      }
	      else 
		  fprintf(fl1,"%20.15Lf %20.15Lf %20.12Lg %4d %10u\n", 
			  0.0L, 0.0L, 0.0L, zbexit[i], steps);
	  }	  
	  fprintf(fl1,"\n"); 
	  
	  sln=routing(sln, qq[ii],pp[ii]);
	  fprintf(fl1,"%15u %15u %4d %4d\n", qq[ii], pp[ii], meql[j], sln);
	  for (i=1;i<=meql[j];i++) {
	      if (zbexit[i]==1) {
		  yorbs[sln][i]=phalf[i];
		  xorbs[sln][i]=findx0(yorbs[sln][i]);
		  fprintf(fl1,"%20.15Lf %20.15Lf %20.12Lg %4d %10u\n", 
			  xorbs[sln][i], yorbs[sln][i], 
			  residue(i, sln, &phalf[i]), zbexit[i], steps);
	      }
	      else 
		  fprintf(fl1,"%20.15Lf %20.15Lf %20.12Lg %4d %10u\n", 
			  0.0L, 0.0L, 0.0L, zbexit[i], steps);
	  }	  
	  fprintf(fl1,"\n"); fflush(fl1);
      }
  }
  
  fclose(fl1);
  printf("%u\n", evalcount);
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
  
  fscanf(fl,"%Lf", &a0); NL(fl, dummy);       /* a */
  fscanf(fl,"%Lf", &b0); NL(fl, dummy);       /* b */
  fscanf(fl,"%Le", &acc); NL(fl, dummy);      /* absolute root acc */
  fscanf(fl,"%Le", &bracs); NL(fl, dummy);    /* bracs */
  fscanf(fl,"%Lf", &brakl0); NL(fl, dummy);   /* lower end of interval for "brak"ket */
  fscanf(fl,"%Lf", &braku0); NL(fl, dummy);   /* upper end of interval for "brak"ket */
  fscanf(fl,"%d", &brakn); NL(fl, dummy);     /* divisions for "brak"ket interval */
  fscanf(fl,"%s", outfilename); NL(fl, dummy);/* output filename */
  fscanf(fl,"%d", &omno); NL(fl, dummy);       /* number of winding numbers below */
  for (i=0;i<omno;i++) fscanf(fl,"%u %u", &qq[i], &pp[i]); NL(fl, dummy);  /* (q,p)'s */
  
  fclose(fl);

  /* adjust initials */
  if (bracs==0.0L) bracs=SMALL;
  if (brakl0==braku0) {brakl0=-0.5L; braku0=0.5L;}
}


void prepoutfiles(void)
{
  fprintf(fl1,"# a:       %10.7Lf , b:      %10.7Lf \n", a0, b0); 
  fprintf(fl1,"# brakl0:  %10.7Lf , braku0:  %10.7Lf , brakn:    %10d\n", 
	  brakl0, braku0, brakn);  
  fprintf(fl1,"# acc:     %10.3Le , bracs:    %10.7Le\n", acc, bracs); 
  fprintf(fl1,"# %13s %15s %4s %4s\n", "qq", "pp", "#n", "sln"); 
  fprintf(fl1,"# %17s %20s %20s %4s%10s\n", "x0", "y0", "res", "ok", "steps"); 
}

/* single iteration of NT-map */
void ntmapstep(long double *x, long double *y, long double *lft)
{

  *y = *y - b0*sinl(2.0L*pi*(*x));
  *x = *x + a0*(1.0L-(*y)*(*y));
  
  /* x modulo 1 */
  while (*x> 0.995L) {*x=*x-1.0L; (*lft)++;}
  while (*x< -0.005L) {*x=*x+1.0L; (*lft)--;}
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
  while (tmp> 0.995L) tmp=tmp-1.0L;
  while (tmp< -0.005L) tmp=tmp+1.0L;
  return tmp;
}

/* sin(2 pi x) after n/2 map iterations for root search */
long double rootfunceven1(long double y0)
{
  int i;
  long double xx, yy, lft;

  xx=findx0(y0); yy=y0; lft=0;
  for (i=1;i<=pp[ii]/2;i++) map(&xx, &yy, &lft);
  evalcount++;
  return sinl(2.0L*pi*xx); /* I0 -> I0 after n/2 */
}

/* sin(2 pi x) after n/2 map iterations for root search */
long double rootfunceven2(long double y0)
{
  int i;
  long double xx, yy, lft;

  xx=findx0(y0); yy=y0; lft=0;
  for (i=1;i<=pp[ii]/2;i++) map(&xx, &yy, &lft);
  evalcount++;
  return sinl(2.0L*pi*(xx-a0/2.0L*(1.0L-yy*yy))); /* I1 -> I1 after n/2 */
}

/* sin(2 pi x) after n/2 map iterations for root search */
long double rootfuncodd1(long double y0)
{
  int i;
  long double xx, yy, lft;

  xx=findx0(y0); yy=y0; lft=0;
  for (i=1;i<=pp[ii]/2+1;i++) map(&xx, &yy, &lft);    /* n odd: n/2+1=(n+1)/2 */
  evalcount++;
  return sinl(2.0L*pi*(xx-a0/2.0L*(1.0L-yy*yy))); /* I0 -> I1 after n/2+1=(n+1)/2 */
}

/* sin(2 pi x) after n/2 map iterations for root search */
long double rootfuncodd2(long double y0)
{
  int i;
  long double xx, yy, lft;

  xx=findx0(y0); yy=y0; lft=0;
  for (i=1;i<=pp[ii]/2;i++) map(&xx, &yy, &lft); /* n odd: n/2=(n-1)/2 */
  evalcount++;
  return sinl(2.0L*pi*xx);                   /* I1 -> I0 after n/2=(n-1)/2 */
}



/* find winding number of known periodic orbit, with optional output to fl3 */
long double wind(long double q, long double p)
{
  long double qqq=q, ppp=p, lft=0.0L;
  unsigned long int i;

  lft=0;
  for (i=1;i<=pp[ii];i++) map(&qqq, &ppp, &lft);
  evalcount+=2;
  return lft;
}

/* find residues of a known periodic orbit, with opt. output to fl2 and *tmp on I0 */
long double residue(int i0, int j0, long double *tmp)
{
  long double qqq, ppp, lft=0.0L, dpq[2][2], lm[2][2], lm0[2][2];
  unsigned long int i, j, k;

  qqq=xorbs[j0][i0];
  ppp=yorbs[j0][i0];
  
  ntdm(qqq,ppp,lm);

  for (k=1;k<pp[ii];k++) {
    map(&qqq, &ppp, &lft);
    if (k==(pp[ii]+1)/2) *tmp=ppp;   /* good for all even orbits, and s1/s2 odd orbits: I0->I1 */
    ntdm(qqq,ppp,dpq);
    for (i=0;i<=1;i++) for (j=0;j<=1;j++) lm0[i][j]=lm[i][j]; 
    for (i=0;i<=1;i++) for (j=0;j<=1;j++) lm[i][j]=dpq[i][0]*lm0[0][j]+dpq[i][1]*lm0[1][j];
  }

  map(&qqq, &ppp, &lft);
  evalcount+=2;
  return (2.0L-lm[0][0]-lm[1][1])/4.0L;
}

int routing(int oldsl, unsigned long int nmr, unsigned long int dnm) 
{
    if (EVEN(dnm)){
	switch (oldsl) {
	    case 1: return 2;
	    case 2: return 1;
	    case 3: return 4;
	    case 4: return 3;
	    default: return 0;
	}
    } else {
	if (EVEN(nmr)){
	    switch (oldsl) {
		case 1: return 3;
		case 2: return 4;
		case 3: return 1;
		case 4: return 2;
		default: return 0;
	    }	    
	} else {
	    switch (oldsl) {
		case 1: return 4;
		case 2: return 3;
		case 3: return 2;
		case 4: return 1;
		default: return 0;
	    }
	}
    }
}


/* slightly modified version of zbrac from Num. Rec. (limits on brackets, x1<x2!) */
int zbrac(long double (*func)(long double), long double *x1, long double x1min, 
	     long double *x2, long double x2max)
{
  const int NTRY=50;
  const long double FACTOR=1.6L;
  const long double TOOBIG=1.0L/SMALL;
  int j;
  long double f1, f2, tmp;

  if (*x1==*x2) *x2+=SMALL;
  f1=(*func)(*x1);
  f2=(*func)(*x2);
  for (j=1;j<=NTRY;j++) { 
      if (f1*f2<0.0) return 1;
      tmp=FACTOR*(*x2 - *x1);
      if ((fabsl(f1) < fabsl(f2) || *x2+tmp>x2max) && *x1-tmp>=x1min)
	  f1=(*func)(*x1 -= tmp);
      else if (*x2+tmp<=x2max)
	  f2=(*func)(*x2 += tmp);
      else break;
  }
  return 0;
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
  const int JMAX=50;
  int j;
  long double dx, f, fmid, xmid, rtb;
  
  f=(*func)(x1);
  fmid=(*func)(x2);
  if (f*fmid >= 0.0) {printf("Root must be bracketed for bisection in rtbis!!!\n"); exit(1);} 
  rtb = f < 0.0 ? (dx=x2-x1, x1) : (dx=x1-x2, x2); 
  for (j=1;j<=JMAX;j++) {
    fmid=(*func)(xmid=rtb+(dx *= 0.5L)); 
    if (fmid <= 0.0) rtb=xmid;
    if (fabsl(dx) < xacc || fmid == 0.0) return rtb;
  }
  printf("Too many bisections in rtbis!!!\n");
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


/* routing pattern of periodic orbits:
   odd/even:    1<->2    3<->4
   odd/odd:     1<->4    2<->3
   even/odd:    1<->3    2<->4 */




