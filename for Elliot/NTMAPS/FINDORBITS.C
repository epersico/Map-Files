#include <stdio.h>
#include <math.h>
#define ROOT_NO 10000

const char *infilename="findorbits.in", 
  *outfilename="/tmp/findorbits.out", 
  *orbitfilename="/tmp/findorbits.data", 
  *allorbitfilename="/tmp/findorbits.all.data";

void makeinitials(void);
void prepoutfiles(void);
void harpermapstep(long double *, long double *, int *);
void ntmapstep(long double *, long double *, int *);
void harperdm(long double, long double, long double[2][2]);
void ntdm(long double, long double, long double[2][2]);
int b0root(long double *, long double *);
long double findx0(long double);
long double rootfunc(long double);
int wind(long double, long double);
long double residue(long double, long double, long double *);
 
void mybrac(long double (*)(long double), long double, long double, long double, 
	    long double, int, long double, long double[], long double[], int *);
int zbracdir(long double (*)(long double), long double, long double *);
void zbrak(long double (*)(long double), long double, long double, int, long double[], long double[], int *);
long double rtbis(long double (*)(long double), long double, long double, long double);

long double pi;
long double a, b;
long double a0, b0, da, db, xa, xb, acc, brakl0, braku0, bracs, bracx;
int  mt, sln, n, m, brakn, outfiles;
int lift=0;
void (*map)(long double *, long double *, int *);
FILE *fl1, *fl2, *fl3;

main(int argc, char **argv)
{
  extern long double pi;
  extern long double a, b, a0, b0, da, db, xa, xb, acc, brakl0, braku0, bracs, bracx;
  extern int  mt, n, m, brakn, lift, outfiles;
  extern FILE *fl1, *fl2, *fl3;

  long double pery, phalf;
  long double brackl[ROOT_NO], bracku[ROOT_NO], yorbs[ROOT_NO];
  int maxbrack=ROOT_NO;
  int i, j, meql=0;
  
  /* read input */
  makeinitials();

  /* prepare output files */
  fl1=fopen(outfilename, "w");
  if (outfiles>1) fl2=fopen(orbitfilename, "w");
  if (outfiles>2) fl3=fopen(allorbitfilename, "w");
  prepoutfiles();
    
  /* main loop */
  for (a=a0;a<=xa;a+=da)
    {
      meql=0;
      if (braku0<=brakl0) if (b0root(&yorbs[1],&yorbs[2])==1) meql=2;
      for (b=b0;b<=xb;b+=db)
	{    
	  fprintf(fl1,"%15.11Lf %15.11Lf", a, b);
 	  if (outfiles>1) fprintf(fl2,"#\n# a=%10.7Lf\n", a);
	  if (outfiles>1) fprintf(fl2,"# %8s%12s%10s    etc.\n", "b", "q0", "p0");
 	  if (outfiles>2) fprintf(fl3,"#\n# a=%10.7Lf\n", a);
	  if (outfiles>2) fprintf(fl3,"# %8s%12s%10s    etc.\n", "b", "q0", "p0");

	  /* find root number and brackets */
	  if (bracs==0.0 || (meql==0 && braku0>brakl0)) {
	    maxbrack=ROOT_NO;
	    zbrak(&rootfunc, brakl0, braku0, brakn, &brackl[0], &bracku[0], &maxbrack);
	  }
	  else {
	    maxbrack=0; 
	    for (i=1;i<=meql;i++) {
	      j=maxbrack+1; 
	      maxbrack=(ROOT_NO/meql)*i;
	      mybrac(&rootfunc, 
		     ((i<meql) ? ((yorbs[i+1]-yorbs[i])/2.0L) : (bracx)),
		     ((i>1) ? ((yorbs[i]-yorbs[i-1])/2.0L) : (bracx)), bracx, 
		     bracs, j, yorbs[i], &brackl[0], &bracku[0], &maxbrack);
	    }

	  }
	  	  
	  fprintf(fl1,"%4d", maxbrack);  /* number of roots found */
	  meql=0;                        /* to count number of roots with lift m */

	  /* examine root brackets */
	  for (i=1;i<=maxbrack;i++) {

	    /* find roots */
	    pery=rtbis(&rootfunc, brackl[i], bracku[i], acc);

	    /* count and store ?/n roots with m/n winding numbers */
	    lift=wind(findx0(pery), pery); /* optional output to fl3 */
	    if (lift==m) yorbs[++meql]=pery;
	  }
	  
	  fprintf(fl1,"%4d", meql);
	  for (i=1;i<=meql;i++) {
	    fprintf(fl1,"%15.7Lf%12.7Lf%8.3Lf", 
		    findx0(yorbs[i]), yorbs[i], residue(findx0(yorbs[i]), yorbs[i], &phalf));
	    fprintf(fl1,"%12.7Lf", phalf);
	  }
      
	  fprintf(fl1,"\n");
	  if (outfiles>1) fprintf(fl2,"\n");
	  if (outfiles>2) fprintf(fl3,"\n");
	}
    }
  
  fclose(fl1);
  if (outfiles>1) fclose(fl2);
  if (outfiles>2) fclose(fl3);
  exit(0); 
} 


void makeinitials(void)      /* read and adjust initial values */
{
  FILE *fl;
  char dummy[80];
  extern long double a0, b0, da, db, xa, xb, acc, brakl0, braku0, bracs, bracx;
  extern int  mt, n, m, brakn, outfiles;
  extern const char *infilename;
  extern void (*map)(long double *, long double *, int *);
  extern void harpermapstep(long double*, long double *, int *);
  extern void ntmapstep(long double*, long double *, int *);  
  
  extern long double pi;
  pi=4.0L*atanl(1.0L);

  /* read initials */
  fl=fopen(infilename, "r");
  
  fscanf(fl,"%ld", &mt); fgets(dummy, 80, fl);       /* maptype */
  fscanf(fl,"%ld", &sln); fgets(dummy, 80, fl);      /* symm. line # */
  fscanf(fl,"%Lf", &a0); fgets(dummy, 80, fl);       /* start a */
  fscanf(fl,"%Lf", &b0); fgets(dummy, 80, fl);       /* start b */
  fscanf(fl,"%ld", &n); fgets(dummy, 80, fl);        /* period n */
  fscanf(fl,"%ld", &m); fgets(dummy, 80, fl);        /* m */
  fscanf(fl,"%Le", &acc); fgets(dummy, 80, fl);      /* absolute root acc */
  fscanf(fl,"%Le", &bracs); fgets(dummy, 80, fl);    /* size of "brac"ket start guess */
  fscanf(fl,"%Lf", &bracx); fgets(dummy, 80, fl);    /* maximal "brac"ket shift in search*/
  fscanf(fl,"%Lf", &brakl0); fgets(dummy, 80, fl);   /* lower end of interval for "brak"ket */
  fscanf(fl,"%Lf", &braku0); fgets(dummy, 80, fl);   /* upper end of interval for "brak"ket */
  fscanf(fl,"%ld", &brakn); fgets(dummy, 80, fl);    /* divisions for "brak"ket interval */
  fscanf(fl,"%Lf", &da); fgets(dummy, 80, fl);       /* da */
  fscanf(fl,"%Lf", &db); fgets(dummy, 80, fl);       /* db */
  fscanf(fl,"%Lf", &xa); fgets(dummy, 80, fl);       /* max a */
  fscanf(fl,"%Lf", &xb); fgets(dummy, 80, fl);       /* max b */
  fscanf(fl,"%ld", &outfiles); fgets(dummy, 80, fl); /* number of output files */
  
  fclose(fl);

  /* adjust initials */
  if (xa<a0) xa=a0;
  if (xb<b0) xb=b0;
  if (da==0.0L) da=0.00000001L;
  if (db==0.0L) db=0.00000001L;
  if ((bracs==0.0L) && (brakl0==braku0)) {brakl0=-0.5L; braku0=0.5L;}
  mt==0 ? (map=&harpermapstep) : (map=&ntmapstep);
}


void prepoutfiles(void)
{
  extern long double a, b, a0, b0, da, db, xa, xb, acc, brakl0, braku0, bracs;
  extern int  mt, n, m, brakn, outfiles;
  extern FILE *fl1, *fl2, *fl3;

  fprintf(fl1,"# Maptype: %2d \n", mt);
  fprintf(fl1,"# Symmline:%2d \n", sln);
  fprintf(fl1,"# a0:      %10.7Lf , b0:      %10.7Lf \n", a, b); 
  fprintf(fl1,"# da:      %10.7Lf , db:      %10.7Lf \n", da, db); 
  fprintf(fl1,"# xa:      %10.7Lf , xb:      %10.7Lf \n", xa, xb); 
  fprintf(fl1,"# Period:   (m =%4d) / (n =%4d) \n", m, n); 
  fprintf(fl1,"# brakl0:  %10.7Lf , braku0:  %10.7Lf , brakn:    %10d\n", 
	  brakl0, braku0, brakn); 
  fprintf(fl1,"# acc:     %10.3Le , bracs:   %10.3Le , bracx:    %10.7Lf\n", acc, bracs, bracx); 
  fprintf(fl1,"# %13s%15s%4s%4s%12s%12s%8s%12s%30s\n", 
	  "a", "b", "#n", "#w", "x0", "y0", "res", "y0", "on I0 for m/n-orbits..."); 

  if (outfiles>1) {
    fprintf(fl2,"# Maptype: %2d \n", mt);
    fprintf(fl2,"# Symmline:%2d \n", sln);
    fprintf(fl2,"# a0:      %10.7Lf , b0:      %10.7Lf \n", a, b); 
    fprintf(fl2,"# da:      %10.7Lf , db:      %10.7Lf \n", da, db); 
    fprintf(fl2,"# xa:      %10.7Lf , xb:      %10.7Lf \n", xa, xb); 
    fprintf(fl2,"# Period:   (m =%4d) / (n =%4d) \n", m, n); 
    fprintf(fl2,"# brakl0:  %10.7Lf , braku0:  %10.7Lf , brakn:    %10d\n", 
	    brakl0, braku0, brakn); 
    fprintf(fl2,"# acc:     %10.3Le , bracs:   %10.3Le , bracx:    %10.7Lf\n", acc, bracs, bracx); 
  }
  
  if (outfiles>2) {
    fprintf(fl3,"# Maptype: %2d \n", mt);
    fprintf(fl3,"# Symmline:%2d \n", sln);
    fprintf(fl3,"# a0:      %10.7Lf , b0:      %10.7Lf \n", a, b); 
    fprintf(fl3,"# da:      %10.7Lf , db:      %10.7Lf \n", da, db); 
    fprintf(fl3,"# xa:      %10.7Lf , xb:      %10.7Lf \n", xa, xb); 
    fprintf(fl3,"# Period:   (m =%4d) / (n =%4d) \n", m, n); 
    fprintf(fl3,"# brakl0:  %10.7Lf , braku0:  %10.7Lf , brakn:    %10d\n", 
	    brakl0, braku0, brakn); 
    fprintf(fl3,"# acc:     %10.3Le , bracs:   %10.3Le , bracx:    %10.7Lf\n", acc, bracs, bracx); 
  }
}


/* single iteration of Harper map */
void harpermapstep(long double *x, long double *y, int *lft)
{
  extern long double a, b;
  extern long double pi;

  *y = *y + b*sinl(2.0L*pi*(*x));
  *x = *x - a*sinl(2.0L*pi*(*y));
  
  /* x modulo 1 */
  while (*x> 0.5001L) {*x=*x-1.0L; (*lft)++;}
  while (*x< -0.5L) {*x=*x+1.0L; (*lft)--;}
  /* y modulo 1 */
  while (*y> 0.5001L) *y=*y-1.0L;
  while (*y< -0.5L) *y=*y+1.0L;
}


/* single iteration of NT-map */
void ntmapstep(long double *x, long double *y, int *lft)
{
  extern long double a, b;
  extern long double pi;

  *y = *y - b*sinl(2.0L*pi*(*x));
  *x = *x + a*(1.0L-(*y)*(*y));
  
  /* x modulo 1 */
  while (*x> 0.5001L) {*x=*x-1.0L; (*lft)++;}
  while (*x< -0.5L) {*x=*x+1.0L; (*lft)--;}
  /* y modulo 1
  while (*y> 0.5) *y=*y-1.0;
  while (*y< -0.5) *y=*y+1.0; */
}


/* functional determinant of Harper map */
void harperdm(long double x, long double y, long double dxy[2][2])
{
  extern long double a, b;
  extern long double pi;
  long double tmp1, tmp2;

  tmp1=-2.0L*pi*a*cosl(2.0L*pi*(y+b*sinl(2.0L*pi*x)));
  tmp2=2.0L*pi*b*cosl(2.0L*pi*x);
  dxy[0][0]=1.0L+tmp1*tmp2;
  dxy[0][1]=tmp1;
  dxy[1][0]=tmp2;
  dxy[1][1]=1.0L;
}


/* functional determinant of NT-map */
void ntdm(long double x, long double y, long double dxy[2][2])
{
  extern long double a, b;
  extern long double pi;
  long double tmp1, tmp2;

  tmp1=-2.0L*a*(y-b*sinl(2.0L*pi*x));
  tmp2=-2.0L*pi*b*cosl(2.0L*pi*x);
  dxy[0][0]=1.0L+tmp1*tmp2;
  dxy[0][1]=tmp1;
  dxy[1][0]=tmp2;
  dxy[1][1]=1.0L;
}


/* directly calculate the root value for b=0 */
int b0root(long double *yl, long double *yu)
{
  extern long double a;
  extern int n, m, mt;
  extern long double pi;

  if (((long double)(m)) / ((long double)(n)) /a>1.0) {*yl=0; *yu=0; return 0;}
  else {
    if (mt==0) {
      *yu=-asinl( ((long double)(m)) / ((long double)(n)) /a)/2.0L/pi; 
      *yl=-0.5L-(*yu);
    }
    else {
      *yu=sqrtl(1.0L- ((long double)(m)) / ((long double)(n)) /a); 
      *yl=-(*yu);
    }
    return 1;
  }  
}


/* find x0(y0) so that point lies on symmetry line */
long double findx0(long double y0)
{
  extern long double a;
  extern int n, mt;
  extern long double pi;
  long double tmp;

  if (mt==0){  
    switch (sln) {
    case 1: {tmp=0.0L; break;}
    case 2: {tmp=0.5L; break;}
    case 3: {tmp=-a/2.0L*sinl(2.0L*pi*y0); break;}
    case 4: {tmp=-a/2.0L*sinl(2.0L*pi*y0)+0.5L; break;}
    default: {tmp=0.0L; printf("Uhoh - sln...\n");}
    }
  }
  else{
    switch (sln) {
    case 1: {tmp=0.0L; break;}
    case 2: {tmp=0.5L; break;}
    case 3: {tmp=a/2.0L*(1.0L-y0*y0); break;}
    case 4: {tmp=a/2.0L*(1.0L-y0*y0)+0.5L; break;}
    default: {tmp=0.0L; printf("Uhoh - sln...\n");}
    }
  }
  
  /* x modulo 1 */
  while (tmp> 0.5001L) tmp=tmp-1.0L;
  while (tmp< -0.5L) tmp=tmp+1.0L;
  return tmp;
}


/* sin(2 pi x) after n/2 map iterations for root search */
long double rootfunc(long double y0)
{
  extern long double pi;
  extern void (*map)(long double *, long double *, int *);
  extern int lift;

  int i;
  long double xx, yy;

  xx=findx0(y0); yy=y0; lift=0;
  for (i=1;i<=n/2;i++) map(&xx, &yy, &lift); /* n odd: n/2=(n-1)/2 */

  if (2*(n/2)==n) {                          /* n even */
    if (sln<=2) return sinl(2.0L*pi*xx);     /* I0 -> I0 after n/2 */
    else {                                   /* I1 -> I1 after n/2 */
      if (mt==0) return sinl(2.0L*pi*(xx+a/2.0L*sinl(2.0L*pi*yy)));
      else return sinl(2.0L*pi*(xx-a/2.0L*(1.0L-yy*yy)));
    }
  }
  else{                                       /* n odd */
    if (sln<=2) {                             /* I0 -> I1 after n/2+1=(n+1)/2 */
      map(&xx, &yy, &lift);
      if (mt==0) return sinl(2.0L*pi*(xx+a/2.0L*sinl(2.0L*pi*yy)));
      else return sinl(2.0L*pi*(xx-a/2.0L*(1.0L-yy*yy)));
    }
    else return sinl(2.0L*pi*xx);             /* I1 -> I0 after n/2=(n-1)/2 */
  }
}


/* find winding number of known periodic orbit, with optional output to fl3 */
int wind(long double q, long double p)
{
  extern int n, outfiles;
  extern long double b;
  extern void (*map)(long double *, long double *, int *);
  extern FILE *fl3;

  long double qq=q, pp=p;
  int lft=0, i;

  lft=0;
  if (outfiles>2) fprintf(fl3,"%10.6Lf", b);
  for (i=1;i<=n;i++) {if (outfiles>2) fprintf(fl3,"%12.6Lf %10.6Lf", qq, pp); map(&qq, &pp, &lft);}
  if (outfiles>2) fprintf(fl3,"%12.6Lf%10.6Lf%5d%10.2Le%10.2Le\n",qq,pp,lft,qq-q, pp-p);
  return lft;
}


/* find residues of a known periodic orbit, with opt. output to fl2 and *tmp on I0 */
long double residue(long double q, long double p, long double *tmp)
{
  extern int n, mt, outfiles;
  extern void (*map)(long double *, long double *, int *);
  extern FILE *fl1,*fl2;

  long double qq=q, pp=p, dpq[2][2], lm[2][2], lm0[2][2];
  int i, j, k, lft=0;

  if (mt==0) harperdm(qq,pp,lm); else ntdm(qq,pp,lm);

  if (outfiles>1) fprintf(fl2,"%10.6Lf", b);

  for (k=1;k<n;k++) {

    if (outfiles>1) fprintf(fl2,"%12.6Lf %10.6Lf", qq, pp); 

    map(&qq, &pp, &lft);
    if (k==n/2) *tmp=pp;
    if (mt==0) harperdm(qq,pp,dpq); else ntdm(qq,pp,dpq);
    for (i=0;i<=1;i++) for (j=0;j<=1;j++) lm0[i][j]=lm[i][j]; 
    for (i=0;i<=1;i++) for (j=0;j<=1;j++) 
      lm[i][j]=dpq[i][0]*lm0[0][j]+dpq[i][1]*lm0[1][j];
  }

  if (outfiles>1) fprintf(fl2,"%12.6Lf%10.6Lf%12.2Le%10.2Le\n", qq, pp, qq-q, pp-p);

  return (2.0L-lm[0][0]-lm[1][1])/4.0L;
}

/* zbrak-similar search for multiple roots with "directional" zbrac  version*/
void mybrac(long double (*func)(long double), 
	    long double maxup, long double maxdown, long double maxout, long double step, int nb0, 
	    long double xb0, long double xb1[], long double xb2[], int *nb)
{
  int nbb=nb0;
  
  xb1[nbb]=xb0;
  xb2[nbb]=xb0+bracs;
  while (1==1) {
    if (zbracdir(&rootfunc, xb1[nbb], &xb2[nbb])==0) break;
    if (xb2[nbb]>xb0+ ((maxup>maxout) ? (maxout) : (maxup))) break;
    if ((++nbb)>nb0+((*nb)-nb0)/2) break;
    xb1[nbb]=xb2[nbb-1];
    xb2[nbb]=xb2[nbb-1]+step;
  }
  xb1[nbb]=xb0-bracs;
  xb2[nbb]=xb0;
  while (1==1) {
    if (zbracdir(&rootfunc, xb2[nbb], &xb1[nbb])==0) break;
    if (xb1[nbb]<xb0- ((maxdown>maxout) ? (maxout) : (maxdown))) break;
    if ((++nbb)>(*nb)) break;
    xb1[nbb]=xb1[nbb-1]-step;
    xb2[nbb]=xb1[nbb-1];
  }
  *nb=nbb-1;
}



/* slightly modified version of zbrac from Num. Rec. (only one direction) */
int zbracdir(long double (*func)(long double), long double x1, long double *x2)
{
  const int NTRY=50;
  const long double FACTOR=1.6L;
  const long double TOOBIG=1.0e03L;
  int j, nrerror;
  long double f1, f2;

  if (x1==*x2) nrerror=1;
  f1=(*func)(x1);
  f2=(*func)(*x2);
  for (j=1;j<=NTRY;j++) {
    if (f1*f2<0.0) return 1;
    f2=(*func)(*x2 += FACTOR*(*x2-x1));
    if (fabs(*x2)>TOOBIG) break;
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
  const int JMAX=40;
  int nrerror;
  int j;
  long double dx, f, fmid, xmid, rtb;
  
  f=(*func)(x1);
  fmid=(*func)(x2);
  if (f*fmid >= 0.0) nrerror=1; /* "Root must be bracketed for bisection in rtbis" */
  rtb = f < 0.0 ? (dx=x2-x1, x1) : (dx=x1-x2, x2); 
  for (j=1;j<=JMAX;j++) {
    fmid=(*func)(xmid=rtb+(dx *= 0.5L)); 
    if (fmid <= 0.0) rtb=xmid;
    if (fabs(dx) < xacc || fmid == 0.0) return rtb;
  }
  nrerror=2; /* "Too many bisections in rtbis" */
  return 0.0L;
}
