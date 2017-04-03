#include <stdio.h>
#include <math.h>
#include <string.h>
#include <float.h>
#include <stdlib.h>

#define NL(flnm,dumch) {do dumch=fgetc(flnm); while (dumch!='\n' && dumch!=EOF);}
#define ntmapstep(x,y,a,b) {y-=(b)*sinl(pi2*(x)); x+=(a)-(a)*(y)*(y);}
#define smoothen(x,y) {if ((x)>0.0L) y=((x)/(1.0L+(x))); else y=((x)/(1.0L-(x)));}
#define SHFT(a,b,c,d) (a)=(b);(b)=(c);(c)=(d); /* for mnbrac & brent */
#define SIGN(a,b) ((b) >= 0.0 ? fabsl(a) : -fabsl(a))
#define FMAX(a,b) ((a) > (b) ? (a) : (b))
#define IPOW3(a) ((int)(pow(3.0, (double)(a))+0.5)) 
#define ERROR 123.456789L
#define UNDEFINED -123456.789L

const char *infilename="bca1.in";
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
int mymnbrac(long double *, long double *, long double *, long double,
	    long double (*)(long double, int, unsigned long int, unsigned long int), 
	    int, unsigned long int, unsigned long int);
long double brent(long double, long double, long double, 
	    long double (*)(long double, int, unsigned long int, unsigned long int), 
		  int, unsigned long int, unsigned long int, long double *);
void mkstr(char *, unsigned long int);
int qsortcompare(const void *, const void *);

long double pi2;
long double omega, a0, a00, a0min, a0max, a0step, b0, b00, ymin, ymax, ystep, y00, bxx;
long double bracdist, bracdistmin, INTERVALS, bsharp;
int sln, slnmin, slnmax, qpi, qp_no, bres, xbranch;
unsigned long int qq[500], pp[500], itstime2mod;
char iterfilename[160], iterfilenamebase[80];
long double braku[3], brakl[3];

FILE *fl1;
long double (*fiterfunc)(long double, long double, int, unsigned long int, unsigned long int);

int main(int argc, char **argv)
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

      maincalc();
      
  }
  }
  
  fclose(fl1);
  exit(0); 
} 


void makeinitials(void)      /* read and adjust initial values */
{
  FILE *fl;
  char dummy;
  int i, iintervals;
  long double tmp;
  
  pi2=8.0L*atanl(1.0L);

  /* read initials */
  fl=fopen(infilename, "r");
  
  fscanf(fl,"%d", &sln); NL(fl, dummy);                /* shinohara point */
  fscanf(fl,"%Lf %Lf", &a00, &a0max);  NL(fl, dummy);  /* a0-range*/
  fscanf(fl,"%Lf", &a0step);  NL(fl, dummy);           /* a0-step */
  fscanf(fl,"%Lf %Lf", &ymin, &ymax);  NL(fl, dummy);  /* y-range */
  fscanf(fl,"%Lf", &ystep);  NL(fl, dummy);            /* ystep */
  fscanf(fl,"%Lf", &b00);  NL(fl, dummy);              /* min b */
  fscanf(fl,"%Lf", &bxx);  NL(fl, dummy);              /* max b */
  fscanf(fl,"%d %d", &bres, &xbranch);  NL(fl, dummy); /* # of brackets, # of extrema branches*/
  fscanf(fl,"%d", &iintervals);  NL(fl, dummy);        /* # of mymnbrak-"intervals" */
  fscanf(fl,"%Lf", &bsharp);  NL(fl, dummy);           /* # sharpness of boundaries */
  fscanf(fl,"%Le %Le", &bracdist, &bracdistmin);  NL(fl, dummy);  /* for boundary estimates */
  fscanf(fl,"%lu", &itstime2mod);  NL(fl, dummy);      /* mod intervals threshold */
  fscanf(fl,"%s", iterfilenamebase); NL(fl, dummy);    /* outfilename */
  fscanf(fl,"%d", &qp_no);  NL(fl, dummy);    /* number of p/q values below */
  for (i=0;i<qp_no;i++) {
      fscanf(fl,"%lu %lu", &qq[i], &pp[i]); NL(fl, dummy);                     /* q/p */
  }

  INTERVALS=(long double)(iintervals);
  if (sln==0) {slnmin=1; slnmax=4;}
  else {slnmin=sln; slnmax=sln;}
  if (bxx<b00) {tmp=b00; b00=bxx; bxx=tmp;}
  if (bres<1) bres=1;
  if (xbranch<0) xbranch=0;
  if (bracdistmin>bracdist) bracdistmin=bracdist;
  fclose(fl);
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

  strcpy(iterfilename,iterfilenamebase);
  strcpy(ctmp,""); mkstr(&ctmp[0],qq[qpi]); strcat(iterfilename,ctmp);
  strcat(iterfilename,"_");
  strcpy(ctmp,""); mkstr(&ctmp[0],pp[qpi]); strcat(iterfilename,ctmp);
  strcat(iterfilename,".s");
  strcpy(ctmp,""); mkstr(&ctmp[0],sln); strcat(iterfilename,ctmp);

  fl1=fopen(iterfilename, "w");
  fprintf(fl1,"# Symm.l.: %2d \n", sln);
  fprintf(fl1,"# q:       %10lu , p:       %10lu \n", qq[qpi], pp[qpi]);
  fprintf(fl1,"# a0:      %10.7Lf , ax:      %10.7Lf \n", a0min, a0max); 
  fprintf(fl1,"# astep:   %10.7Lf , ystep:   %10.7Lf \n", a0step, ystep); 
  fprintf(fl1,"# ymin:    %10.7Lf , ymax:    %10.7Lf \n", ymin, ymax); 
  fprintf(fl1,"# b0:      %10.7Lf , bx:      %10.7Lf \n", b0, bxx); 
  fprintf(fl1,"# bres:    %10d , xbranch: %10d\n", bres, xbranch); 
  fprintf(fl1,"# inters:  %10.0Lf , bsharp:  %10.7Lf \n", INTERVALS, bsharp); 
  fprintf(fl1,"# ydist:   %10.3Le , ydistmin:%10.3Le\n", bracdist, bracdistmin); 
  fprintf(fl1,"# mod:     %10lu \n", itstime2mod); 
  fprintf(fl1,"# %23s%25s%25s%15s\n\n", "a", "b", "y", "# sol");
}

void maincalc() 
{
/*    long double braca, bracb, bracc, sol[10][2], nsol[10][2], sol[2*xbranch+2][2], nsol[2*xbranch+2][2], checker, bracd, mindist, intervals=INTERVALS; */

    long double braca, bracb, bracc, sol[IPOW3(xbranch)+3][2], nsol[IPOW3(xbranch)+3][2], 
	checker, inters=INTERVALS, boundl, boundu, tmp[2]; 
    int xcount, solno=IPOW3(xbranch), i, j, k, jold, okbrac, esti[2];

    /* find initial extrema */
    a0=a0min;
    if (fabsl(a0min-omega)<sqrt(LDBL_EPSILON)) {
	xcount=1;
	sol[xcount][0]=0.0L; sol[xcount][1]=0.0L;
	nsol[xcount][0]=0.0L; nsol[xcount][1]=0.0L;
	boundl=0.0L-bracdist; boundu=0.0L+bracdist;
    } else {	
	
	xcount=1;
	bracb=sqrtl(1.0L-omega/a0min); /* initial guess */
	if (rootbyfunc(bracb, sln, qq[qpi], pp[qpi])>=0.0L) { bracb=-bracb;
	if (rootbyfunc(bracb, sln, qq[qpi], pp[qpi])>=0.0L) { bracb=0.0L;
	if (rootbyfunc(bracb, sln, qq[qpi], pp[qpi])>=0.0L) {
	    printf("(%d) Warning: No valid initial guess found, aborting...\n", sln);
	    exit(-10); }}}
	
	/* estimate bounds on domain */
	boundl=bracb-1.0L/bracdist; boundu=bracb+1.0L/bracdist;
	tmp[0]=boundl; tmp[1]=bracb; nsol[0][1]=-1.0L;
	if (rootbyfunc(boundl, sln, qq[qpi], pp[qpi])>=0.0L) {
	  esti[0]=1; while (fabsl(tmp[0]-tmp[1])>bsharp && esti[0]<INTERVALS)
	    {
	      esti[0]++; boundl=(tmp[0]+tmp[1])/2.0L;
	      if ((braca=rootbyfunc(boundl, sln, qq[qpi], pp[qpi]))>=0.0L) tmp[0]=boundl;
	      else {tmp[1]=boundl; nsol[0][1]=-braca;}
	    }
	    boundl=tmp[1]; nsol[0][0]=boundl;
	    if (nsol[0][1]<0.0L) nsol[0][1]=-rootbyfunc(boundl, sln, qq[qpi], pp[qpi]);
	}
	
	tmp[0]=bracb; tmp[1]=boundu; nsol[xcount+1][1]=-1.0L;
	if (rootbyfunc(boundu, sln, qq[qpi], pp[qpi])>=0.0L) {
	  esti[1]=1; while (fabsl(tmp[0]-tmp[1])>bsharp && esti[1]<INTERVALS)
	    {
	      esti[1]++; boundu=(tmp[0]+tmp[1])/2.0L;
		if ((bracc=rootbyfunc(boundu, sln, qq[qpi], pp[qpi]))>=0.0L) tmp[1]=boundu;
		else {tmp[0]=boundu; nsol[xcount+1][1]=-bracc;}
	    }
	    boundu=tmp[0]; nsol[xcount+1][0]=boundu;
	    if (nsol[xcount+1][1]<0.0L) nsol[xcount+1][1]=-rootbyfunc(boundu, sln, qq[qpi], pp[qpi]);
	}

	/* find first maximum */
	braca=boundl; bracc=boundu; bracb=(boundl+boundu)/2.0L;
	if((okbrac=mymnbrac(&braca, &bracb, &bracc, inters, rootbyfunc, sln, qq[qpi], pp[qpi]))==1)
	    nsol[xcount][1]=-brent(braca, bracb, bracc, rootbyfunc, sln, qq[qpi], pp[qpi], 
				   &nsol[xcount][0]);
	else { printf("(%d) Warning: No start maximum found, aborting...\n", sln); exit(-1); }

	if (rootbyfunc(nsol[xcount][0], sln, qq[qpi], pp[qpi])>=0.0L) {
	  printf("(%d) Warning: No start maximum found, aborting...\n", sln);
	  exit(-2); }

	/* find higher branch maxima */
	i=1; while (i<=xbranch && xcount+2<=solno){
	    j=1; while (j<=xcount+2 && xcount+2<=solno){
		braca=(j>=3? nsol[j-2][0]:boundl); 
		bracc=(j<=xcount? nsol[j][0]: boundu);
		/* bracb=(bracc-braca)/inters/10.0L; /* correction factor to exclude bounds */
		/* braca+=bracb; bracc-=bracb; */
		bracb=braca-1.0L; 
		if((okbrac=mymnbrac(&braca, &bracb, &bracc, inters, rootbyfunc, sln, qq[qpi], pp[qpi]))==1)
		    tmp[1]=-brent(braca, bracb, bracc, rootbyfunc, sln, qq[qpi], pp[qpi], &tmp[0]);
		
		if (okbrac==1 && rootbyfunc(tmp[0], sln, qq[qpi], pp[qpi])<0.0L) {
		    for (k=xcount+1; k>=j; k--) {nsol[k+2][0]=nsol[k][0]; nsol[k+2][1]=nsol[k][1];}
		    nsol[j][0]=tmp[0]; nsol[j][1]=tmp[1]; xcount+=2;
		    if (nsol[j-1][0]<nsol[j][0] && j+2<=xcount) {
			braca=nsol[j][0]; bracb=braca-1.0L; bracc=nsol[j+2][0];
			if((okbrac=mymnbrac(&braca, &bracb, &bracc, inters, rootbynegfunc, sln, qq[qpi], pp[qpi]))==1)
			  nsol[j+1][1]=brent(braca, bracb, bracc, rootbynegfunc, sln, qq[qpi], pp[qpi], 
					     &nsol[j+1][0]);
			else printf("(%d) Warning: No right minimum between start maxima %d...\n", sln, j+1);
		    }
		    else if (nsol[j-1][0]>=nsol[j][0] && j>=3) {
		        nsol[j+1][0]=nsol[j-1][0]; nsol[j+1][1]=nsol[j-1][1];
			braca=nsol[j-2][0]; bracb=braca-1.0L; bracc=nsol[j][0];
			if((okbrac=mymnbrac(&braca, &bracb, &bracc, inters, rootbynegfunc, sln, qq[qpi], pp[qpi]))==1)
			  nsol[j-1][1]=brent(braca, bracb, bracc, rootbynegfunc, sln, qq[qpi], pp[qpi], 
					   &nsol[j-1][0]);
			else printf("(%d) Warning: No left minimum between start maxima %d...\n", sln, j-1);
		    }
		    else printf("(%d) Warning: No minimum between start maxima %d...\n", sln , j);
		    
		    j+=2;
		}
	    j+=2;}
	i++;}
    } /* end else */
    
    for (i=0; i<=xcount+1; i++) {
	fprintf(fl1, "%25.18Lf%25.18Lf%25.18Lf%15d\n", 
		a0, nsol[i][1], nsol[i][0], (i<1 || i>xcount ? esti[i/(xcount+1)]*100+xcount: 2*(i%2)-1));  
	sol[i][0]=nsol[i][0]; sol[i][1]=nsol[i][1];
    }
    fprintf(fl1,"\n");
    

    /* find solutions at a0>a0min */
    for (a0=a0min+a0step; a0<=a0max; a0+=a0step) {

      /* reestimate lower bound on domain */
      bracb=bracdistmin;
      tmp[0]=boundl-bracb; tmp[1]=boundl+bracb; okbrac=0;
      esti[0]=1; while (okbrac==0 && bracb<=bracdist && esti[0]<INTERVALS) {
	esti[0]++; if (rootbyfunc(tmp[0], sln, qq[qpi], pp[qpi])>=0.0L) { 
	  if (rootbyfunc(tmp[1], sln, qq[qpi], pp[qpi])<0.0L) okbrac=1;
	  else {
	    bracb*=5.0L;
	    tmp[0]=tmp[1]; tmp[1]+=bracb;	    
	  }
	}
	else {
	  bracb*=5.0L;
	  tmp[1]=tmp[0]; tmp[0]-=bracb;
	}
      }
      if (!okbrac) printf("(%d) Warning, lower bound lost at a=%Lf.\n", sln, a0);
      nsol[0][1]=-1.0L;
      esti[0]*=1000; while (fabsl(tmp[0]-tmp[1])>bsharp && esti[0]%1000<INTERVALS)
	{
	  esti[0]++; boundl=(tmp[0]+tmp[1])/2.0L;
	  if ((braca=rootbyfunc(boundl, sln, qq[qpi], pp[qpi]))>=0.0L) tmp[0]=boundl;
	  else {tmp[1]=boundl; nsol[0][1]=-braca;}
	}
      boundl=tmp[1]; nsol[0][0]=boundl;
      if (nsol[0][1]<0.0L) nsol[0][1]=-rootbyfunc(boundl, sln, qq[qpi], pp[qpi]);
      

      /* reestimate upper bound on domain */
      bracb=bracdistmin;
      tmp[0]=boundu-bracb; tmp[1]=boundu+bracb; okbrac=0;
      esti[1]=1; while (okbrac==0 && bracb<=bracdist && esti[1]<INTERVALS) {
	esti[1]++; if (rootbyfunc(tmp[1], sln, qq[qpi], pp[qpi])>=0.0L) { 
	  if (rootbyfunc(tmp[0], sln, qq[qpi], pp[qpi])<0.0L) okbrac=1;
	  else {
	    bracb*=5.0L;
	    tmp[1]=tmp[0]; tmp[0]-=bracb;	    
	  }
	}
	else {
	  bracb*=5.0L;
	  tmp[0]=tmp[1]; tmp[1]+=bracb;
	}
      }
      if (!okbrac) printf("(%d) Warning, upper bound lost at a=%Lf.\n", sln, a0);
      nsol[xcount+1][1]=-1.0L;
      esti[1]*=1000; while (fabsl(tmp[0]-tmp[1])>bsharp && esti[1]%1000<INTERVALS)
	{
	  esti[1]++; boundu=(tmp[0]+tmp[1])/2.0L;
	  if ((bracc=rootbyfunc(boundu, sln, qq[qpi], pp[qpi]))>=0.0L) tmp[1]=boundu;
	  else {tmp[0]=boundu; nsol[xcount+1][1]=-bracc;}
	}
      boundu=tmp[0]; nsol[xcount+1][0]=boundu;
      if (nsol[xcount+1][1]<0.0L) nsol[xcount+1][1]=-rootbyfunc(boundu, sln, qq[qpi], pp[qpi]);

      /* find new maxima from old ones */
      i=1; j=1; while (j<xcount && xcount+2<=solno){
	braca=nsol[j-1][0]; 
	bracb=sol[j][0]; 
	bracc=sol[j+1][0];
	if((okbrac=mymnbrac(&braca, &bracb, &bracc, inters, rootbyfunc, sln, qq[qpi], pp[qpi]))==1)
	  tmp[1]=-brent(braca, bracb, bracc, rootbyfunc, sln, qq[qpi], pp[qpi], &tmp[0]);
	
	if (okbrac==1 && rootbyfunc(tmp[0], sln, qq[qpi], pp[qpi])<0.0L) {
	  nsol[j][0]=tmp[0]; nsol[j][1]=tmp[1];
	  if (nsol[j-1][0]<nsol[j][0] && j+2<=xcount) {
	    braca=nsol[j][0]; bracb=sol[j+1][0]; bracc=sol[j+2][0];
	    if((okbrac=mymnbrac(&braca, &bracb, &bracc, inters, rootbynegfunc, sln, qq[qpi], pp[qpi]))==1)
	      nsol[j+1][1]=brent(braca, bracb, bracc, rootbynegfunc, sln, qq[qpi], pp[qpi], 
				 &nsol[j+1][0]);
	    else printf("(%d) Warning: No right minimum between maxima %d at a0=%Lf...\n", sln, j+1, a0);
	  }
	  else if (nsol[j-1][0]>=nsol[j][0] && j>=3) {
	    nsol[j+1][0]=nsol[j-1][0]; nsol[j+1][1]=nsol[j-1][1];
	    braca=nsol[j-2][0]; bracb=sol[j-1][0]; bracc=sol[j][0];
	    if((okbrac=mymnbrac(&braca, &bracb, &bracc, inters, rootbynegfunc, sln, qq[qpi], pp[qpi]))==1)
	      nsol[j-1][1]=brent(braca, bracb, bracc, rootbynegfunc, sln, qq[qpi], pp[qpi], 
				 &nsol[j-1][0]);
	    else printf("(%d) Warning: No left minimum between maxima %d at a0=%Lf...\n", sln, j-1, a0);
	  }
	  else printf("(%d) Warning: No minimum between maxima %d... at a0=%Lf\n", sln , j, a0);
	  
	}
	j+=2;}
      if (j==xcount) {	      
	  braca=nsol[j-1][0]; bracb=sol[j][0]; bracc=sol[j+1][0];
	  if((okbrac=mymnbrac(&braca, &bracb, &bracc, inters, rootbyfunc, sln, qq[qpi], pp[qpi]))==1)
	      tmp[1]=-brent(braca, bracb, bracc, rootbyfunc, sln, qq[qpi], pp[qpi], &tmp[0]);
	  if (okbrac==1 && rootbyfunc(tmp[0], sln, qq[qpi], pp[qpi])<0.0L) 
	  { nsol[j][0]=tmp[0]; nsol[j][1]=tmp[1];}
      }

      /* find new maxima branches in between */
      i++; while (i<=xbranch && xcount+2<=solno){
	  j=1; while (j<=xcount+2 && xcount+2<=solno){
	      braca=(j>=3? nsol[j-2][0]:boundl); 
	      bracc=(j<=xcount? nsol[j][0]: boundu);
	      bracb=braca-1.0L; 
	      if((okbrac=mymnbrac(&braca, &bracb, &bracc, inters, rootbyfunc, sln, qq[qpi], pp[qpi]))==1)
		  tmp[1]=-brent(braca, bracb, bracc, rootbyfunc, sln, qq[qpi], pp[qpi], &tmp[0]);
	      
	      if (okbrac==1 && rootbyfunc(tmp[0], sln, qq[qpi], pp[qpi])<0.0L) {
		  for (k=xcount+1; k>=j; k--) {nsol[k+2][0]=nsol[k][0]; nsol[k+2][1]=nsol[k][1];}
		  nsol[j][0]=tmp[0]; nsol[j][1]=tmp[1]; xcount+=2;
		  if (nsol[j-1][0]<nsol[j][0] && j+2<=xcount) {
		      braca=nsol[j][0]; bracb=braca-1.0L; bracc=nsol[j+2][0];
		      if((okbrac=mymnbrac(&braca, &bracb, &bracc, inters, rootbynegfunc, sln, qq[qpi], pp[qpi]))==1)
			  nsol[j+1][1]=brent(braca, bracb, bracc, rootbynegfunc, sln, qq[qpi], pp[qpi], 
					     &nsol[j+1][0]);
		      else printf("(%d) Warning: No right minimum between maxima %d at a0=%Lf...\n", sln, j+1, a0);
		  }
		  else if (nsol[j-1][0]>=nsol[j][0] && j>=3) {
		      nsol[j+1][0]=nsol[j-1][0]; nsol[j+1][1]=nsol[j-1][1];
		      braca=nsol[j-2][0]; bracb=braca-1.0L; bracc=nsol[j][0];
		      if((okbrac=mymnbrac(&braca, &bracb, &bracc, inters, rootbynegfunc, sln, qq[qpi], pp[qpi]))==1)
			  nsol[j-1][1]=brent(braca, bracb, bracc, rootbynegfunc, sln, qq[qpi], pp[qpi], 
					     &nsol[j-1][0]);
		      else printf("(%d) Warning: No left minimum between maxima %d at a0=%Lf...\n", sln, j-1, a0);
		  }
		  else printf("(%d) Warning: No minimum between maxima %d at a0=%Lf...\n", sln , j, a0);
		  
		  j+=2;
	      }
	      j+=2;}
	  i++;}
      
      for (i=0; i<=xcount+1; i++) {
	fprintf(fl1, "%25.18Lf%25.18Lf%25.18Lf%15d\n", 
		a0, nsol[i][1], nsol[i][0], (i<1 || i>xcount ? esti[i/(xcount+1)]*100+xcount: 2*(i%2)-1));   
	sol[i][0]=nsol[i][0]; sol[i][1]=nsol[i][1];
      }
      fprintf(fl1,"\n");
      
    }
    
    /* find first solutions at a0=a0min
    mindist=0.0L; a0=a0min;
    if (fabsl(a0min-omega)<sqrt(LDBL_EPSILON)) {
	xcount=0;
	sol[0][0]=0.0L; sol[0][1]=0.0L;
	nsol[0][0]=0.0L; nsol[0][1]=0.0L;
    } else {	
	
	xcount=0;
	braca=sqrtl(1.0L-omega/a0min);
	if (rootbyfunc(braca, sln, qq[qpi], pp[qpi])>0.0L) braca=-braca;
	bracd=bracdist;
	bracb=braca-bracd; 
	mnbrac(&braca, &bracb, &bracc, rootbyfunc, sln, qq[qpi], pp[qpi]);
	nsol[xcount][1]=-brent(braca, bracb, bracc, rootbyfunc, sln, qq[qpi], pp[qpi], 
			      &nsol[xcount][0]);
	
	(xcount)++;
	braca=nsol[0][0]; bracb=braca+bracd;
	if ((checker=rootbyfunc(bracb, sln, qq[qpi], pp[qpi]))<0.0L) {
	    mnbrac(&braca, &bracb, &bracc, rootbynegfunc, sln, qq[qpi], pp[qpi]);
	    nsol[xcount][1]=brent(braca, bracb, bracc, rootbynegfunc, sln, qq[qpi], pp[qpi], 
				 &nsol[xcount][0]);
	    checker=rootbyfunc(nsol[xcount][0]+fabs(bracd*(braca-bracb))/(bracb-braca), 
			       sln, qq[qpi], pp[qpi]);
	}
	if (checker>0.0L) (xcount)--;
      
	while (checker <0.0L && (xcount)<2*xbranch) {
	    (xcount)++;
	    braca=nsol[xcount-1][0]; bracb=braca+bracd;
	    mnbrac(&braca, &bracb, &bracc, rootbyfunc, sln, qq[qpi], pp[qpi]);
	    nsol[xcount][1]=-brent(braca, bracb, bracc, rootbyfunc, sln, qq[qpi], pp[qpi], 
				  &nsol[xcount][0]);
	    
	    (xcount)++;
	    braca=nsol[xcount-1][0]; bracb=braca+bracd;
	    if ((checker=rootbyfunc(bracb, sln, qq[qpi], pp[qpi]))<0.0L) {
		mnbrac(&braca, &bracb, &bracc, rootbynegfunc, sln, qq[qpi], pp[qpi]);
		nsol[xcount][1]=brent(braca, bracb, bracc, rootbynegfunc, sln, qq[qpi], pp[qpi], 
				     &nsol[xcount][0]);
		checker=rootbyfunc(nsol[xcount][0]+fabs(bracd*(braca-bracb))/(bracb-braca), 
				   sln, qq[qpi], pp[qpi]);
	    }
	    if (checker>0.0L) (xcount)--;
	}    
	
	if ((xcount)<2*xbranch) {
	    (xcount)++;
	    braca=nsol[0][0]; bracb=braca-bracd;
	    mnbrac(&braca, &bracb, &bracc, rootbynegfunc, sln, qq[qpi], pp[qpi]);
	    nsol[xcount][1]=brent(braca, bracb, bracc, rootbynegfunc, sln, qq[qpi], pp[qpi], 
				 &nsol[xcount][0]);
	    checker=rootbyfunc(nsol[xcount][0]+fabs(bracd*(braca-bracb))/(bracb-braca), 
			       sln, qq[qpi], pp[qpi]);
	    if (checker>0.0L) (xcount)--;
	    
	    while (checker <0.0L && (xcount)<2*xbranch) {
		(xcount)++;
		braca=nsol[xcount-1][0]; bracb=braca-bracd;
		mnbrac(&braca, &bracb, &bracc, rootbyfunc, sln, qq[qpi], pp[qpi]);
		nsol[xcount][1]=-brent(braca, bracb, bracc, rootbyfunc, sln, qq[qpi], pp[qpi], 
				      &nsol[xcount][0]);
		
		(xcount)++;
		braca=nsol[xcount-1][0]; bracb=braca-bracd;
		if ((checker=rootbyfunc(bracb, sln, qq[qpi], pp[qpi]))<0.0L) {
		    mnbrac(&braca, &bracb, &bracc, rootbynegfunc, sln, qq[qpi], pp[qpi]);
		    nsol[xcount][1]=brent(braca, bracb, bracc, rootbynegfunc, sln, qq[qpi], pp[qpi], 
					 &nsol[xcount][0]);
		    checker=rootbyfunc(nsol[xcount][0]+fabs(bracd*(braca-bracb))/(bracb-braca), 
				       sln, qq[qpi], pp[qpi]);
		}
		if (checker>0.0L) (xcount)--;
	    }    
	}
	
	qsort(nsol, xcount+1, 2*sizeof(long double), qsortcompare);
	j=0; for (i=1; i<=xcount; i++) {
	    if (nsol[j][0]-nsol[i][0]>2.5L*sqrtl(LDBL_EPSILON)) {
		nsol[++j][0]=nsol[i][0];
		nsol[j][1]=nsol[i][1];
		mindist=(mindist<nsol[j-1][0]-nsol[j][0] ?
			 mindist : nsol[j-1][0]-nsol[j][0]);
	    }	  
	} 
	if (xcount!=j) printf("Warning, found %d times same roots twice at a=%Lf.\n", xcount-j, a0min); 
	xcount=j;
	
    }
    for (i=0; i<=xcount; i++) {
	fprintf(fl1, "%25.18Lf%25.18Lf%25.18Lf%10d\n", a0min, nsol[i][1], nsol[i][0], xcount+1);  
	sol[i][0]=nsol[i][0]; sol[i][1]=nsol[i][1];
    }	
    if (mindist<2.5L*LDBL_EPSILON) mindist=1.0L/bracdist;
    /* found first solutions at a0=a0min

    /* find further solutions close to the previous ones 
    for (a0=a0min+a0step; a0<=a0max; a0+=a0step) {
	bracd=(mindist/3.0L<bracdist ? mindist/3.0L : bracdist);
	bracd=(bracd>bracdistmin ? bracd : bracdistmin);
	j=0; jold=0; for (i=0;i<=xcount;i+=2) {

	    /* search above previous maximum
	    braca=sol[i][0]+bracd; bracb=braca+bracd/10.0L;
	    mnbrac(&braca, &bracb, &bracc, rootbyfunc, sln, qq[qpi], pp[qpi]);
	    nsol[j][1]=-brent(braca, bracb, bracc, rootbyfunc, sln, qq[qpi], pp[qpi], 
				   &nsol[j][0]);
	    if (fabsl(nsol[j][1]+ERROR)>SMALL) 
		if (nsol[j][0]<sol[i][0]+mindist/2.0L && nsol[j][0]>sol[i][0]-mindist/2.0L) 
		  j+=2;
	    
	    /* search below previous maximum
	    braca=sol[i][0]-bracd; bracb=braca-bracd/10.0L;
	    mnbrac(&braca, &bracb, &bracc, rootbyfunc, sln, qq[qpi], pp[qpi]);
	    nsol[j][1]=-brent(braca, bracb, bracc, rootbyfunc, sln, qq[qpi], pp[qpi], 
				   &nsol[j][0]);
	    if (fabsl(nsol[j][1]+ERROR)>SMALL) 
		if (nsol[j][0]<sol[i][0]+mindist/2.0L && nsol[j][0]>sol[i][0]-mindist/2.0L) 
		  j++;

	    /* find intermediate minimum if applicable
	    if (j==jold+2) j-=1;
	    else if (j==jold+3) {
  	        if (nsol[jold][0]-nsol[jold+2][0]>bracdistmin) {
		  checker=rootbynegfunc((nsol[jold][0]+nsol[jold+2][0])/2.0L, 
					sln, qq[qpi], pp[qpi]);
		  if (rootbynegfunc(nsol[jold][0], sln, qq[qpi], pp[qpi])>checker &&
		      rootbynegfunc(nsol[jold+2][0], sln, qq[qpi], pp[qpi])>checker)
		    nsol[j-2][1]=brent(nsol[jold][0], (nsol[jold][0]+nsol[jold+2][0])/2.0L,
				       nsol[jold+2][0], rootbynegfunc, sln, qq[qpi], pp[qpi], 
				       &nsol[j-2][0]);
		  else {
		    braca=sol[jold+2][0]; bracb=braca-(nsol[jold][0]-braca)/10.0L;
		    mnbrac(&braca, &bracb, &bracc, rootbynegfunc, sln, qq[qpi], pp[qpi]);
		    nsol[j-2][1]=brent(braca, bracb, bracc, rootbynegfunc, sln, qq[qpi], pp[qpi], 
				       &nsol[j-2][0]);
		    if (nsol[j-2][0]<nsol[jold+2][0] || nsol[j-2][0]>nsol[jold][0]
			|| fabsl(nsol[j-2][1]-ERROR)<SMALL) 
		      {nsol[j-2][0]=nsol[jold+2][0]; nsol[j-2][1]=nsol[jold+2][1];
		      printf("Uhoh... No min. between 2 max.es at a0=%Lf?!\n", a0);}
		  }
		}
		else j-=2;
		
	    }
	    jold=j;
	    
	    if (i+1<=xcount) 
		{
		    /* search above previous minimum
		    braca=sol[i+1][0]+bracd; bracb=braca+bracd/10.0L;
		    mnbrac(&braca, &bracb, &bracc, rootbynegfunc, sln, qq[qpi], pp[qpi]);
		    nsol[j][1]=brent(braca, bracb, bracc, rootbynegfunc, sln, qq[qpi], pp[qpi], 
				     &nsol[j][0]);
		    checker=rootbyfunc(nsol[j][0]+fabs(bracd*(braca-bracb))/(bracb-braca), 
				       sln, qq[qpi], pp[qpi]);
		    if (checker>0.0L) nsol[j][1]=ERROR;
		    if (fabsl(nsol[j][1]-ERROR)>SMALL) 
		      if (nsol[j][0]<sol[i+1][0]+mindist/2.0L && nsol[j][0]>sol[i+1][0]-mindist/2.0L) 
			j+=2;
		    
		    /* search below previous minimum
		    braca=sol[i+1][0]-bracd; bracb=braca+bracd/10.0L;
		    mnbrac(&braca, &bracb, &bracc, rootbynegfunc, sln, qq[qpi], pp[qpi]);
		    nsol[j][1]=brent(braca, bracb, bracc, rootbynegfunc, sln, qq[qpi], pp[qpi], 
				     &nsol[j][0]);
		    checker=rootbyfunc(nsol[j][0]+fabs(bracd*(braca-bracb))/(bracb-braca), 
				       sln, qq[qpi], pp[qpi]);
		    if (checker>0.0L) nsol[j][1]=ERROR;
		    if (fabsl(nsol[j][1]-ERROR)>SMALL) 
		      if (nsol[j][0]<sol[i+1][0]+mindist/2.0L && nsol[j][0]>sol[i+1][0]-mindist/2.0L) 
			j++;
		    
		    /* find intermediate maximum if applicable
		    if (j==jold+2) j-=1;
		    else if (j==jold+3) {
		        if (nsol[jold][0]-nsol[jold+2][0]>bracdistmin) {
			  checker=rootbynegfunc((nsol[jold][0]+nsol[jold+2][0])/2.0L, 
						sln, qq[qpi], pp[qpi]);
			  if (rootbyfunc(nsol[jold][0], sln, qq[qpi], pp[qpi])>checker &&
			      rootbyfunc(nsol[jold+2][0], sln, qq[qpi], pp[qpi])>checker)
			    nsol[j-2][1]=-brent(nsol[jold][0], (nsol[jold][0]+nsol[jold+2][0])/2.0L,
						nsol[jold+2][0], rootbyfunc, sln, qq[qpi], pp[qpi], 
						&nsol[j-2][0]);
			  else {
			    braca=sol[jold+2][0]; bracb=braca+(nsol[jold][0]-braca)/10.0L;
			    mnbrac(&braca, &bracb, &bracc, rootbyfunc, sln, qq[qpi], pp[qpi]);
			    nsol[j-2][1]=-brent(braca, bracb, bracc, rootbyfunc, sln, qq[qpi], pp[qpi], 
						&nsol[j-2][0]);
			    if (nsol[j-2][0]<nsol[jold+2][0] || nsol[j-2][0]>nsol[jold][0]
				|| fabsl(nsol[j-2][1]+ERROR)<SMALL) 
			      {nsol[j-2][0]=nsol[jold+2][0]; nsol[j-2][1]=nsol[jold+2][1];
			      printf("Uhoh... No min. between 2 max.es at a0=%Lf?!\n", a0);}
			  }
			} else j-=2;
		    }
		    jold=j;
		}
	}
	j-=1;

	xcount=j;
	fprintf(fl1, "%25.18Lf%25.18Lf%25.18Lf%10d\n", a0, nsol[0][1], nsol[0][0], xcount+1);  
	if (xcount>=0) {sol[0][0]=nsol[0][0]; sol[0][1]=nsol[0][1];}
	for (i=1; i<=xcount; i++) {
	    fprintf(fl1, "%25.18Lf%25.18Lf%25.18Lf%10d\n", a0, nsol[i][1], nsol[i][0], xcount+1);  
	    sol[i][0]=nsol[i][0]; sol[i][1]=nsol[i][1];
	    if (nsol[i-1][0]-nsol[i][0]>bracdistmin)	
	      mindist=(mindist<nsol[i-1][0]-nsol[i][0] ? mindist : nsol[i-1][0]-nsol[i][0]);	
	}
	}*/   
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
    long double EPS=LDBL_EPSILON, ztol=LDBL_EPSILON;
    int ITMAX=100;

    int i;
    long double a=x1,b=x2,c=x2,d,e,min1,min2;
    long double fa, fb, fc,p,q,r,s,tol1,xm;

    fa = (*func)(a, par0, par1, par2, par3); fb = (*func)(b, par0, par1, par2, par3);
    if ((fa > 0.0L && fb > 0.0L) || (fa < 0.0L && fb < 0.0L)) return ERROR;
    fc=fb;
    for (i=1;i<=ITMAX;i++) {
	if ((fb > 0.0L && fc > 0.0L) || (fb < 0.0L && fc < 0.0L)) {
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
	tol1=2.0L*EPS*fabsl(b)+0.5L*ztol;
	xm=0.5L*(c-b);
	if (fabsl(xm) <= tol1 || fb == 0.0L) return b;
	if (fabsl(e) >= tol1 && fabsl(fa) > fabsl(fb)) {
	    s=fb/fa;
	    if (a == c) {
		p=2.0L*xm*s;
		q=1.0L-s;
	    } else {
		q=fa/fc;
		r=fb/fc;
		p=s*(2.0L*xm*q*(q-r)-(b-a)*(r-1.0L));
		q=(q-1.0L)*(r-1.0L)*(s-1.0L);
	    }
	    if (p > 0.0L) q = -q;
	    p=fabsl(p);
	    min1=3.0L*xm*q-fabsl(tol1*q);
	    min2=fabsl(e*q);
	    if (2.0L*p < (min1 <min2 ? min1 : min2)) {
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
	    b += (xm > 0.0L ? fabsl(tol1) : -fabsl(tol1));
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

int mymnbrac(long double *ax, long double *bx, long double *cx, long double distf,
	    long double (*func)(long double, int, unsigned long int, unsigned long int), 
		   int par1, unsigned long int par2, unsigned long int par3)
/* find minimum bracket between *ax and *cx, i.e. f(*bx)<f(*ax),
f(*bx)<f(*cx), where *ax and *cx must not be exceeded */
{
  long double GLIMIT=5.0L, EPS=LDBL_EPSILON, TINY=EPS;
  long double GOLD=1.6180339887498948482045868343656381177203091798058L;
  int NTRY=10;
 
  long double ulim, u, r, q, fu, dum;
  long double fa, fb, fc, a0=*ax, b0=*bx, c0=*cx, distx;
  int i;

  if (fabsl(*ax-*cx)<5.0L*TINY) return 0;          /* bracket size too small */ 
  distx=(*cx-*ax)/distf; if (fabsl(distx)<5.0L*TINY) distx=5.0L*TINY;
  if ((*ax-*bx)*(*bx-*cx)<0) *bx=(*ax + *cx)/2.0L; b0=*bx; /* need *bx inside bracket */
    
  fa = (*func)(a0, par1, par2, par3); 
  fb = (*func)(b0, par1, par2, par3);
  fc = (*func)(c0, par1, par2, par3);
  if (fa>fb && fc>fb) return 1;                    /* bracket already found */

  *ax=b0; *bx=b0+distx/2.0L; *cx=b0+distx;         /* try search from center of bracket */
  fa = (*func)(*ax, par1, par2, par3); 
  fb = (*func)(*bx, par1, par2, par3);
  fc = (*func)(*cx, par1, par2, par3);
  if (fa>fb && fc>fb) return 1;                    /* bracket inside start interval */

  *bx=*cx; fb=fc;
  if (fb > fa) {
    SHFT(dum, *ax, *bx, dum);
    SHFT(dum, fb, fa, dum);
    SHFT(dum, a0, c0, dum);
    distx=-distx;
  }
  
  *cx = (*bx) + GOLD*(*bx-*ax);
  fc = (*func)(*cx, par1, par2, par3);

  for (i=1;i<=NTRY;i++) {
      while (fb > fc && (c0-*cx)*(*cx-*bx)>0.0L) {
	  r = (*bx - *ax)*(fb-fc);
	  q = (*bx - *cx)*(fb-fa);
	  u = (*bx) - 
	      ((*bx - *cx)*q - (*bx - *ax)*r)/(2.0*SIGN(FMAX(fabsl(q-r), TINY),q-r));
	  ulim = (*bx) + GLIMIT*(*cx-*bx);               /* allowed limit for interpolation */
	  if ((c0-ulim)*(*cx-*bx)<0.0L) ulim=c0-distx;
	  if ((*bx-u)*(u-*cx) > 0.0) {                   /* parabolic u between *bx and *cx */
	      fu = (*func)(u, par1, par2, par3);
	      if (fu < fc) {
		  *ax = (*bx);
		  *bx = u;
		  fa = fb;
		  return 1;
	      } else if (fu > fb) {
		  *cx = u;
		  fc = fu;
		  return 1;
	      }
	      u = (*cx) + GOLD*(*cx-*bx);
	      fu = (*func)(u, par1, par2, par3);
	  } else if ((*cx-u)*(u-ulim) > 0.0) {           /* parabolic u between *cx and limit */
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
      if (fb < fc && (c0-*cx)*(*cx-*bx)>0.0L) return 1;
      if ((c0-*cx)*(*cx-*bx)<=0.0L) {
	  *cx=c0; 
	  fc = (*func)(*cx, par1, par2, par3);
	  if (fb < fc) return 1;
      }
      
      c0=b0; b0=(a0+c0)/2.0L;
      if (i==NTRY) {                                   /* last try */ 
	  *ax=a0; *bx=a0+distx/2.0L; *cx=a0+distx;
	  if ((*cx-c0)*(*cx-*bx)>0.0) {
	      *cx=c0; *bx=(*ax+*cx)/2.0L;
	  }
	  fa = (*func)(*ax, par1, par2, par3); 
	  fb = (*func)(*bx, par1, par2, par3);
	  fc = (*func)(*cx, par1, par2, par3);
	  if (fa>fb && fc>fb) return 1;                /* bracket inside start interval */
      
	  *bx=*cx; fb=fc;
	  if (fb>fa) return -i;                        /* wrong bracket direction */
      }
      
      if (i<NTRY) {                                    /* try again, in ignored end of bracket */ 
	  *ax=b0; *bx=*ax+distx/2.0L; *cx=*ax+distx;
	  if ((*cx-c0)*(*cx-*bx)>0.0) {
	      i=NTRY;
	      *cx=c0; *bx=(*ax+*cx)/2.0L;
	  }
	  fa = (*func)(*ax, par1, par2, par3); 
	  fb = (*func)(*bx, par1, par2, par3);
	  fc = (*func)(*cx, par1, par2, par3);
	  if (fa>fb && fc>fb) return 1;                /* bracket inside start interval */
	  *bx=*cx; fb=fc;
	  if (fb > fa) {
	      SHFT(dum, *ax, *bx, dum);
	      SHFT(dum, fb, fa, dum);
	      SHFT(dum, a0, c0, dum);
	      distx=-distx;
	  }
      }
  
      *cx = (*bx) + GOLD*(*bx-*ax);
      fc = (*func)(*cx, par1, par2, par3);
  }
  if (fb < fc && (c0-*cx)*(*cx-*bx)>0.0L) return 1;
  else return -1;
}

long double brent(long double ax, long double bx, long double cx, 
	    long double (*f)(long double, int, unsigned long int, unsigned long int), 
		   int par1, unsigned long int par2, unsigned long int par3, long double *xmin)
{
  long double ITMAX=100, ZEPS=LDBL_EPSILON, tol=sqrtl(ZEPS);
  long double CGOLD=0.3819660112501051517954131656343618822796908201942L;
  int iter;
  long double a,b,d,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm,e=0.0L;

  a=(ax < cx ? ax : cx);
  b=(ax > cx ? ax : cx);
  x=w=v=bx;
  fw=fv=fx=(*f)(x, par1, par2, par3); 
  for (iter=1; iter<=ITMAX; iter++) {
    xm=0.5L*(a+b);
    tol2= 2.0L*(tol1=tol*fabsl(x)+ZEPS);
    if (fabsl(x-xm) <= (tol2-0.5L*(b-a))) {
      *xmin = x;
      return fx;
    }
    if (fabsl(e) > tol1) {
      r = (x-w)*(fx-fv);
      q = (x-v)*(fx-fw);
      p = (x-v)*q - (x-w)*r;
      q = 2.0L*(q-r);
      if (q > 0.0L) p = -p;
      q = fabsl(q);
      etemp=e;
      e=d;
      if (fabsl(p) >= fabsl(0.5L*q*etemp) || p <= q*(a-x) || p >= q*(b-x))
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

int qsortcompare(const void *a, const void *b)
{
  return (((long double *)a)[0]<((long double *)b)[0]? 
	  1 : (((long double *)a)[0]>((long double *)b)[0] ? -1 : 0));
}
