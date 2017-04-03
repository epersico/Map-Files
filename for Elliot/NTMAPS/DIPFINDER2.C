#include <stdio.h>
#include <math.h>
#include <string.h>
#define NL(flnm,dumch) {do dumch=fgetc(flnm); while (dumch!='\n' && dumch!=EOF);}

const char *infilename="dipfinder2.in";
const long double SMALL=1.0e-10L;
const int SUPPERROR=15;

void makeinitials(void);
void prepoutfiles(void);
void result(void);
void assignshino(long double, int);
void ntmapstep(long double *, long double *, long double);
long double xrootfunc(long double);
long double yrootfunc(long double);
 
int rootbybrak(long double (*)(long double), long double *, long double, long double, int, int);
int rootbybrak2(long double (*)(long double), long double *, long double, long double, int, int);

int zbracup(long double (*)(long double), long double, long double *);
int zbrac(long double (*func)(long double), long double *,  long double *);
void zbrak(long double (*)(long double), long double, long double, int, long double[], long double[], int *);
long double rtbis(long double (*)(long double), long double, long double, long double);
void mkstr(char *, unsigned long int);

long double pi;
long double a0, ax, axx, da, b0, bx, bxx, xacc, yacc, bracs, mfac, comptol;
long double a, b, b2, bifca, bifcb;
long double olda[2], oldb[2], oldm[2], oldb2[2], oldm2[2];

long double shinox, shinoyf, shinopx, shinopyf;
int shp=0, root_no, brakno, sheqbc=1, rootexit=1, rootexit2=1;
unsigned long int qq, pp;
char shinofilename[160], shinofilenamebase[60], shinofilenameend[60], bifcfilename[60], compfilename[60];

FILE *fl1,*fl2,*fl3;

main(int argc, char **argv)
{
    long double tmp1, tmp2, bxint;
    int oldpt=1, errorcount=1;
    
    /* read input */
    makeinitials();

    /* prepare output files */
    prepoutfiles();

    /* calculation */
    shp=0;
    if (bracs<SMALL/2.0)   /* "brute force" "zbrak"eting of a fixed interval*/ 
    {
	b=b0; b2=b0;
	for (a=a0+da;a<=ax;a+=da)
	{
	    if (shp==0) {
		assignshino(a,1);
		rootexit=rootbybrak(&xrootfunc, &b, b, ((b+bx<bxx)? (b+bx): bxx), 
				    brakno, root_no); 
		assignshino(a,2);
		rootexit2=rootbybrak(&xrootfunc, &b2, b2, ((b2+bx<bxx)? (b2+bx): bxx), 
				    brakno, root_no);

		result();
		if (rootexit<=0 || rootexit2<=0) {
		    if (errorcount<=SUPPERROR)
			printf("(%d) Warning, no yroot found at a=%10.5Lf among %d xroot candidates!\n",
			       errorcount++, a, -rootexit);
		} 
	    } 
	    else {
		assignshino(a,shp);
		rootexit=rootbybrak(&xrootfunc, &b, b, ((b+bx<bxx)? (b+bx): bxx), 
				    brakno, root_no); 
		result();
		if (rootexit<=0 && errorcount<=SUPPERROR)
		    printf("(%d) Warning, no yroot found at a=%10.5Lf among %d xroot candidates!\n",
			   errorcount++, a, -rootexit);
		
	    }
	}
    }
    

    else {                 /* "zbrac"eting with interpolation and error bridging */

	a=a0;
	/* find valid starting value */
	if (a<=((long double)(qq))/((long double)(pp))+SMALL) {b=0; rootexit=1; b2=0; rootexit2=1;}
	else {
	    assignshino(a,1);
	    rootexit=rootbybrak(&xrootfunc, &b, b0, bx, brakno, root_no);
	    assignshino(a,2);
	    rootexit=rootbybrak(&xrootfunc, &b2, b0, bx, brakno, root_no);
	}
	result();
	if ((rootexit<=0 || rootexit2<=0) && errorcount<=SUPPERROR)
	    printf("(%d) Warning, no start yroot found at a=%10.5Lf among %d xroot candidates!\n",
		   errorcount++, a, -rootexit);

	a=a0+da;
	while ((rootexit<=0 || rootexit2<=0) && a<=ax) {
	    assignshino(a,1);
	    rootexit=rootbybrak(&xrootfunc, &b, b0, bx, brakno, root_no); 
	    assignshino(a,2);
	    rootexit=rootbybrak(&xrootfunc, &b2, b0, bx, brakno, root_no); 
	    result();
	    if ((rootexit<=0 || rootexit2<=0) && errorcount<=SUPPERROR)
		printf("(%d) Warning, still no start yroot found at a=%10.5Lf among %d xroot candidates!\n",
		       errorcount++, a, -rootexit);
	    a+=da;
	}

	oldpt=1-oldpt; olda[oldpt]=a-da; oldb[oldpt]=b; oldb2[oldpt]=b2;

	/* find second value */
	do {
	    assignshino(a,1);
	    tmp1=b;
	    tmp2=b+bracs;
	    if ((rootexit=zbrac(&xrootfunc, &tmp1, &tmp2))==1)
		tmp1=rtbis(&xrootfunc, tmp1, tmp2, xacc);
	    else if (errorcount<=SUPPERROR)
		printf("(%d) Warning, problems with zbrac at a=%10.5Lf!\n",errorcount++,a);
	    if (rootexit!=1 || fabs(yrootfunc(tmp1))>yacc) {
		if (errorcount<=SUPPERROR)
		    printf("(%d) Warning, found xroot is no yroot at a=%10.5Lf: yrtf(b)=%10.5Lg!\n",
			   errorcount++,a,yrootfunc(tmp1));
		rootexit=rootbybrak(&xrootfunc, &b, b, ((b+bx-b0<bxx)? (b+bx-b0): bxx), 
							brakno, root_no); 
		if (rootexit<=0 && errorcount<=SUPPERROR)
		    printf("(%d) Warning, still no yroot found at a=%10.5Lf among %d xroot candidates!\n",
			   errorcount++, a, -rootexit);
	    } 
	    else b=tmp1;

	    assignshino(a,2);
	    tmp1=b2;
	    tmp2=b2+bracs;
	    if ((rootexit2=zbrac(&xrootfunc, &tmp1, &tmp2))==1)
		tmp1=rtbis(&xrootfunc, tmp1, tmp2, xacc);
	    else if (errorcount<=SUPPERROR)
		printf("(%d) Warning, problems with zbrac at a=%10.5Lf!\n",errorcount++,a);
	    if (rootexit2!=1 || fabs(yrootfunc(tmp1))>yacc) {
		if (errorcount<=SUPPERROR)
		    printf("(%d) Warning, found xroot is no yroot at a=%10.5Lf: yrtf(b)=%10.5Lg!\n",
			   errorcount++,a,yrootfunc(tmp1));
		rootexit2=rootbybrak(&xrootfunc, &b2, b2, ((b2+bx-b0<bxx)? (b2+bx-b0): bxx), 
							brakno, root_no); 
		if (rootexit2<=0 && errorcount<=SUPPERROR)
		    printf("(%d) Warning, still no yroot found at a=%10.5Lf among %d xroot candidates!\n",
			   errorcount++, a, -rootexit2);
	    } 
	    else b2=tmp1;
	   
	    result();
	    a+=da;
	} 
	while ((rootexit<=0 || rootexit2<=0) && a<=ax);

	oldm[oldpt]=(b-oldb[oldpt])/(a-da-olda[oldpt]);
	oldm2[oldpt]=(b2-oldb2[oldpt])/(a-da-olda[oldpt]);
	oldpt=1-oldpt; olda[oldpt]=a-da; oldb[oldpt]=b;	oldb2[oldpt]=b2;	

	/* find third value */
	do {
	    bxint=b+oldm[1-oldpt]*(a-olda[oldpt])*mfac;
	    assignshino(a,1);
	    tmp2=b+bracs;
	    if ((rootexit=zbrac(&xrootfunc, &b, &tmp2))==1)
		tmp1=rtbis(&xrootfunc, b, tmp2, xacc);
	    else if (errorcount<=SUPPERROR)
		printf("(%d) Warning, problems with zbrac at a=%10.5Lf!\n",errorcount++,a);
	    if (rootexit!=1 || fabs(yrootfunc(tmp1))>yacc) {
		if (errorcount<=SUPPERROR)
		    printf("(%d) Warning, found xroot is no yroot at a=%10.5Lf: yrtf(b)=%10.5Lg!\n",
			   errorcount++, a,yrootfunc(tmp1));
		rootexit=rootbybrak(&xrootfunc, &b, b, ((bxx>bxint)? bxint : bxx), 
				    brakno, root_no); 
		if (rootexit<=0 && errorcount<=SUPPERROR)
		    printf("(%d) Warning, still no yroot found at a=%10.5Lf among %d xroot candidates!\n",
			   errorcount++, a, -rootexit);
	    }
	    else b=tmp1;

	    bxint=b2+oldm2[1-oldpt]*(a-olda[oldpt])*mfac;
	    assignshino(a,2);
	    tmp2=b2+bracs;
	    if ((rootexit2=zbrac(&xrootfunc, &b2, &tmp2))==1)
		tmp1=rtbis(&xrootfunc, b2, tmp2, xacc);
	    else if (errorcount<=SUPPERROR)
		printf("(%d) Warning, problems with zbrac at a=%10.5Lf!\n",errorcount++,a);
	    if (rootexit2!=1 || fabs(yrootfunc(tmp1))>yacc) {
		if (errorcount<=SUPPERROR)
		    printf("(%d) Warning, found xroot is no yroot at a=%10.5Lf: yrtf(b)=%10.5Lg!\n",
			   errorcount++, a,yrootfunc(tmp1));
		rootexit2=rootbybrak(&xrootfunc, &b2, b2, ((bxx>bxint)? bxint : bxx), 
				    brakno, root_no); 
		if (rootexit2<=0 && errorcount<=SUPPERROR)
		    printf("(%d) Warning, still no yroot found at a=%10.5Lf among %d xroot candidates!\n",
			   errorcount++, a, -rootexit2);
	    }
	    else b2=tmp1;

	    result();
	    a+=da;
	} 
	while ((rootexit<=0 || rootexit2<=0) && a<=ax);

	oldm[oldpt]=(b-oldb[oldpt])/(a-da-olda[oldpt]);
	oldm2[oldpt]=(b2-oldb2[oldpt])/(a-da-olda[oldpt]);
	oldpt=1-oldpt; olda[oldpt]=a-da; oldb[oldpt]=b;	oldb2[oldpt]=b2;	

	/* main loop */
	a0=a;
	for (a=a0;a<=ax;a+=da)
	{
	    if (shp==0){
		
		tmp2=oldm[1-oldpt]+(oldm[1-oldpt]-oldm[oldpt])*(a-olda[oldpt]); /* interpol. slope */
		bxint=b+tmp2*(a-olda[oldpt])*mfac;
		assignshino(a,1);
		tmp1=b+tmp2*(a-olda[oldpt]);
		tmp2=tmp1+bracs;
		if ((rootexit=zbrac(&xrootfunc, &tmp1, &tmp2))==1)
		    tmp1=rtbis(&xrootfunc, tmp1, tmp2, xacc);
		else if (errorcount<=SUPPERROR)
		    printf("(%d) Warning, problems with zbrac at a=%10.5Lf!\n",errorcount++,a);
		if (rootexit!=1 || fabs(yrootfunc(tmp1))>yacc) {
		    if (errorcount<=SUPPERROR)
			printf("(%d) Warning, found xroot is no yroot at a=%10.5Lf: yrtf(b)=%10.5Lg!\n",
			       errorcount++,a,yrootfunc(tmp1));
		    rootexit=rootbybrak(&xrootfunc, &b, b, ((bxx>bxint)? bxint : bxx), 
					brakno, root_no); 
		    if (rootexit<=0 && errorcount<=SUPPERROR)
			printf("(%d) Warning, still no yroot found at a=%10.5Lf among %d xroot candidates!\n",
			       errorcount++,a, -rootexit);
		}
		else b=tmp1;

		tmp2=oldm2[1-oldpt]+(oldm2[1-oldpt]-oldm2[oldpt])*(a-olda[oldpt]); /* interpol. slope */
		bxint=b2+tmp2*(a-olda[oldpt])*mfac;
		assignshino(a,2);
		tmp1=b2+tmp2*(a-olda[oldpt]);
		tmp2=tmp1+bracs;
		if ((rootexit2=zbrac(&xrootfunc, &tmp1, &tmp2))==1)
		    tmp1=rtbis(&xrootfunc, tmp1, tmp2, xacc);
		else if (errorcount<=SUPPERROR)
		    printf("(%d) Warning, problems with zbrac at a=%10.5Lf!\n",errorcount++,a);
		if (rootexit2!=1 || fabs(yrootfunc(tmp1))>yacc) {
		    if (errorcount<=SUPPERROR)
			printf("(%d) Warning, found xroot is no yroot at a=%10.5Lf: yrtf(b)=%10.5Lg!\n",
			       errorcount++,a,yrootfunc(tmp1));
		    rootexit2=rootbybrak(&xrootfunc, &b2, b2, ((bxx>bxint)? bxint : bxx), 
					brakno, root_no); 
		    if (rootexit2<=0 && errorcount<=SUPPERROR)
			printf("(%d) Warning, still no yroot found at a=%10.5Lf among %d xroot candidates!\n",
			       errorcount++,a, -rootexit);
		}
		else b2=tmp1;

		result();
		
		if (rootexit>0) {
		    oldm[oldpt]=(b-oldb[oldpt])/(a-olda[oldpt]);
		    oldm2[oldpt]=(b2-oldb2[oldpt])/(a-olda[oldpt]);
		    oldpt=1-oldpt; olda[oldpt]=a; oldb[oldpt]=b; oldb2[oldpt]=b2;
		}
		
	    } 
	    
	    else {
		tmp2=oldm[1-oldpt]+(oldm[1-oldpt]-oldm[oldpt])*(a-olda[oldpt]); /* interpol. slope */
		bxint=b+tmp2*(a-olda[oldpt])*mfac;
		assignshino(a,shp);
		tmp1=b+tmp2*(a-olda[oldpt]);
		tmp2=tmp1+bracs;
		if ((rootexit=zbrac(&xrootfunc, &tmp1, &tmp2))==1)
		    tmp1=rtbis(&xrootfunc, tmp1, tmp2, xacc);
		else if (errorcount<=SUPPERROR)
		    printf("(%d) Warning, problems with zbrac at a=%10.5Lf!\n",errorcount++,a);
		if (rootexit!=1 || fabs(yrootfunc(tmp1))>yacc) {
		    if (errorcount<=SUPPERROR)
			printf("(%d) Warning, found xroot is no yroot at a=%10.5Lf: yrtf(b)=%10.5Lg!\n",
			       errorcount++,a,yrootfunc(tmp1));
		    rootexit=rootbybrak(&xrootfunc, &b, b, ((bxx>bxint)? bxint : bxx), 
					brakno, root_no); 
		    if (rootexit<=0 && errorcount<=SUPPERROR)
			printf("(%d) Warning, still no yroot found at a=%10.5Lf among %d xroot candidates!\n",
			       errorcount++,a, -rootexit);
		}
		else b=tmp1;
		result();
		
		if (rootexit>0) {
		    oldm[oldpt]=(b-oldb[oldpt])/(a-olda[oldpt]);
		    oldpt=1-oldpt; olda[oldpt]=a; oldb[oldpt]=b;	
		}

	    }	
	}
    }
    
    if (errorcount>SUPPERROR) 
	printf("Error count exceeds %d, further warnings suppressed...\n", SUPPERROR); 
    fclose(fl1); fclose(fl2); fclose(fl3);
    exit(0); 
} 


void makeinitials(void)      /* read and adjust initial values */
{
  FILE *fl;
  char dummy;
  
  pi=4.0L*atanl(1.0L);

  /* read initials */
  fl=fopen(infilename, "r");
  
  fscanf(fl,"%Lf", &a0); NL(fl, dummy);       /* a0 */
  fscanf(fl,"%Lf", &ax); NL(fl, dummy);       /* amax */
  fscanf(fl,"%Lf", &da); NL(fl, dummy);       /* da */
  fscanf(fl,"%Lf", &b0); NL(fl, dummy);       /* b0 */
  fscanf(fl,"%Lf", &bx); NL(fl, dummy);       /* bmax */
  fscanf(fl,"%Lf", &bxx); NL(fl, dummy);      /* bmaxx */
  fscanf(fl,"%d", &brakno); NL(fl, dummy);    /* # of zbrak intervals */
  fscanf(fl,"%Le", &xacc); NL(fl, dummy);     /* absolute root xacc */
  fscanf(fl,"%Le", &yacc); NL(fl, dummy);     /* absolute root yacc */
  fscanf(fl,"%Le", &bracs); NL(fl, dummy);    /* size of "brac"ket start interval */
  fscanf(fl,"%d", &root_no); NL(fl, dummy);   /* max "brak" intervals stored */
  fscanf(fl,"%Lf", &mfac); NL(fl, dummy);     /* factor for later "brak"ket interval */
  fscanf(fl,"%Le", &comptol); NL(fl, dummy);  /* tolerance for sh-bifc comparison */
  fscanf(fl,"%s %s", shinofilenamebase, shinofilenameend); NL(fl, dummy); /* shinohara output filename */
  fscanf(fl,"%s", bifcfilename); NL(fl, dummy);  /* bifc input filename */
  fscanf(fl,"%s", compfilename); NL(fl, dummy);  /* comparison outfilename */
  
  fclose(fl);
}


void prepoutfiles(void)
{
    char dummy, ctmp[40];

    fl2=fopen(bifcfilename, "r");
    fl3=fopen(compfilename, "a");
    fscanf(fl2, "# bif curve for m/n = %u/%u",&qq,&pp);  NL(fl2,dummy);
    NL(fl2,dummy);
    NL(fl2,dummy);
    NL(fl2,dummy);
    NL(fl2,dummy);
    NL(fl2,dummy);

    /* adjust initials */
    if (a0<=((long double)(qq))/((long double)(pp))) {
	a0=((long double)(qq))/((long double)(pp));
	b0=0.0L;
    }
    if (b0<=0.0L) b0=0.0L;  
    if (ax<a0) ax=a0; 
    axx=ax; /* for stopping comparison at EOF of bifcurve */
    if (bx<b0) bx=b0;
    if (da<=0.0L) da=SMALL;
    if (brakno<=1) brakno=1;
    if (yacc<xacc) yacc=10.0L*xacc;
    if (mfac<1.0L) mfac=1.0L;
    if (bx>bxx) bx=bxx;
    strcpy(shinofilename,shinofilenamebase);
    strcpy(ctmp,""); mkstr(&ctmp[0],qq); strcat(shinofilename,ctmp);
    strcat(shinofilename,"_");
    strcpy(ctmp,""); mkstr(&ctmp[0],pp); strcat(shinofilename,ctmp);
    strcat(shinofilename,shinofilenameend);

    fl1=fopen(shinofilename, "w");
    fprintf(fl1,"# q:       %10u , p:       %10u \n", qq, pp);
    fprintf(fl1,"# a0:      %10.7Lf , b0:      %10.7Lf \n", a0, b0); 
    fprintf(fl1,"# ax:      %10.7Lf \n", ax); 
    fprintf(fl1,"# bx:      %10.7Lf , bxx:     %10.7Lf \n", bx, bxx); 
    fprintf(fl1,"# da:      %10.7Lf , brakno:  %10d \n", da, brakno); 
    fprintf(fl1,"# xacc:    %10.3Le , yacc:    %10.3Le \n", xacc, yacc); 
    fprintf(fl1,"# bracs:   %10.3Le , mfac:    %10.7Lf \n", bracs, mfac); 
    fprintf(fl1,"# root#:   %10d \n", root_no); 
    fprintf(fl1,"# comptol: %10.3Le \n", comptol); 
    fprintf(fl1,"# %18s%20s%5s%5s%10s%20s\n\n", "a", "b", "shp", "s=b","exitcode","yrootfct");
}

void result(void)
{
    char dummy;

    if (shp==0) {
	if (a<=axx && rootexit>0 && rootexit2>0) {
	    do {	    
		fscanf(fl2,"%Le", &bifca); 
		if (&bifca!=NULL) {
		    fscanf(fl2,"%Le", &bifcb); 
		    NL(fl2,dummy); 
		    if (dummy==EOF) axx=a-SMALL;		
		}   
		else axx=a-SMALL;
	    }
	    while (a-bifca>comptol && a<=axx);
	    if (fabs(a-bifca)<=comptol){
		if (fabs(b-bifcb)<=comptol && fabs(b2-bifcb)<=comptol) sheqbc=1;
		else {
		    if (fabs(b-bifcb)<=comptol) {
			sheqbc=1;
			shp=1;
		    }
		    else if (fabs(b2-bifcb)<=comptol) {
			sheqbc=1;
			shp=2;
			b=b2; rootexit=rootexit2;
			oldb[0]=oldb2[0]; oldb[1]=oldb2[1];
			oldm[0]=oldm2[0]; oldm[1]=oldm2[1];
		    } 
		    else {
			if (sheqbc==1)
			    fprintf(fl3,"%20.15Lf %20.15Lf %10u %10u %12.8Le\n", 
				    a, b, qq, pp, b-bifcb);
			sheqbc=0;
		    }
		}
	    }
	    else sheqbc=0;
	}
	else sheqbc=0;
 
	fprintf(fl1,"%20.15Lf%20.15Lf%5d%5d%10d%20.10Lg\n", a, b, shp, sheqbc, rootexit, yrootfunc(b));
    } 
    else {
	if (a<=axx && rootexit>0) {
	    do {	    
		fscanf(fl2,"%Le", &bifca); 
		if (&bifca!=NULL) {
		    fscanf(fl2,"%Le", &bifcb); 
		    NL(fl2,dummy); 
		    if (dummy==EOF) axx=a-SMALL;		
		}   
		else axx=a-SMALL;
	    }
	    while (a-bifca>comptol && a<=axx);
	    if (fabs(a-bifca)<=comptol){
		if (fabs(b-bifcb)<=comptol) sheqbc=1;
		else {
		    if (sheqbc==1)
			fprintf(fl3,"%20.15Lf %20.15Lf %10u %10u %12.8Le\n", 
				a, b, qq, pp, b-bifcb);
		    sheqbc=0;	    
		}
	    }
	    else sheqbc=0;
	}
	else sheqbc=0;

	fprintf(fl1,"%20.15Lf%20.15Lf%5d%5d%10d%20.10Lg\n", a, b, shp, sheqbc, rootexit, yrootfunc(b));
    }
}



/* assign shinohara points or prefactors depending on symmetry line and a */
void assignshino(long double aa, int ss) 
{
	switch (ss) {
	    case 1: {
		shinox=-0.25L; shinoyf=-0.5L; 
		shinopx=-0.25L+((long double)(qq))/2.0L; shinopyf=((2*(qq/2)==qq) ? -0.5L : 0.5L); 
		break;}
	    case 2: {
		shinox=aa/2.0L-0.25L; shinoyf=0.0L; 
		shinopx=aa/2.0L-0.25L+((long double)(qq))/2.0L; shinopyf=0.0L; 
		break;}
	    case 3: {
		shinox=0.25L; shinoyf=0.5L; 
		shinopx=0.25L+((long double)(qq))/2.0L; shinopyf=((2*(qq/2)==qq) ? -0.5L : 0.5L);
		break;}
	    case 4: {
		shinox=aa/2.0L+0.25L; shinoyf=0.0L; 
		shinopx=aa/2.0L+0.25L+((long double)(qq))/2.0L; shinopyf=0.0L; 
		break;}
	    default: {
		shinox=-0.25L; shinoyf=-0.5L;
		shinopx=aa/2.0L+0.25L+((long double)(qq))/2.0L; shinopyf=0.0L;}
	}
}



/* single iteration of NT-map */
void ntmapstep(long double *x, long double *y, long double bb)
{

  *y = *y - bb*sinl(2.0L*pi*(*x));
  *x = *x + a*(1.0L-(*y)*(*y));
  
  /* x modulo 1 
  while (*x> 0.5001L) {*x=*x-1.0L; (*lft)++;}
  while (*x< -0.5L) {*x=*x+1.0L; (*lft)--;} */
}


/*  p/2 map iterations minus shino(m+q) for root search */
long double xrootfunc(long double bb)
{
  int i;
  long double xx, yy;

  xx=shinox; yy=shinoyf*bb;
  for (i=1;i<=pp/2;i++) ntmapstep(&xx, &yy, bb);
  return (xx-shinopx); 
}

/*  p/2 map iterations minus shino(m+q) for root search */
long double yrootfunc(long double bb)
{
  int i;
  long double xx, yy;

  xx=shinox; yy=shinoyf*bb;
  for (i=1;i<=pp/2;i++) ntmapstep(&xx, &yy, bb);
  return (yy-shinopyf*bb); 
}

/* zbrak bracketing and root search and test */
int rootbybrak2(long double (*func)(long double), long double *bb,
	       long double x1, long double x2, int intervals, int maxroots)
{
  long double brackl[maxroots], bracku[maxroots];
  int i, maxbrack=maxroots;

  zbrak(func, x1, x2, intervals, &brackl[0], &bracku[0], &maxbrack);
  i=1; 
  do {
      *bb=rtbis(func, brackl[i], bracku[i], xacc);
      i++;
  }
  while (i<=maxbrack && fabs(yrootfunc(*bb))>yacc);
  if (fabs(yrootfunc(*bb))>yacc) return -i;
  else return i;
}

/* "inline" version of zbrak bracketing 
   with immediate root search and stop at first positive test */
int rootbybrak(long double (*func)(long double), long double *bb,
	       long double x1, long double x2, int intervals, int maxroots)
{
  int j, maxbrack=maxroots;
  
  /* variables for num. rec. zbrak */
  long double x, fp, fc, dx; 
  int nn, nb, nbb, i;
  long double xb1, xb2;
  long double (*fx)(long double);
  

  fx=func; nn=intervals; nb=maxbrack;
  /* num. rec. zbrak */
  nbb=0;
  dx=(x2-x1)/nn;
  fp=(*fx)(x=x1);
  for (i=1;i<=nn;i++) {
    fc=(*fx)(x += dx);
    if (fc*fp <= 0.0) {
      nbb++;
      xb1=x-dx;
      xb2=x;
      
      /* root search and test */
      *bb=rtbis(func, xb1, xb2, xacc);
      i++;
      if (fabs(yrootfunc(*bb))<=yacc) return nbb;
      
      /* back to zbrak */
      if (nb==nbb) break;
    }
    fp=fc;
  }

  /* exit value */
  if (fabs(yrootfunc(*bb))<=yacc) return (nbb);
  else return -nbb;
}



/* slightly modified version of zbrac from Num. Rec. (only searching up, up to bx) */
int zbracup(long double (*func)(long double), long double x1, long double *x2)
{
  const int NTRY=50;
  const long double FACTOR=1.6L;
  const long double TOOBIG=1.0L/SMALL;
  int j;
  long double f1, f2;

  if (x1>=*x2) *x2+=SMALL;
  f1=(*func)(x1);
  f2=(*func)(*x2);
  for (j=1;j<=NTRY;j++) {
    if (f1*f2<0.0) return 1;
    f2=(*func)(*x2 += FACTOR*(*x2-x1));
    if (fabs(*x2)>TOOBIG || *x2>bxx) break;
  }
  return 0;
}

/* zbrac from Num. Rec. */
int zbrac(long double (*func)(long double), long double *x1,  long double *x2)
{
  const int NTRY=50;
  const long double FACTOR=1.6L;
  const long double TOOBIG=1.0L/SMALL;
  int j;
  long double f1, f2;

  if (*x1 == *x2) (*x2 += SMALL);
  f1=(*func)(*x1);
  f2=(*func)(*x2);
  for (j=1;j<=NTRY;j++) {
    if (f1*f2 < 0.0) return 1;
    if (fabs(f1) < fabs(f2))
      f1=(*func)(*x1 += FACTOR*(*x1 - *x2));
    else
      f2=(*func)(*x2 += FACTOR*(*x2 - *x1));
    if (fabs(*x2)>TOOBIG || *x1<b0 || *x2>bxx) break;
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
long double rtbis(long double (*func)(long double), long double x1, long double x2, long double xxacc)
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
    if (fabs(dx) < xxacc || fmid == 0.0) return rtb;
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

