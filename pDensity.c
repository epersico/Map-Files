//TO DO LIST

//a,b loops		done?
//m loops		 	done?
//Variable Bracket size
//ask about shinohara
//set up output file architecture





#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>

#include <mpi.h>
#include <omp.h>

#define ROOT_NO 100000
#define SMALL 1.0e-10L
#define itstime2mod 1000
#define SMOOTH(fs) (fs >0) ? (fs)/(1.0+fabsl(fs)): (fs)/(1.0-fabsl(fs))
#define NL(flnm,dumch) {do dumch=fgetc(flnm); while (dumch!='\n' && dumch!=EOF);}

const char *infilename="Density.in";
char outfilename[4][180];

void makeinitials(void);
void prepoutfiles(void);
void ntmapstep(long double *, long double *, long double, long double);
void invntmapstep(long double *, long double *, long double, long double);
void ntdm(long double, long double, long double[2][2]);
long double findx0(long double, long double, int);
void assignshino(int, long double *, long double *);
long double fiter2funcoe0(long double);
long double fiter2funcoo0(long double);
long double fiter2funceo0(long double);
long double fiter2funcoe1(long double);
long double fiter2funcoo1(long double);
long double fiter2funceo1(long double);
int findshinofp(int, long double *, long double *, long double *, long double *);
int wind(long double, long double);
long double residue(int, int *, int *, long double *,long double *, long double *);
 
void zbrak(long double (*)(long double), long double, long double, int, long double[], long double[], int *);
long double rtbis(long double (*)(long double), long double, long double, long double);
long int gcd(long int, long int);
long double removeZeroes(long double *arr,int);
double mysecond();
void initialize(long double *arr, int);

long double pi;
long double a0, aStepSize, b0, bStepSize, acc, lacc, eqd, brakl0, braku0;
long double xorbs[ROOT_NO], yorbs[ROOT_NO];
int  sln, outfiles,n_aSteps, n_bSteps;
long int n, m, n_BaseBraks, n_Braks, nMax;
void (*map)(long double *, long double *, long double, long double);
FILE *fl1, *fl2;

int main(int argc, char **argv)
{
      

//MPI STUFF

int ierr, irank,nranks;

ierr = MPI_Init(&argc,&argv);
ierr = MPI_Comm_size(MPI_COMM_WORLD,&nranks);
ierr = MPI_Comm_rank(MPI_COMM_WORLD,&irank);

//END OF MPI SET UP

	long double phalf, shlamx, shlamy;
	long double brackl[ROOT_NO], bracku[ROOT_NO];
	initialize(brackl,ROOT_NO);
	initialize(bracku,ROOT_NO);
	int lift=0, maxbrack=(ROOT_NO-8)/4;
	int i, j, meql, meqlold=0, idorb, idorbiter;
	long double (*rootfunc)(long double);
	double t0,t1,t2;
	int i_a,i_b,i_m;  //Indexes to cycle through the a and b values
	
	int nt = omp_get_max_threads();
	int id;
	char astr[50],bstr[50];
	long double *pery;
	pery = (long double *)malloc( sizeof(long double)*ROOT_NO );
	/* read input */
	makeinitials();
	/* prepare output files */

 
  long double *qq, *pp;
  qq = (long double *)malloc( sizeof(long double)*(nMax+1) );
  pp = (long double *)malloc( sizeof(long double)*(nMax+1) );
  if(qq==NULL || pp==NULL){
    printf("Unable to allocate memory");
    exit(1);
  } 

aStepSize = (1/(long double)n_aSteps);
bStepSize = (1/(long double)n_bSteps);

printf("There are going to be %d OMP threads\n",omp_get_max_threads()); 
t2=mysecond();
for (i_a=irank;i_a<=n_aSteps;i_a+=nranks) { //Start of "a" loop
	a0 = i_a * aStepSize;

	for (i_b=0; i_b<=n_bSteps; i_b++){// Start of "b" loop
		b0 = i_b * bStepSize;
		t1 = mysecond();
		 // a0 =  0.94;
		 // b0 = 0.2059;
		printf("MPI Process %d has (a,b)=(%Lf,%Lf)\n",irank,a0,b0);
		printf("------ \n\n\n\n NEW (a,b) VALUES: (%Lf,%Lf) \n\n\n\n -------\n",a0,b0);
		
		sprintf(astr,"a_%Lf_",a0);
		sprintf(bstr,"b_%Lf_",b0);
	
		strcpy(outfilename[3],"./data/");
		strcat(outfilename[3],astr);
		strcat(outfilename[3],bstr);
		strcat(outfilename[3],outfilename[0]);
		fl1=fopen(outfilename[3], "w");
		if (outfiles>1) {
		  strcpy(outfilename[3],"./data/");
		  strcat(outfilename[3],astr);
		  strcat(outfilename[3],bstr);
		  strcat(outfilename[3],outfilename[1]);
		  fl2=fopen(outfilename[3], "w");
		}

		//   /* limits of shino point orbits */
		//   map=&ntmapstep;
		//   for (sln=1;sln<=4;sln++) {
		//       ++meql;
		//       maxbrack=findshinofp(sln, &xorbs[meql], &yorbs[meql], &shlamx, &shlamy);
		//       lift=wind(xorbs[meql], yorbs[meql]); /* optional output to fl3 */
		//       fprintf(fl1,"# shp=%-2d, lambdax=%-15.7Lg, lambday=%-15.7Lg\n", sln, shlamx, shlamy);  
		//       if (lift!=m || maxbrack<=0) meql--;
		//       else {
		//    fprintf(fl1,"%15.11Lf %15.11Lf %10.3Lg ", 
		//        xorbs[meql], yorbs[meql], residue(meql, &idorb, &idorbiter, &phalf));
		//    fprintf(fl1,"%10.5Lf %6d %6d %6d\n\n", phalf, meql, idorb, idorbiter);
		//       }      
		//   }
		//   printf("75\n");
		//   /*  limits of shino point inverse orbits */
		//   map=&invntmapstep;
		//   for (sln=1;sln<=4;sln++) {
		//       ++meql;
		//       maxbrack=findshinofp(sln, &xorbs[meql], &yorbs[meql], &shlamx, &shlamy);
		//       lift=wind(xorbs[meql], yorbs[meql]); /* optional output to fl3 */
		//       fprintf(fl1,"# shp=-%-2d, lambdax=%-15.7Lg, lambday=%-15.7Lg\n", sln, shlamx, shlamy);  
		//       if (lift!=-m+1 || maxbrack<=0) meql--;
		//       else {
		//    fprintf(fl1,"%15.11Lf %15.11Lf %10.3Lg ", 
		//        xorbs[meql], yorbs[meql], residue(meql, &idorb, &idorbiter, &phalf));
		//    fprintf(fl1,"%10.5Lf %6d %6d %6d\n\n", phalf, meql, idorb, idorbiter);
		//       }      
		//   }       

		//Here is where the loop over periodic orbits needs to happen. Fix denominator, iterate over numerator 
		for (i_m=1; i_m<nMax; i_m++){
   			meql=0;			
			//if(i_m==3) exit(0);

			n = nMax/gcd(i_m,nMax);
			m = i_m/gcd(i_m,nMax);
			n_Braks = n_BaseBraks * n;
			printf("Now doing %ld / %ld \n",m,n);
			prepoutfiles();
			fprintf(fl1,"# a=%15.11Lf, b=%15.11Lf\n", a0, b0);
			if (outfiles>1) fprintf(fl2,"# a=%15.11Lf, b=%15.11Lf\n", a0, b0);
			if (outfiles>1) fprintf(fl2,"# %12s %14s\n", "q0", "p0");
		//Write a function to make output files for each winding number.  Will make a file tree probably.
		//File for given a,b value.  Then make a file for each orbit in there.  

		/* fixed points on symmetry lines 1-4.  I set the top of the for loop to 1 so it only
		 searches the root symmetry line.  */
			map=&ntmapstep;
			for (sln=1;sln<=1;sln++) {
				meqlold=meql;
				/* find root number and brackets */
				maxbrack=(ROOT_NO-8)/4;
				if (2*(n/2)==n) {        /* p even (-> q odd) */
					if (sln<=2) rootfunc=&fiter2funcoe0;
					else rootfunc=&fiter2funcoe1; 
				}
				else {                               /* p odd */
					if (sln<=2) {
						if (2*(m/2)==m)      /* q even */
							rootfunc=&fiter2funceo0;        /* I_0^i -> I_1^i after (p+1)/2, q/2 */
						else rootfunc=&fiter2funcoo0;       /* I_0^i -> I_1^(1-i) after (p+1)/2, (q+-1)/2 */
					}
					else {
						if (2*(m/2)==m)      /* q even */
							rootfunc=&fiter2funceo1;        /* I_0^1 -> I_1^1 after (p+1)/2, q/2 */
						else rootfunc=&fiter2funcoo1;       /* I_0^1 -> I_1^0 after (p+1)/2, (q+1)/2 */
					}
				}
					
					long double intervalLength=braku0-brakl0, intervalStep=intervalLength/nt,braklower,brakupper;
					t0=mysecond();			
					maxbrack=0;
					initialize(brackl,ROOT_NO);
					initialize(bracku,ROOT_NO);
					#pragma omp parallel private(braklower,brakupper,id) reduction(+:maxbrack)
					{
					id = omp_get_thread_num();		
					braklower = brakl0+id*intervalStep;
					brakupper = brakl0+(id+1)*intervalStep;
		//			printf("Thread %d is finding roots in range (%Lf,%Lf)\n",id,braklower,brakupper);
					zbrak(rootfunc, braklower,brakupper, n_Braks, &brackl[n_Braks*id], &bracku[n_Braks*id], &maxbrack);
					}
					
					fprintf(fl1,"# sln=%-2d, maxbrack=%-4d, ", sln, maxbrack);  /* sln, number of roots found */
					printf("Finding %d brackets took %f (sec)\n",maxbrack,mysecond()-t0);

					t0=mysecond();
					//for(i=0;i<=ROOT_NO-1;i++) printf("%d=(%Lf,%Lf)\n",i,brackl[i],bracku[i]);
					removeZeroes(brackl,maxbrack);
					removeZeroes(bracku,maxbrack);
					printf("Removing zeroes took %f(sec)\n",mysecond()-t0);

					/* examine root brackets */
					//for(i=0;i<=ROOT_NO-1;i++) printf("%d=(%Lf,%Lf)\n",i,brackl[i],bracku[i]);
					t0=mysecond();
					#pragma omp parallel for
					for (i=0;i<=maxbrack-1;i++) {
					  pery[i]=rtbis(rootfunc, brackl[i], bracku[i], acc);
					  }
					printf("Finding roots took %f(sec)\n",mysecond()-t0);
					t0=mysecond();
					for (i=0;i<maxbrack;i++) {
			 	          xorbs[meql]=findx0(a0, pery[i], sln);yorbs[meql]=pery[i];
					  meql++;  
					}
			
					fprintf(fl1,"meql=%-4d\n", meql-meqlold);
					
					
					for (i=meqlold;i<meql;i++) {
				   	  fprintf(fl1,"%15.11Lf %15.11Lf %10.3Lg ", 
					    xorbs[i], yorbs[i], residue(i, &idorb, &idorbiter, &phalf, qq,pp));
					  fprintf(fl1,"%10.5Lf %6d %6d %6d\n", phalf, i, idorb, idorbiter);
					}
					
					fprintf(fl1,"\n");   
					printf("Writing to file and resiude calc took %f sec\n",mysecond()-t0);
			}

			

		  //fprintf(fl1,"Calculation took %f sec\n\n\n",t1);
				


		} //end of "m" loop

	t1 = mysecond()-t1;
	printf("root finding took %f sec for (a,b)=(%Lf,%Lf)\n",t1,a0,b0);

	fclose(fl1);
	if (outfiles>1) fclose(fl2);	
	} //end of "b" loop

} //end of "a" loop
printf("Whole program took %f sec\n",mysecond()-t2);
MPI_Finalize();
return 0;
} 


void makeinitials(void)      /* read and adjust initial values */
{
	FILE *fl;
	char dummy;
	int i;

	pi=4.0L*atanl(1.0L);

	/* read initials */
	fl=fopen(infilename, "r");

	fscanf(fl,"%d",&n_aSteps);  NL(fl,dummy);  //a steps
	fscanf(fl,"%d",&n_bSteps);  NL(fl,dummy);  //b steps
	fscanf(fl,"%ld",&nMax); NL(fl,dummy);  //highest periodic orbit
	fscanf(fl,"%Lf",&acc);  NL(fl,dummy);  //absolute root acc
	//Number of output files

	fscanf(fl,"%Le", &eqd); NL(fl, dummy);      /* absolute root acc */
	fscanf(fl,"%Lf", &brakl0); NL(fl, dummy);   /* lower end of interval for "brak"ket */
	fscanf(fl,"%Lf", &braku0); NL(fl, dummy);   /* upper end of interval for "brak"ket */
	fscanf(fl,"%ld",&n_BaseBraks);  NL(fl,dummy);  //divisions for "brak"ket interval baseline
	fscanf(fl,"%d", &outfiles); NL(fl, dummy); /* number of output files */

	for (i=0;i<outfiles;i++) {
	    fscanf(fl,"%s", &outfilename[i][0]); NL(fl, dummy);  /* bifc input filename */
	}

	//put the output filename loop here.  Or don't.  Make then Procedurally. 

	
	// fscanf(fl,"%Lf", &a0); NL(fl, dummy);       /* start a */
	// fscanf(fl,"%Lf", &b0); NL(fl, dummy);       /* start b */
	// fscanf(fl,"%ld", &m); NL(fl, dummy);        /* m */
	// fscanf(fl,"%ld", &n); NL(fl, dummy);        /* period n */
	// fscanf(fl,"%Le", &acc); NL(fl, dummy);      /* absolute root acc */
	// fscanf(fl,"%Le", &lacc); NL(fl, dummy);     /* shino exp acc */
	// fscanf(fl,"%Le", &eqd); NL(fl, dummy);      /* absolute root acc */
	// fscanf(fl,"%Lf", &brakl0); NL(fl, dummy);   /* lower end of interval for "brak"ket */
	// fscanf(fl,"%Lf", &braku0); NL(fl, dummy);   /* upper end of interval for "brak"ket */
	// fscanf(fl,"%ld", &n_Braks); NL(fl, dummy);    /* divisions for "brak"ket interval */
	// fscanf(fl,"%d", &outfiles); NL(fl, dummy); /* number of output files */

	// for (i=0;i<outfiles;i++) {
	//     printf("168\n");
	//     fscanf(fl,"%s", &outfilename[i][0]); NL(fl, dummy);  /* bifc input filename */
	// }
	fclose(fl);
	/* adjust initials */
	if (brakl0==braku0) {brakl0=-0.5L; braku0=0.5L;}
}


void prepoutfiles(void)
{
	fprintf(fl1,"# a:       %10.7Lf , b:      %10.7Lf \n", a0, b0); 
	fprintf(fl1,"# omega:   (m =%4ld) / (n =%4ld) \n", m, n); 
	fprintf(fl1,"# brakl0:  %10.7Lf , braku0:  %10.7Lf , n_Braks:    %10ld\n", 
		brakl0, braku0, n_Braks); 
	fprintf(fl1,"# acc:     %10.3Le , lacc:      %10.3Le , eqd:     %10.3Le\n", acc, lacc, eqd); 
	fprintf(fl1,"# %13s %15s %10s %10s %6s %6s %6s\n", 
		"x0", "y0", "res", "yhalf", "orb#", "eq to", "@ i="); 

	if (outfiles>1) {
		fprintf(fl2,"# a:       %10.7Lf , b:      %10.7Lf \n", a0, b0); 
		fprintf(fl2,"# omega:   (m =%4ld) / (n =%4ld) \n", m, n); 
		fprintf(fl2,"# brakl0:  %10.7Lf , braku0:  %10.7Lf , n_Braks:    %10ld\n", 
			brakl0, braku0, n_Braks); 
		fprintf(fl2,"# acc:     %10.3Le , lacc:      %10.3Le , eqd:     %10.3Le\n", acc, lacc, eqd); 
	}
	
}

void ntmapstep(long double *xxx, long double *yyy, long double aaa, long double bbb) 
{
	(*yyy)-=bbb*sinl(2.0L*pi*(*xxx)); 
	(*xxx)+=aaa-aaa*(*yyy)*(*yyy);
}

// void invntmapstep(long double *xxx, long double *yyy, long double aaa, long double bbb) 
// {
//   (*xxx)-=aaa-aaa*(*yyy)*(*yyy);
//   (*yyy)+=bbb*sinl(2.0L*pi*(*xxx)); 
// }


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


// /* assign shinohara points */
// void assignshino(int ss, long double *xval, long double *yval)
// {
//   switch (sln) {
//       case 1: {*xval=-0.25L; *yval=-0.5L*b0; break;}
//       case 2: {*xval=a0/2.0L-0.25L; *yval=0.0L; break;}
//       case 3: {*xval=0.25L; *yval=0.5L*b0;break;}
//       case 4: {*xval=a0/2.0L-0.25L; *yval=0.0L; break;}
//       default: {*xval=-0.25L; *yval=-0.5L*b0; printf("Uhoh - shinopoint...\n");}
//   }

//   *xval=remainderl(*xval,1.0L);
// }

//new fiter functions

		
long double fiter2funcoe0(long double yy0) 
{
	unsigned long int i, j, jj, qf, pf;
	long double xxtmp, yytmp;

	pf=n; if (sln==1) qf=m-1; else qf=m+1;
	xxtmp=findx0(a0, yy0, sln); yytmp=yy0; //jj=itstime2mod; j=jj*pf/qf;
	for (i=1;i<=n/2;i++) {
		map(&xxtmp, &yytmp, a0, b0); 
		// if (i>=j){ xxtmp-=itstime2mod; /*jj+=itstime2mod; j=jj*pf/qf; */} /* mod enat to have subtracted (q+-1)/2 in the end */ 
	}
	xxtmp-=((long double)(qf/2/*-(jj-itstime2mod)*/));
	if (sln==1) return SMOOTH(xxtmp-0.5L);
	else return SMOOTH(xxtmp);
}

long double fiter2funcoe1(long double yy0) 
{
	unsigned long int i, j, jj, qf, pf;
	long double xxtmp, yytmp;

	pf=n; if (sln==3) qf=m-1; else qf=m+1;
	xxtmp=findx0(a0, yy0, sln); yytmp=yy0; /*jj=itstime2mod; j=jj*pf/qf*/;
	for (i=1;i<=n/2;i++) {
		map(&xxtmp, &yytmp, a0, b0); 
		// if (i>=j){ xxtmp-=itstime2mod; jj+=itstime2mod; j=jj*pf/qf; } /* mod enat to have subtracted (q+-1)/2 in the end */ 
	}
	xxtmp-=((long double)(qf/2/*-(jj-itstime2mod)*/));
	if (sln==3) return SMOOTH(xxtmp-a0/2.0L*(1.0L-yytmp*yytmp)-0.5L);
	else return SMOOTH(xxtmp-a0/2.0L*(1.0L-yytmp*yytmp));
}

long double fiter2funcoo0(long double yy0) 
{
	unsigned long int i, j, jj, qf, pf;
	long double xxtmp, yytmp;

	pf=n+1; if (sln==1) qf=m-1; else qf=m+1;	
	xxtmp=findx0(a0, yy0, sln); yytmp=yy0; /*jj=itstime2mod; j=jj*pf/qf;*/	
	for (i=1;i<=(n+1)/2;i++) {
		map(&xxtmp, &yytmp, a0, b0);
		// if (i>=j){ xxtmp-=itstime2mod; jj+=itstime2mod; j=jj*pf/qf; } /* mod enat to have subtracted (q+-1)/2 in the end */ 
	}
	xxtmp-=((long double)(qf/2/*-(jj-itstime2mod)*/));
	if (sln==1) return SMOOTH(xxtmp-a0/2.0L*(1.0L-yytmp*yytmp)-0.5L);
	else return SMOOTH(xxtmp-a0/2.0L*(1.0L-yytmp*yytmp));
}

long double fiter2funcoo1(long double yy0) 
{
	unsigned long int i, j, jj, qf, pf;
	long double xxtmp, yytmp;

	pf=n-1; if (sln==3) qf=m-1; else qf=m+1;
	xxtmp=findx0(a0, yy0, sln); yytmp=yy0; /*jj=itstime2mod; j=jj*pf/qf;*/
	for (i=1;i<=(n-1)/2;i++) {
		map(&xxtmp, &yytmp, a0, b0); 
		// if (i>=j){ xxtmp-=itstime2mod; jj+=itstime2mod; j=jj*pf/qf; } /* mod enat to have subtracted (q+-1)/2 in the end */ 
	}
	xxtmp-=((long double)(qf/2/*-(jj-itstime2mod)*/));
	if (sln==3) return SMOOTH(xxtmp-0.5L);
	else return SMOOTH(xxtmp);
}

long double fiter2funceo0(long double yy0) 
{
	unsigned long int i, j, jj, qf, pf;
	long double xxtmp, yytmp;

	pf=n+1; qf=m;
	xxtmp=findx0(a0, yy0, sln); yytmp=yy0;/* jj=itstime2mod; j=jj*pf/qf;*/
	for (i=1;i<=(n+1)/2;i++) {
		map(&xxtmp, &yytmp, a0, b0); 
		// if (i>=j){ xxtmp-=itstime2mod; jj+=itstime2mod; j=jj*pf/qf; } /* mod enat to have subtracted q/2 in the end */ 
	}
	xxtmp-=((long double)(qf/2/*-(jj-itstime2mod)*/));
	if (sln==1) return SMOOTH(xxtmp-a0/2.0L*(1.0L-yytmp*yytmp));
	else return SMOOTH(xxtmp-a0/2.0L*(1.0L-yytmp*yytmp)-0.5L);
}

long double fiter2funceo1(long double yy0) 
{
	unsigned long int i, j, jj, qf, pf;
	long double xxtmp, yytmp;

	pf=n-1; qf=m;
	xxtmp=findx0(a0, yy0, sln); yytmp=yy0; /*jj=itstime2mod; j=jj*pf/qf;*/
	for (i=1;i<=(n-1)/2;i++) {
		map(&xxtmp, &yytmp, a0, b0); 
		// if (i>=j){ xxtmp-=itstime2mod; jj+=itstime2mod; j=jj*pf/qf; } /* mod enat to have subtracted q/2 in the end */ 
	}
	xxtmp-=((long double)(qf/2/*-(jj-itstime2mod)*/));
	if (sln==3) return SMOOTH(xxtmp);
	else return SMOOTH(xxtmp-0.5L);
}

double mysecond()
{
   struct timeval tp;
   struct timezone tzp;
   int i;
   i = gettimeofday(&tp,&tzp);
   return ( (double) tp.tv_sec + (double) tp.tv_usec * 1.e-6 );
}

// /* find limiting fixed point of shinohara orbit */
// int findshinofp(int shp, long double *shx, long double *shy, long double *shlx, long double *shly)
// {
//     long double xx0, yy0, xx, yy, dxx, dyy;
//     int i, j;

//     assignshino(shp, &xx0, &yy0); xx0=remainderl(xx0,1.0L);
//     xx=xx0; yy=yy0;

//     for (i=1;i<=9*n;i++) map(&xx,&yy,a0,b0); xx=remainderl(xx,1.0L);
//     xx0=xx; yy0=yy;
//     for (i=1;i<=n;i++) map(&xx,&yy,a0,b0); xx=remainderl(xx,1.0L);
//     dxx=remainderl(xx-xx0,1.0L); dyy=yy-yy0;
//     *shlx=exp(log(fabsl(dxx))/((long double) (10)));
//     *shly=exp(log(fabsl(dyy))/((long double) (10)));
//     xx0=xx; yy0=yy;
//     for (i=1;i<=n;i++) map(&xx,&yy,a0, b0); xx=remainderl(xx,1.0L);
//     j=11;
		
//     while (fabsl(remainderl(xx-xx0,1.0L))<=fabsl(dxx)+acc 
//     && fabsl(yy-yy0)<=fabsl(dyy)+acc && fabsl(yy-yy0)>acc && 
//     /*(fabsl(*shlx-exp(log(fabsl(remainderl(xx-xx0,1.0L)))/((long double) (j))))>lacc ||
//       fabsl(*shly-exp(log(fabsl(yy-yy0))/((long double) (j))))>lacc) &&*/ j<ROOT_NO) {
//  dxx=(fabsl(remainderl(xx-xx0,1.0L))>acc? remainderl(xx-xx0,1.0L):acc); dyy=yy-yy0;
//  if (fabsl(remainderl(dxx,1.0L))<fabsl(dyy)) dxx=dyy;
//  *shlx=exp(log(fabsl(dxx))/((long double) (j)));
//  *shly=exp(log(fabsl(dyy))/((long double) (j)));
//  xx0=xx; yy0=yy;
//  for (i=1;i<=n;i++) map(&xx,&yy,a0,b0); xx=remainderl(xx,1.0L);
//  j++;
//     }
//     *shx=xx; *shy=yy;
//     if (fabsl(yy-yy0)<acc) {
//  *shx=xx; 
//  *shy=yy;  
//  *shlx=exp(log(fabsl(dxx))/((long double) (j)));
//  *shly=exp(log(fabsl(dyy))/((long double) (j)));
//  return 1;
//     }
//     else if (fabsl(remainderl(xx-xx0,1.0L))<=fabsl(dxx)+5.0L*acc 
//       && fabsl(yy-yy0)<=fabsl(dyy)+5.0L*acc && j<ROOT_NO) {
//  *shlx=exp(log(fabsl(remainderl(xx-xx0,1.0L)))/((long double) (j)));
//  *shly=exp(log(fabsl(yy-yy0))/((long double) (j)));
//  if (remainderl(xx-xx0,1.0L)>0) xx=xx0+pow(*shlx,((long double) (j)))/(1.0L-(*shlx));
//  else xx=xx0-pow(*shlx,((long double) (j)))/(1.0L-(*shlx));
//  if (yy>yy0) yy=yy0+pow(*shly,((long double) (j)))/(1.0L-(*shly));
//  else yy=yy0-pow(*shly,((long double) (j)))/(1.0L-(*shly));
//  xx=remainderl(xx,1.0L); xx0=xx; yy0=yy;
//  for (i=1;i<=n;i++) map(&xx,&yy,a0,b0); xx=remainderl(xx,1.0L);
//  if (fabsl(yy-yy0)<acc) {*shx=xx; *shy=yy; return 1;} else {
//      while (fabsl(remainderl(xx-xx0,1.0L))<=fabsl(dxx)+acc 
//       && fabsl(yy-yy0)<=fabsl(dyy)+acc && fabsl(yy-yy0)>acc && j<ROOT_NO) {
//    dxx=(fabsl(remainderl(xx-xx0,1.0L))>acc? remainderl(xx-xx0,1.0L):acc); dyy=yy-yy0;
//    if (fabsl(remainderl(dxx,1.0L))<fabsl(dyy)) dxx=dyy;
//    *shlx=exp(log(fabsl(dxx))/((long double) (j)));
//    *shly=exp(log(fabsl(dyy))/((long double) (j)));
//    xx0=xx; yy0=yy;
//    for (i=1;i<=n;i++) map(&xx,&yy,a0,b0); xx=remainderl(xx,1.0L);
//    j++;
//      }
//      if (fabsl(yy-yy0)<acc) {
//    *shx=xx; 
//    *shy=yy;  
//    *shlx=exp(log(fabsl(remainderl(xx-xx0,1.0L)))/((long double) (j)));
//    *shly=exp(log(fabsl(yy-yy0))/((long double) (j)));
//    return 1;
//      } 
//      else return 0;
//  } 
//     }
//     else return 0;    
// }

// /* find winding number of known periodic orbit, with optional output to fl3 */
// int wind(long double q, long double p)
// {
//   long double qq=q, pp=p;
//   int i;

//   for (i=1;i<=n;i++) map(&qq, &pp, a0, b0);
//   return ((int)(qq-q+SMALL));
// }

/* find residues of a known periodic orbit, with opt. output to fl2 and *tmp on I0 */
long double residue(int ii, int *iieq, int *iiter, long double *tmp, long double *qq, long double *pp)
{
  long double dpq[2][2], lm[2][2], lm0[2][2];
  int i, j, k, lft=0;
  *iieq=ii;
  *iiter=0;
  qq[0]=xorbs[ii]; qq[0]=remainderl(qq[0],1.0L);
  pp[0]=yorbs[ii];
  ntdm(qq[0],pp[0],lm);

  for (k=1;k<n;k++) {

    if (outfiles>1) fprintf(fl2,"%14.8Lf %14.8Lf\n", qq[k-1], pp[k-1]); 
    qq[k]=qq[k-1];
    pp[k]=pp[k-1];
    map(&qq[k], &pp[k], a0, b0); qq[k]=remainderl(qq[k],1.0L);
    if (k==n/2) *tmp=pp[k];
    ntdm(qq[k],pp[k],dpq);
    for (i=0;i<=1;i++) for (j=0;j<=1;j++) lm0[i][j]=lm[i][j]; 
    for (i=0;i<=1;i++) for (j=0;j<=1;j++) lm[i][j]=dpq[i][0]*lm0[0][j]+dpq[i][1]*lm0[1][j];
    i=1; while (i<ii) {
 	if (fabsl(pp[k]-yorbs[i])<eqd && fabsl(qq[k]-xorbs[i])<eqd && i<(*iieq)) {*iieq=i; *iiter=k;}
	 i++;
    }    
  }

  if (outfiles>1) fprintf(fl2,"%14.8Lf %14.8Lf ", qq[k], pp[k]);
  map(&qq[k], &pp[k], a0, b0); qq[k]=remainderl(qq[k],1.0L);
  fprintf(fl2,"%14.4Le %14.4Le\n\n", qq[k]-findx0(yorbs[i], a0, sln), pp[k]-yorbs[i]);
  return (2.0L-lm[0][0]-lm[1][1])/4.0L;
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
		if (fabsl(dx) < xacc || fmid == 0.0) return rtb;
	}
	nrerror=2; /* "Too many bisections in rtbis" */
	return 0.0L;
}

long int gcd(long int n, long int m)
{
  int remainder;

  while(n != 0)
  {
    remainder = m % n;
    m = n;
    n = remainder;
  }

  return m;
}

long double removeZeroes(long double *arr,int n)
{
int current=0,next=0;
for (current=0; current<=n-1; current++) 
    {
	while(arr[next] == 0) {
            next++;
        }
        arr[current] = arr[next];
//	printf("arr[%d]=%Lf next = %d n = %d\n",current,arr[current],next,n);
	next++;
    }
}

void initialize(long double *arr, int len)
{
int i;
for (i=0;i<len-1;i++) arr[i]=0;
}
