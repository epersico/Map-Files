#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>

#define ROOT_NO 10000
#define SMALL 1.0e-10L
#define itstime2mod 1000
#define SMOOTH(fs) (fs >0) ? (fs)/(1.0+fabsl(fs)): (fs)/(1.0-fabsl(fs))
#define NL(flnm,dumch) {do dumch=fgetc(flnm); while (dumch!='\n' && dumch!=EOF);}

const char *infilename="temp.txt";

void makeinitials(void);
void prepoutputfiles(void);
void ntmapstep(long double *, long double *, long double, long double);
long double windingNumber(long double, long double);

//OUTLINE

//Read in file that will specify 
//	1. Filename of orbits to find.
//	2. Pattern to search for base root?
//	3. Step Size to probe up

//For given orbit
//	1. Find distance to Shinohara point
//	2. Increase y0 and find finding number
//	3. Plot

void makeinitials();
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