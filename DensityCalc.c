//
//  main.c
//  ResonanceWidths
//
//  Created by Elliot Persico on 5/10/17.
//  Copyright © 2017 Elliot Persico. All rights reserved.
//

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

//
#define NL(flnm,dumch) {do dumch=fgetc(flnm); while (dumch!='\n' && dumch!=EOF);}
#define UNDEFINED -123456789.0
#define nParameters 1
#define nMax 20
#define pi  3.1415926535897932384626433832795028841971693993751058209749445923078164062862089986280348253421170679L


void ntmapstep(long double *, long double *, long double, long double);
void standardmapstep(long double *, long double *, long double,long double);
void (*map)(long double *, long double *, long double, long double);
void find_om(long double, long double, long double, long double, long double *, long double *,
             int *, long double *, int *,
             long double *, int *, long double *, int *);
void assignshino(long double, long double, int, long double *, long double *);
void findShinoOrbitYonSLN(long double, long double, long double *,long double *);
void ResonanceWidths(long double, long double, long double, long double,int, int, int, long double *);



long double cutoff,dn,speedup;
//long double pi;  //This was from before I did a #define of pi.

char infilename[64],outfilename[64], densityFilename[60],distanceFilename[60];
FILE *fl, *fld, *fldist; //fl is the file for the individual winding number, flw is the one just for the widths.

int main(int argc, const char * argv[]) {
    
    //Start of the program.  Load in the data from Density.c  This has a particular format.
    sprintf(infilename,"./data/orbits_3_n_20.out");
    char dummy;
    int n,m,nOrbits[nParameters];
    long double a0[nParameters],b0[nParameters],x0[nParameters][2],y0[nParameters][2];
    int sln, nRoots,i,j,updown[nParameters][2];
    long double distance[nParameters];
    long double shinox[nParameters],shinoy[nParameters],dum1,dum2;
    
    // Uncomment the map I want
    map = &ntmapstep;
    //map = &standardmapstep;
    
    printf("The filename opened was %s\n",infilename);
    fl = fopen(infilename,"r");

    //Loop through the input file.  The organization is that for each (a,b) parameter pair there will be a number of roots to iterate through.
        for(i=0;i<nParameters;i++){
          nOrbits[i] = 0;
          for(j=0;j<=nMax;j++){  
            
            fscanf(fl,"# a: %La (%Le) , b: %La (%Le)",&a0[i],&dum1,&b0[i],&dum2);   NL(fl,dummy);    
            if(j==0) printf("(a,b) = (%Lf,%Lf)\n",a0[i],b0[i]);
        
            fscanf(fl, "# omega:   (m =   %d) / (n =  %d)",&m,&n); NL(fl,dummy);
            printf("m/n = %d/%d\n",m,n);
            fscanf(fl,"# shinox: %La (%Le) , shinoy: %La (%Le)",&shinox[i],&dum1,&shinoy[i],&dum2); NL(fl,dummy);
            if(j==0)printf("(shinox,shinoy) = (%Lf,%Lf)\n",shinox[i],shinoy[i]);
            //Initialization of the widths and density file.
            if(i==0){
            sprintf(densityFilename,"./density/a_%Lf_density.out",a0[0]);
            fld = fopen(densityFilename,"w");
            sprintf(distanceFilename,"./density/a_%Lf_distance.out",a0[0]);
            fldist = fopen(distanceFilename,"w");
            }
            
            NL(fl,dummy);
            NL(fl,dummy);
            NL(fl,dummy);

            fscanf(fl,"# sln=%d , maxbrack=%d",&sln,&nRoots); NL(fl,dummy);
            //printf("(a,b) = (%Lf,%Lf) # sln=%d , maxbrack=%d\n",a0[i],b0[i],sln,nRoots);
            if (nRoots>2) {
                for(j=0;j<nRoots;j++){
                    NL(fl,dummy);
                }
                printf("There were more than 2 roots for (a,b)=(%Lf,%Lf) m/n = %d/%d\n",a0[i],b0[i],m,n);
                break;
            }

            nOrbits[i] += nRoots;
            for(j=0;j<nRoots;j++)
            {
            //    NL(fl,dummy);
                fscanf(fl,"%La (%Le) %La (%Le)",&x0[i][j],&dum1,&y0[i][j],&dum2);   NL(fl,dummy);
                distance[i] += fabsl(y0[i][j] - shinoy[i]);
            //     printf("Root was (x,y)=(%4Lf,%4Lf)\n",x0[i][j],y0[i][j]);
            //     updown[i][j] = (y0[i][j]>shinoy[i]) ? 1 : 0;
            //     //Now to find the winding number as we move up and down.
            //     //Most of the calculation is contained within ResonanceWidths
            }
            
            NL(fl,dummy);
          }
            printf("There were %d orbits\n",nOrbits[i]);
            distance[i] = distance[i]/nOrbits[i];
            printf("The average distance from the shearless curve was %Lf\n",distance[i]);
            fprintf(fld,"%21.17Lf %d\n",b0[i],nOrbits[i]);
        }
    


    fclose(fl);
    fclose(fld);

    
    return 0;
}

void ntmapstep(long double *xxx, long double *yyy, long double aaa, long double bbb)
{
    (*yyy)-=bbb*sinl(2.0L*pi*(*xxx));
    (*xxx)+=aaa-aaa*(*yyy)*(*yyy);
}

void standardmapstep(long double *xxx, long double *yyy, long double aaa, long double kkk)
{
    (*yyy)+=kkk*sinl(2.0L*pi*(*xxx));
    (*xxx)+=(*yyy);
}

/* assign shinohara points */
void assignshino(long double a0, long double b0, int shinoNum, long double *xval, long double *yval){
    switch (shinoNum) {
       case 1: {*xval=-0.25L; *yval=-0.5L*b0; break;}
       case 2: {*xval=a0/2.0L-0.25L; *yval=0.0L; break;}
       case 3: {*xval=0.25L; *yval=0.5L*b0;break;}
       case 4: {*xval=a0/2.0L-0.25L; *yval=0.0L; break;}
       default: {*xval=-0.25L; *yval=-0.5L*b0; printf("Uhoh - shinopoint...\n");}
}

*xval=remainderl(*xval,1.0L);
}
// 