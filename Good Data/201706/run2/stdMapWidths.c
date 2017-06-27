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
#define nParameters 100

void standardmapstep(long double *, long double * ,long double);
void find_om(long double, long double, long double, long double *, long double *,
             int *, long double *, int *,
             long double *, int *, long double *, int *);
void ResonanceWidths(long double, long double, long double, int, int, long double *);



long double  pi,cutoff,dn,speedup;
int nMax;
int n_ysteps=4000;
long double yRange = .5;
long double windAccuracy = 2e-5;
long double bMax = 0.35;
char outfilename[64], widthFilename[60];
FILE *fl,*flw;

//int n_ysteps=1000;
//long double yRange = 2e-6;
//long double windAccuracy = 2e-7;
//long double bMax = 0.0

int main(int argc, const char * argv[]) {
    pi=4.0L*atanl(1.0L);
    int n,m,i;//nParameters=n_bsteps;
    long double b;
    long double x,y;
    long double widths[nParameters];
    nMax = 5.0e5L; //Total number of iterations allowed for the winding number to converge
    cutoff = 1.0e-4L; //The winding number needs to stay within this amount to call it converged
    dn=1.0e4L;  //How long it needs to stay the same to be called convergent

    x=0.0L; y=0.5L; m=1; n=2;
    for(i=1;i<=nParameters;i++){
        b = i * bMax/nParameters;
        printf("Working on value b = %Lf\n",b);
        if(i==1){
        sprintf(widthFilename,"./widths/m_%d_n_%d_widths.out",m,n);
        flw = fopen(widthFilename,"w");
        }
        ResonanceWidths(b,x,y,m,n,&widths[i-1]);
    }

    fclose(flw);

    return 0;
}

void ntmapstep(long double *xxx, long double *yyy, long double aaa, long double bbb)
{
    (*yyy)-=bbb*sinl(2.0L*pi*(*xxx));
    (*xxx)+=aaa-aaa*(*yyy)*(*yyy);
}

void standardmapstep(long double *xxx, long double *yyy, long double kkk)
{
    (*yyy)+=kkk*sinl(2.0L*pi*(*xxx));
    (*xxx)+=(*yyy);
}

void find_om(long double b0, long double x0, long double y0, long double *yf, long double *om,
             int *ncut, long double *diffom, int *ndiffom,
             long double *maxom, int *nmaxom, long double *minom, int *nminom)
{
    int i;
    long double omo, x=x0, y=y0;
    
    *ncut=nMax+1; *diffom=(UNDEFINED); *ndiffom=nMax+1;
    *maxom=-1.0e11L; *nmaxom=0;
    *minom=1.0e11L; *nminom=0;
    
    standardmapstep(&x, &y,b0);
    omo = (x-x0);
    
    i=2;
    do           /* iteration of nt-map */
    {
        standardmapstep(&x, &y, b0);
        *om = (x-x0)/((long double)(i));
        if (*om>*maxom) {*maxom=*om; *nmaxom=i;}     /* store max. of om */
        if (*om<*minom) {*minom=*om; *nminom=i;}     /* store min. of om */
        if (*ncut<=nMax) {                              /* check if dom remains small */
            if (fabsl(omo-*om)<cutoff) {                /* store max. dom after conv. */
                if (fabsl(omo-*om)>fabsl(*diffom)) {*diffom=*om-omo; *ndiffom=i;}
            }
            else {
                *ncut=nMax+1;                               /* left convergence range */
                if (fabsl(omo-*om)>fabsl(*diffom)) {*diffom=*om-omo; *ndiffom=i;}
            }
        }
        else if (fabsl(omo-*om)<cutoff) {             /* entered convergence range */
            *ncut=i;
            *diffom=*om-omo; *ndiffom=i;
        }
        omo=*om;
        i++;
    }
    while ((i<=nMax)&&((i<=*ncut+dn)||(i<=*nmaxom+dn)||(i<=*nminom+dn)));
    if ((i>nMax) && ((i<=*ncut+dn)||(i<=*nmaxom+dn)||(i<=*nminom+dn))) *om=(UNDEFINED);
    *yf=y;
}


void ResonanceWidths(long double b, long double x0, long double y0, int m, int n,long double *width){
    int i,ncutoff,ndom,nomin,nomax,error=0,widthCount=0;
    long double stepsize = yRange/n_ysteps,yf,omega[2*n_ysteps],omega0,dom,omax,omin,y=y0-stepsize*n_ysteps;
    char filename[60]="./widths/",astr[60],bstr[60],mstr[60],nstr[60],updownstr[60];
    FILE *fl;
    
    long double yLow=y, yHigh=y0+stepsize*n_ysteps;
    
    sprintf(bstr,"b_%Lf_",b);
    sprintf(mstr,"m_%d_",m);
    sprintf(nstr,"n_%d_",n);

    strcat(filename,bstr);
    strcat(filename,mstr);
    strcat(filename,nstr);
    strcat(filename,updownstr);
    strcat(filename, "width.out");
    fl = fopen(filename,"w");
    
    find_om(b, x0, y0, &yf, &omega0, &ncutoff, &dom, &ndom, &omax, &nomax, &omin, &nomin);
    
    for(i=0;i<2*n_ysteps;i++){
        y += stepsize;
        find_om(b, x0, y, &yf, &omega[i], &ncutoff, &dom, &ndom, &omax, &nomax, &omin, &nomin);
        if(fabsl(omega[i] - omega0)< 1){
            fprintf(fl,"%21.17Lf %21.17Lf\n",i*stepsize-stepsize*n_ysteps, omega[i]);
            if(  fabsl(omega[i]-omega0) < windAccuracy  &&  i<=n_ysteps && error == 0 ){
                yLow = y; error=1; widthCount++;
                //printf("Bottom bound was %Lf\n",yLow);
            }
            if(  fabsl(omega[i]-omega0) < windAccuracy  &&  i>n_ysteps && error != 0){
                yHigh = y;
                error+=1;
                //printf("Upper bound was %Lf\n",yHigh);
            }
            *width = fabsl(yHigh-yLow);
            //*width = (long double)stepsize*widthCount;
        }
    }
    if (error >= 2) {
        printf("Width was, %Lf\n",*width);
        fprintf(flw,"%21.17Lf %21.17Lf\n",b,*width);
        error +=1;
    }
    fclose(fl);
}



