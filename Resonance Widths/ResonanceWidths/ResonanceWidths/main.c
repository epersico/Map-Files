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
#define nParameters 1000

void ntmapstep(long double *, long double *, long double, long double);
void find_om(long double, long double, long double, long double, long double *, long double *,
             int *, long double *, int *,
             long double *, int *, long double *, int *);
void assignshino(long double, long double, int, long double *, long double *);
void findShinoOrbitYonSLN(long double, long double, long double *,long double *);
void ResonanceWidths(long double, long double, long double, long double,int, int, int, long double *);


long double  pi,cutoff,dn,speedup;
int nMax;
int n_ysteps=5000;
long double yRange = 3e-3;
long double windAccuracy = 1.0e-4L;
int n_bsteps=300;
char infilename[64],outfilename[64], widthFilename[60];
FILE *fl,*flw;

int main(int argc, const char * argv[]) {
    // insert code here...
    printf("Hello, World!\n");

    sprintf(infilename,"./data/orbits_3.out");
    char dummy;
    pi=4.0L*atanl(1.0L);
    int n,m;//nParameters=n_bsteps;
    long double a0[nParameters],b0[nParameters];
    long double x0[nParameters][2],y0[nParameters][2],residue[nParameters][2];
    long double widths[nParameters][2];
    int sln, nRoots,i,j,updown[nParameters][2];
    long double shinox,shinoy;
    nMax = 1.0e6L; //Total number of iterations allowed for the winding number to converge
    cutoff = 1.0e-3L; //The winding number needs to stay within this amount to call it converged
    dn=1.0e3L;  //How long it needs to stay the same to be called convergent
    
    //Rigorous values are nMax=1.0e-8L, cutoff - 1.0e-4L, dn =1.0e5L
    printf("The filename opened was %s\n",infilename);
    fl = fopen(infilename,"r");

    for(i=0;i<nParameters;i++){
        fscanf(fl,"# a: %Lf , b: %Lf",&a0[i],&b0[i]);   NL(fl,dummy);
        printf("(a,b) = (%Lf,%Lf)\n",a0[i],b0[i]);
        fscanf(fl, "# omega:   (m =   %d) / (n =  %d) ",&m,&n);  NL(fl,dummy);
        //printf("(a,b)=(%Lf,%Lf)\n",a0[i],b0[i]);
        
        if(i==0){
        sprintf(widthFilename,"./widths/a_%Lf_m_%d_n_%d_widths.out",a0[0],m,n);
        flw = fopen(widthFilename,"w");
        }
        
        //printf("(shinox,shinoy)=(%Lf,%Lf)\n",shinox,shinoy);
        NL(fl,dummy);
        NL(fl,dummy);
        NL(fl,dummy);
        fscanf(fl,"# sln=%d , maxbrack=%d",&sln,&nRoots); NL(fl,dummy);
        //printf("(a,b) = (%Lf,%Lf) # sln=%d , maxbrack=%d\n",a0[i],b0[i],sln,nRoots);
        if (nRoots>2) {
            for(j=0;j<nRoots;j++){
                NL(fl,dummy);
            }
            printf("There were more than 2 roots for (a,b)=(%Lf,%Lf)\n",a0[i],b0[i]);
            break;
        }
        else{
            findShinoOrbitYonSLN(a0[i], b0[i], &shinox, &shinoy);
        }
        for(j=0;j<nRoots;j++)
        {
            fscanf(fl,"%Lf  %Lf %Lf",&x0[i][j],&y0[i][j],&residue[i][j]);   NL(fl,dummy);
            updown[i][j] = (y0[i][j]>shinoy) ? 1 : 0;
            //Now to find the winding number as we move up and down.
            if (residue[i][j]<1 && residue[i][j] > 0)   ResonanceWidths(a0[i],b0[i],x0[i][j],y0[i][j],m,n,updown[i][j],&widths[i][updown[i][j]]);
        }
        
        NL(fl,dummy);

    }
    
    

    
    return 0;
}

void ntmapstep(long double *xxx, long double *yyy, long double aaa, long double bbb)
{
    (*yyy)-=bbb*sinl(2.0L*pi*(*xxx));
    (*xxx)+=aaa-aaa*(*yyy)*(*yyy);
}


void find_om(long double a0, long double b0, long double x0, long double y0, long double *yf, long double *om,
             int *ncut, long double *diffom, int *ndiffom,
             long double *maxom, int *nmaxom, long double *minom, int *nminom)
{
    int i;
    long double omo, x=x0, y=y0;
    
    *ncut=nMax+1; *diffom=(UNDEFINED); *ndiffom=nMax+1;
    *maxom=-1.0e11L; *nmaxom=0;
    *minom=1.0e11L; *nminom=0;
    
    ntmapstep(&x, &y, a0, b0);
    omo = (x-x0);
    
    i=2;
    do           /* iteration of nt-map */
    {
        ntmapstep(&x, &y, a0, b0);
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

//note that this code assumes we are looking at the dominant symmetry line.  To get other information need to change where you are looking to stop near a different x value
void findShinoOrbitYonSLN(long double a0, long double b0, long double *x ,long double *y)
{
    long double accuracy = 1e-3L;
    int i=0;
    assignshino(a0, b0, 1, x, y);
    //iterate shino until x is close to 0
    while((fabsl(*x) > accuracy) && fabsl(1-*x) > accuracy){
        ntmapstep(x, y, a0, b0);
        *x = remainderl(*x, 1.0L);
        i++;
        if(i>1e4) {
            printf("Shino broke!\n");
            break;
        }
    }
}

void ResonanceWidths(long double a, long double b, long double x0, long double y0, int m, int n, int updown,long double *width){
    int i,ncutoff,ndom,nomin,nomax,error=0,widthCount=0;
    long double stepsize = yRange/n_ysteps,yf,omega[2*n_ysteps][2],omega0,dom,omax,omin,y=y0-stepsize*n_ysteps;
    char filename[60]="./widths/",astr[60],bstr[60],mstr[60],nstr[60],updownstr[60];
    FILE *fl;
    
    long double yLow=y, yHigh=y0+stepsize*n_ysteps;
    
    sprintf(astr,"a_%Lf_",a);
    sprintf(bstr,"b_%Lf_",b);
    sprintf(mstr,"m_%d_",m);
    sprintf(nstr,"n_%d_",n);
    sprintf(updownstr,"updown_%d_",updown);
    
    strcat(filename,astr);
    strcat(filename,bstr);
    strcat(filename,mstr);
    strcat(filename,nstr);
    strcat(filename,updownstr);
    strcat(filename, "width.out");
    fl = fopen(filename,"w");
    
    find_om(a, b, x0, y0, &yf, &omega0, &ncutoff, &dom, &ndom, &omax, &nomax, &omin, &nomin);
    
    for(i=0;i<2*n_ysteps;i++){
        y += stepsize;
        find_om(a, b, x0, y, &yf, &omega[i][updown], &ncutoff, &dom, &ndom, &omax, &nomax, &omin, &nomin);
        
        if(omega[i][updown] > -1 && omega[i][updown] < 1){
            fprintf(fl,"%21.17Lf %21.17Lf\n",y, omega[i][updown]);
            if(  fabsl(omega[i][updown]-omega0) < windAccuracy  &&  i<=n_ysteps && error == 0 ){
                yLow = y; error=1; widthCount++;
                //printf("Bottom bound was %Lf\n",yLow);
            }
            if(  fabsl(omega[i][updown]-omega0) < windAccuracy  &&  i>n_ysteps && error != 0){
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



