#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

int main(int argc, char *argv[])
{
    double *a;
    int test = atoi(argv[1]);
    double divisor = (double)atoi(argv[2]);
    int n = atoi(argv[3]);
    a = (double*)malloc( sizeof(double)*n );

    int i;
    for(i=0;i<n;i++){
        a[i]=(double)rand()%10000
        if(i%(n/10)==0) printf("%f\n",a[i]);
    }
 float t;

    t=clock();
    if (test == 0) {
        for (int i = 0; i < n; i++)
            a[i] = fmod(a[i],divisor);
    } else if (test == 1) {
        for (int i = 0; i < n; i++)
            while (a[i]>divisor) a[i]=a[i]-divisor;
            while (a[i]<0) a[i]=a[i]+divisor;
    }
    t = (clock()-t)/(float)CLOCKS_PER_SEC;


    printf("The calculation took %f seconds",t);
}
