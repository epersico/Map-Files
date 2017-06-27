#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#define NL(flnm,dumch) {do dumch=fgetc(flnm); while (dumch!='\n' && dumch!=EOF);}
#define pi  3.1415926535897932384626433832795028841971693993751058209749445923078164062862089986280348253421170679L

int main()
{
//	long double pi=4.0L*atanl(1.0L);
	char dummy;
	char filename[60]="pitest.txt";
	FILE *fl;

	fl=fopen(filename,"w");

	fprintf(fl,"no specification %Lf\n",pi);
	fprintf(fl,"20.16: %20.16Lf\n",pi);
	fprintf(fl,"30.24: %30.24Lf\n",pi);
	fprintf(fl,"40.34: %40.34Lf\n",pi);
//	fprintf(fl,"a: %a\n",pi);
	fprintf(fl,"La: %La\n",pi);


	fclose(fl);

	fopen(filename,"r");

	long double tmp,tmp1;

	fscanf(fl,"no specification %Lf\n",&tmp); 

	printf("nospecification %Lf\n",tmp);
	printf("nospecification %Lf\n",pi);

	//tmp=0L;

	fscanf(fl,"20.16: %Lf\n",&tmp1); 

	printf("20.16: %20.16Lf\n",tmp1);
	printf("20.16: %20.16Lf\n",pi);

	fscanf(fl,"30.24: %Lf\n",&tmp);

	printf("20.24: %Lf\n",tmp);
	printf("20.24: %Lf\n",pi);

	fscanf(fl,"40.34: %Lf\n",&tmp); 

	printf("40.34 %44.34Lf\n",tmp);
	printf("40.34 %50.40Lf\n",pi);


	fscanf(fl,"La: %La\n",&tmp);

	printf("La: %100.80Lf\n",tmp);
	printf("La: %100.62Lf\n",pi);

	return 0;
}