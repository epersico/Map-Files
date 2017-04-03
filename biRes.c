#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>
#include <time.h>

#define PI 3.14159265358979323846264338327950288
#define MAX_ORBIT 12
#define MAX_ORBIT_RES 26
#define LARGE_RESIDUE 60.0
#define SYM1 0.0
#define SYM2 0.5
#define SYM3(at,yt) 0.5*(at)*(1.0 - ((yt)*(yt)))
#define SYM4(at,yt) 0.5*(at)*(1.0 - ((yt)*(yt)))+0.5
#define SMOOTH(fs) (fs)/(1.0+fabsl(fs))
#define ERROR 12345678
#define TOLL 1.0e-16
#define PER_RESIDUE 3.0

long q[MAX_ORBIT_RES];
long p[MAX_ORBIT_RES];

long double bmin = 0.0, bmax = 2.0;

long double minmax = -1.0; //-1 if a max, 1 if a min

long double zbrentY(long double (*func)(long double, long double, long double, long double, int),
					long double x1, long double x2, long double tol,
					long double param1, long double param2, long double param3, int param4);

long double zbrentB(long double (*func)(long double, long double, long double, long double, int),
					long double x1, long double x2, long double tol,
					long double param1, long double param2, long double param3, int param4);

void mnbracY(long double *ax, long double *bx, long double *cx,
			 long double (*func)(long double, long double, long double, long double, int),
			 long double param1, long double param2, long double param3, int param4);

void mnbracB(long double *ax, long double *bx, long double *cx,
			 long double (*func)(long double, long double, long double, long double, int),
			 long double param1, long double param2, long double param3, int param4);

long double brentY(long double ax, long double bx, long double cx,
				   long double (*f)(long double, long double, long double, long double, int),
				   long double tol, long double *xmin, long double param1, long double param2, long double param3, int param4);

long double brentB(long double ax, long double bx, long double cx,
				   long double (*f)(long double, long double, long double, long double, int),
				   long double tol, long double *xmin, long double param1, long double param2, long double param3, int param4);

void zbrakY(long double (*fx)(long double, long double, long double, long double, int), long double x1, long double x2, int n,
			long double xb1[], long double xb2[], int *nb, long double param1, long double param2, long double param3, int param4);

void zbrakB(long double (*fx)(long double, long double, long double, long double, int), long double x1, long double x2, int n,
			long double xb1[], long double xb2[], int *nb, long double param1, long double param2, long double param3, int param4);

int zbracY(long double (*func)(long double, long double, long double, long double, int), long double *x1, long double *x2,
		   long double param1, long double param2, long double param3, int param4, long double);

int zbracB(long double (*func)(long double, long double, long double, long double, int), long double *x1, long double *x2,
		   long double param1, long double param2, long double param3, int param4);

long double sol[MAX_ORBIT][2]; // 0 is the y value 1 is the b value
long double oldSol[MAX_ORBIT][2];

long double residueSol[MAX_ORBIT_RES][2]; //0 is down 1 is up
long double residueSolOld[MAX_ORBIT_RES][2]; //0 is down 1 is up
long double residueOrbit[MAX_ORBIT_RES][2];//0 is down 1 is up
long double errorResidue[MAX_ORBIT_RES][2]; //0 is down 1 is up

FILE *data, *resDataUp, *resDataDown, *breakUp, *resDataY;

clock_t endclock, startclock;

long double (*pfunc)(long double, long double, long double, long double, int);
//y, a, b, c, orbit number


void iter(long double a, long double b, long double c, long double *x, long double *y, long n){
	long count;
	for(count = 0; count<n; count++){
		*y = (*y)-b*sinl(2*PI*(*x))-c*sinl(6*PI*(*x));
		*x = (*x)+a*(1-(*y)*(*y));
	}
}

//Symmetry line 1 q odd and p even
long double funcs1OE(long double y, long double a, long double b, long double c, int n){
	long double x;
	x = SYM1;
	iter(a,b,c,&x,&y,p[n]/2);
	//Do R transform
	x-=(q[n]-1)/2;
	//subtract out symmetry line s2
	x-=SYM2;
	x=SMOOTH(x);
	return x;
}

//Symmetry line 1 q odd and p odd
long double funcs1OO(long double y, long double a, long double b, long double c, int n){
	long double x;
	x = SYM1;
	iter(a,b,c,&x,&y,(p[n]+1)/2);
	//Do R transform
	x-=(q[n]-1)/2;
	//subtract out symmetry line s4
	x-=SYM4(a,y);
	x=SMOOTH(x);
	return x;
}

//Symmetry line 1 q even and p odd
long double funcs1EO(long double y, long double a, long double b, long double c, int n){
	long double x;
	x = SYM1;
	iter(a,b,c,&x,&y,(p[n]+1)/2);
	//Do R transform
	x-=q[n]/2;
	//subtract out symmetry line s3
	x-=SYM3(a,y);
	x=SMOOTH(x);
	return x;
}

//Symmetry line 2 q odd and p even
long double funcs2OE(long double y, long double a, long double b, long double c, int n){
	long double x;
	x = SYM2;
	iter(a,b,c,&x,&y,p[n]/2);
	//Do R transform
	x-=(q[n]+1)/2;
	//subtract out symmetry line s1
	x-=SYM1;
	x=SMOOTH(x);
	return x;
}

//Symmetry line 2 q odd and p odd
long double funcs2OO(long double y, long double a, long double b, long double c, int n){
	long double x;
	x = SYM2;
	iter(a,b,c,&x,&y,(p[n]+1)/2);
	//Do R transform
	x-=(q[n]+1)/2;
	//subtract out symmetry line s3
	x-=SYM3(a,y);
	x=SMOOTH(x);
	return x;
}

//Symmetry line 2 q even and p odd
long double funcs2EO(long double y, long double a, long double b, long double c, int n){
	long double x;
	x = SYM2;
	iter(a,b,c,&x,&y,(p[n]+1)/2);
	//Do R transform
	x-=q[n]/2;
	//subtract out symmetry line s4
	x-=SYM4(a,y);
	x=SMOOTH(x);
	return x;
}

//Symmetry line 3 q odd and p even
long double funcs3OE(long double y, long double a, long double b, long double c, int n){
	long double x;
	x = SYM3(a,y);
	iter(a,b,c,&x,&y,p[n]/2);
	//Do R transform
	x-=(q[n]-1)/2;
	//subtract out symmetry line s4
	x-=SYM4(a,y);
	x=SMOOTH(x);
	return x;
}

//Symmetry line 3 q odd and p odd
long double funcs3OO(long double y, long double a, long double b, long double c, int n){
	long double x;
	x = SYM3(a,y);
	iter(a,b,c,&x,&y,(p[n]-1)/2);
	//Do R transform
	x-=(q[n]-1)/2;
	//subtract out symmetry line s2
	x-=SYM2;
	x=SMOOTH(x);
	return x;
}

//Symmetry line 3 q even and p odd
long double funcs3EO(long double y, long double a, long double b, long double c, int n){
	long double x;
	x = SYM3(a,y);
	iter(a,b,c,&x,&y,(p[n]-1)/2);
	//Do R transform
	x-=q[n]/2;
	//subtract out symmetry line s1
	x-=SYM1;
	x=SMOOTH(x);
	return x;
}

//Symmetry line 4 q odd and p even
long double funcs4OE(long double y, long double a, long double b, long double c, int n){
	long double x;
	x = SYM4(a,y);
	iter(a,b,c,&x,&y,p[n]/2);
	//Do R transform
	x-=(q[n]+1)/2;
	//subtract out symmetry line s3
	x-=SYM3(a,y);
	x=SMOOTH(x);
	return x;
}

//Symmetry line 4 q odd and p odd
long double funcs4OO(long double y, long double a, long double b, long double c, int n){
	long double x;
	x = SYM4(a,y);
	iter(a,b,c,&x,&y,(p[n]-1)/2);
	//Do R transform
	x-=(q[n]+1)/2;
	//subtract out symmetry line s1
	x-=SYM1;
	x=SMOOTH(x);
	return x;
}

//Symmetry line 4 q even and p odd
long double funcs4EO(long double y, long double a, long double b, long double c, int n){
	long double x;
	x = SYM4(a,y);
	iter(a,b,c,&x,&y,(p[n]-1)/2);
	//Do R transform
	x-=q[n]/2;
	//subtract out symmetry line s2
	x-=SYM2;
	x=SMOOTH(x);
	return x;
}

//choose the correct function based on the symmetry line and the orbit number
void chooseFunc(int symline, int orbitN){
	switch(symline){
		case 1:
			if(q[orbitN]%2==1 && p[orbitN]%2==0) pfunc = &funcs1OE;
			if(q[orbitN]%2==1 && p[orbitN]%2==1) pfunc = &funcs1OO;
			if(q[orbitN]%2==0 && p[orbitN]%2==1) pfunc = &funcs1EO;
			break;
		case 2:
			if(q[orbitN]%2==1 && p[orbitN]%2==0) pfunc = &funcs2OE;
			if(q[orbitN]%2==1 && p[orbitN]%2==1) pfunc = &funcs2OO;
			if(q[orbitN]%2==0 && p[orbitN]%2==1) pfunc = &funcs2EO;
			break;
		case 3:
			if(q[orbitN]%2==1 && p[orbitN]%2==0) pfunc = &funcs3OE;
			if(q[orbitN]%2==1 && p[orbitN]%2==1) pfunc = &funcs3OO;
			if(q[orbitN]%2==0 && p[orbitN]%2==1) pfunc = &funcs3EO;
			break;
		case 4:
			if(q[orbitN]%2==1 && p[orbitN]%2==0) pfunc = &funcs4OE;
			if(q[orbitN]%2==1 && p[orbitN]%2==1) pfunc = &funcs4OO;
			if(q[orbitN]%2==0 && p[orbitN]%2==1) pfunc = &funcs4EO;
			break;
	}
}

void init(){
	int n;
	q[0]=1;
	p[0]=2;
	q[1]=2;
	p[1]=3;
	for(n=2; n<MAX_ORBIT_RES; n++){
		q[n] = q[n-1]+q[n-2];
		p[n] = p[n-1]+p[n-2];
	}
	//for(n=0; n<MAX_ORBIT; n++) fprintf(data,"%4d %d/%d\n",n,q[n],p[n]);

	for(n=0; n<MAX_ORBIT_RES; n++){
		residueSol[n][0] = 0.0;
		residueSol[n][1] = 0.0;
		residueSolOld[n][0] = 0.0;
		residueSolOld[n][1] = 0.0;
	}
}


long double funcyb(long double ymiddle, long double bracsize, long double a, long double b, long double c, int orb, int *err){
	//find the zeros of a fixed a and b => correspond to y value at which the orbit exists
	long double ybound, ybound1;
	int numroots=2;
	long double ybrakmin[3], ybrakmax[3];
	long double tbracsize;
	int good;
	*err = 0;
	if (bracsize<0.0){
		ybound = ymiddle + bracsize;
		ybound1 = ymiddle; 
	}else{
		ybound = ymiddle;
		ybound1 = ymiddle + bracsize;
	}
	if((*pfunc)(ybound,a,b,c,orb)*(*pfunc)(ybound1,a,b,c,orb)>0.0){
		printf("try growing y bounds out\n");
		good=zbracY(*pfunc,&ybound, &ybound1,a,b,c, orb, bracsize/10.0);
		printf("%Lf %Lf\n", ybound, ybound1);
		if(good==0){
			printf("y bounds too large\n");
			tbracsize= bracsize/10.0;
			if (tbracsize<0.0){
				ybound = ymiddle + tbracsize;
				ybound1 = ymiddle; 
			}else{
				ybound = ymiddle;
				ybound1 = ymiddle + tbracsize;
			}
			good=zbracY(*pfunc,&ybound, &ybound1,a,b,c, orb, bracsize/100);
			if(good == 0){ 
				ybound = ymiddle - fabsl(2*bracsize);
				ybound1 = ymiddle + fabsl(2*bracsize); 
				zbrakY(pfunc,ybound,ybound1,1000,ybrakmin, ybrakmax,&numroots,a,b,c,orb); 
				if(numroots == 2){
					if(bracsize<0){
						printf("returning smaller root");
						return zbrentY(*pfunc, ybrakmin[1], ybrakmax[1], TOLL,a,b,c,orb);
					}else{
						printf("returning larger root");
						return zbrentY(*pfunc, ybrakmin[2], ybrakmax[2], TOLL,a,b,c,orb);
					}
				}else{
					printf("extra tries did not work\n");
					printf("%d numroots\n", numroots);
					*err = 1;
					return 0;
				}
			}
		}
	}
	return zbrentY(*pfunc, ybound, ybound1, TOLL,a,b,c,orb);
}

long double funcby(long double y, long double a, long double b, long double c, int orb){
	//find the zeros of a fixed a and y => correspond to b value at which the orbit exists
	int numRoots=2;
	int good;
	int intervals = 50;
	long double btmin, btmax;
	long double bbrakmin[3], bbrakmax[3];
	if((*pfunc)(y,a,bmin,c,orb)*(*pfunc)(y,a,bmax,c,orb)>0.0){
		zbrakB(*pfunc,bmin,bmax,intervals, bbrakmin,bbrakmax,&numRoots,y,a,c,orb);
		printf("b not bounded try zbrak %d\n", numRoots);
		if(numRoots>0){
			return minmax*zbrentB(*pfunc, bbrakmin[numRoots], bbrakmax[numRoots], TOLL,y,a,c,orb);;
		}else{			
			btmin = bmin; btmax = bmax;
			good = zbracB(*pfunc,&btmin, &btmax, y,a,c,orb);
			if(good == 1){
				bmin = btmin;
				bmax = btmax;
				return minmax*zbrentB(*pfunc, btmin, btmax, TOLL,y,a,c,orb);
			}else{
				intervals = 500;
				numRoots = 2;
				zbrakB(*pfunc,bmin,bmax,intervals, bbrakmin,bbrakmax,&numRoots,y,a,c,orb);			
				printf("number of roots %d\n",numRoots);
				if(numRoots>0){
					return minmax*zbrentB(*pfunc, bbrakmin[numRoots], bbrakmax[numRoots], TOLL,y,a,c,orb);
				}else{
					printf("could not bracket b\n");
					exit(EXIT_FAILURE);
				}
			}
		}
	}else{
		return minmax*zbrentB(*pfunc, bmin, bmax, TOLL,y,a,c,orb);
	}
}

long double characterize(long double a, long double b, long double c, long double y, int orb, int symline){
	long double x;
	long double dx1dx, dx1dy, dy1dx, dy1dy, y1;
	long double M11=1, M12=0, M21=0, M22=1;
	long double M11T=1, M12T=0, M21T=0, M22T=1;
	long double resid;
	int count=0;
	switch(symline){
		case 1:
			x = SYM1;
			break;
		case 2:
			x = SYM2;
			break;
		case 3:
			x = SYM3(a,y);
			break;
		case 4:
			x = SYM4(a,y);
			break;
                default:
                        x = SYM1;
                        break;
	};
	for(count = 0; count<p[orb]; count++){

		y1 = y - b*sinl(2*PI*x)-c*sinl(6*PI*x);
		dy1dx = -2.0*PI*b*cosl(2*PI*x)-6.0*PI*c*cosl(6*PI*x);
		dy1dy = 1;
		dx1dx = 1 - 2*a*y1*dy1dx;
		dx1dy = -2*a*y1;
		M11T = dx1dx*M11+dx1dy*M21;
		M12T = dx1dx*M12+dx1dy*M22;
		M21T = dy1dx*M11+dy1dy*M21;
		M22T = dy1dx*M12+dy1dy*M22;
		M11=M11T;
		M12=M12T;
		M21=M21T;
		M22=M22T;
		iter(a,b,c,&x,&y,1);
		//while(x>.5) x = x-1.0;
		//while(x<-.5) x = x+1.0;
		//printf("%f %f\n%f %f\n\n",M11,M12,M21,M22);
	}
	resid = .25*(2-(M11+M22));
	/*
	printf("%.12f\n", residue);
	if(residue==0||residue==1){
	printf("parabolic orbit\n");
	}else if(0<residue && residue<1){
	printf("elliptic orbit\n");
	}else if(residue<0||residue>1){
	printf("hyperbolic orbit\n");
	}
	*/
	return resid;
}

long double errorOrbit(long double a, long double b, long double c, long double y, int orb, int symline){
	long double x;
	switch(symline){
		case 1:
			x = SYM1;
			break;
		case 2:
			x = SYM2;
			break;
		case 3:
			x = SYM3(a,y);
			break;
		case 4:
			x = SYM4(a,y);
			break;
                default:
                        x = SYM1;
                        break;
	};
	iter(a,b,c,&x,&y,p[orb]);
	x-=q[orb];
	switch(symline){
		case 1:
			x -= SYM1;
			break;
		case 2:
			x -= SYM2;
			break;
		case 3:
			x -= SYM3(a,y);
			break;
		case 4:
			x -= SYM4(a,y);
			break;
	};
	return x;
}

void initFiles(char *nme){
	int n;
	char orb[100] ="orb_";
	char bre[100] ="breakUp_";
	char resu[100] ="resUp_";
  char resd[100] ="resDown_";
  char resy[100] ="resY_";
	data=fopen(strcat(orb, nme), "w");
	if(data == NULL){
		fprintf(stderr, "Error: could not open file");
		exit(EXIT_FAILURE);
	}

	breakUp=fopen(strcat(bre, nme), "w");
	if(breakUp == NULL){
		fprintf(stderr, "Error: could not open file");
		exit(EXIT_FAILURE);
	}

	resDataUp = fopen(strcat(resu, nme), "w");
	if(resDataUp == NULL){
		fprintf(stderr, "Could not open file");
		exit(EXIT_FAILURE);
	}

	resDataDown = fopen(strcat(resd, nme), "w");
	if(resDataDown == NULL){
		fprintf(stderr, "Could not open file");
		exit(EXIT_FAILURE);
	}

	resDataY = fopen(strcat(resy, nme), "w");
	if(resDataY == NULL){
		fprintf(stderr, "Could not open file");
		exit(EXIT_FAILURE);
	}

	//initialize breakUp file
	fprintf(breakUp,"#%15s %15s %15s %15s\n", "a","b","c","y");
	//initialize data file
	fprintf(data,"#            a");
	for(n=2; n<MAX_ORBIT; n++){
		fprintf(data, " %14dy %13db",n+1,n+1);
	}
	fprintf(data, "\n");}

int main(int argc, char *argv[])
{
	int err = 0;
	int n;
	int goOn=0;
	int innercurve;
	int symLineNumber=1;
	long double radiusy;
	//long double y = 0.0;
	long double ymin=-0.3, ymax=0.2;
	long double ymint, ymaxt;
	long double a=0.7, b=0.265, c = 0.00;
	long double ainit, afinal, ainc;
	long double bound1y, bound2y, bound3y;
	long double rootRadius;
	long double yd, yu;
  long double bminS=0.0, bmaxS=2.0;
  long double scaledB;
	char fname[100] = "test.txt";
	ainit = 0.66;
	afinal = 0.75;
	ainc = 1.0e-2;
	bound1y = -0.3;
	bound2y = -0.2;
	bound3y = 0.3;
	switch(argc){
    case 10:
			sscanf(argv[9], "%Lf", &bmaxS);
			sscanf(argv[8], "%Lf", &bminS);
		case 8:
			sscanf(argv[7], "%Lf", &ymax);
			sscanf(argv[6], "%Lf", &ymin);
		case 6:
			sscanf(argv[5], "%Lf", &ainc);		
		case 5:
			sscanf(argv[4], "%Lf", &ainit);
		case 4:
			sscanf(argv[3], "%d", &symLineNumber);
		case 3:
			sscanf(argv[2], "%Lf", &c);
		case 2:
			strcpy(fname, argv[1]);
			break;
		case 0:
		case 1:
			break;
		default:
			fprintf(stderr, "Incorrect Usage");
			exit(EXIT_FAILURE);
			break;
	}
	startclock = clock();
	init();

	initFiles(fname);
  fprintf(resDataUp, "# c = %.12Lf\n", c);
  fprintf(resDataDown, "# c = %.12Lf\n", c);

	for(a=ainit; a<afinal; a+=ainc){
    printf("\na = %.12Lf\n\n", a);
		//find the first two orbits
		bmin=bminS;
		bmax=bmaxS;
		n=2;
    printf("obit 2\n");
		chooseFunc(symLineNumber, n);
		bound1y = ymin;
		bound3y = ymax;
		bound2y = -0.1;
		//mnbracY(&bound1y,&bound2y,&bound3y, funcby, a, b, c, n);
		sol[n][1]=minmax*brentY(bound1y,bound2y,bound3y,funcby,TOLL,&sol[n][0],a,b,c,n);
		//fprintf(data,"%.12f ", errorOrbit(a,sol[n][1],c,sol[n][0],n,symLineNumber));					

    printf("obit 3\n");
		n=3;
		chooseFunc(symLineNumber, n);
		bound1y = ymin;
		bound3y = ymax;
		bound2y = -0.1;
		//mnbracY(&bound1y,&bound2y,&bound3y, funcby, a, b, c, n);
		sol[n][1]=minmax*brentY(bound1y,bound2y,bound3y,funcby,TOLL,&sol[n][0],a,b,c,n);
		//fprintf(data,"%.12f ", errorOrbit(a,sol[n][1],c,sol[n][0],n,symLineNumber));					

		//determine which curve is the inner curve
		//expand outward from the inner curve to the outer curve
		//finding the bounds on y

		if(sol[2][1]>sol[3][1]){
			//curve 3 is the inner curve
			innercurve = 1; //inner curve always odd
			chooseFunc(symLineNumber,2);
			ymaxt = funcyb(sol[3][0],0.4,a,sol[3][1],c,2, &err);
			if(err == 1 && ainit == a && ainit>=.4) {
				printf("changing ainit\n");
				a=ainit -= .01;
				a-=ainc;
				continue;
			}else if(err ==1){
				fprintf(breakUp, "never could find good starting point.\n");
				break;
			}
			ymint = funcyb(sol[3][0],-0.4,a,sol[3][1],c,2, &err);
			if(err == 1 && ainit == a && ainit>=.4) {
				printf("changing ainit\n");
				a=ainit -= .01;
				a-=ainc;
				continue;
			}else if(err == 1){
				fprintf(breakUp, "never could find good starting point.\n");
				break;
			}
			bmin = sol[3][1];
			bmax = sol[2][1];
		}else{
			//curve 2 is the inner curve
			innercurve = 0; //inner curve always even
			chooseFunc(symLineNumber,3);
			ymaxt= funcyb(sol[2][0],0.4,a,sol[2][1],c,3, &err);
			if(err == 1 && ainit == a && ainit>=.4) {
				printf("changing ainit\n");
				a=ainit -= .01;
				a-=ainc;
				continue;
			}else if(err==1){
				fprintf(breakUp, "never could find good starting point.\n");
				break;
			}
			ymint=funcyb(sol[2][0],-0.4,a,sol[2][1],c,3, &err);
			if(err == 1 && ainit == a && ainit>=.4) {
				printf("changing ainit\n");
				a=ainit -= .01;
				a-=ainc;
				continue;
			}else if(err==1){
				fprintf(breakUp, "never could find good starting point.\n");
				break;
			}
			bmin = sol[2][1];
			bmax = sol[3][1];
		}
    printf("obit 4 through n\n");
		if(err==0){
			for(n = 4; n<MAX_ORBIT; n++){
        printf("orbit %d\n", n);
				chooseFunc(symLineNumber, n);
				bound1y = ymint;
				bound3y = ymaxt;
				bound2y = sol[n-1][0];
				//mnbracY(&bound1y,&bound2y,&bound3y, funcby, a, b, c, n);
				sol[n][1]=minmax*brentY(bound1y,bound2y,bound3y,funcby,TOLL,&sol[n][0],a,b,c,n);
				//fprintf(data,"%.12f ", errorOrbit(a,sol[n][1],c,sol[n][0],n,symLineNumber));					

				if(n%2 == innercurve){
					//calculate new y bounds by expanding from innercurve to outer curve
					chooseFunc(symLineNumber,n+1);
					radiusy = (ymaxt-ymint);
					ymaxt= funcyb(sol[n][0],radiusy,a,sol[n][1],c,n+1, &err);
					if(err==1) break;
					ymint=funcyb(sol[n][0],-radiusy,a,sol[n][1],c,n+1, &err);
					if(err==1) break;
					bmin = sol[n][1];
				}else{
					bmax = sol[n][1];
				}
				if(err==1) break;
			}
		}
		//fprintf(data, "\n");

		if(err==0){
/*      scaledB = (sol[MAX_ORBIT-13][1]*sol[MAX_ORBIT-2][1]-sol[MAX_ORBIT-14][1]*sol[MAX_ORBIT-1][1])/
      ((sol[MAX_ORBIT-13][1]-sol[MAX_ORBIT-14][1])-(sol[MAX_ORBIT-1][1]-sol[MAX_ORBIT-2][1]));
*/
      scaledB = sol[MAX_ORBIT-1][1];
			rootRadius = (sqrtl(1.0-(q[n]/(p[n]*a))));
			for(n=0; n<MAX_ORBIT_RES; n=n+2){
				chooseFunc(symLineNumber, n);
        if(MAX_ORBIT>=n){
				  //find down residue
				  yd = funcyb(sol[MAX_ORBIT-1][0],-rootRadius,a,scaledB,c,n, &err);
          residueOrbit[n][0]=yd;			  
          if(err==1){ printf("problem in finding down residue\n"); break;}
				  residueSol[n][0] = characterize(a,scaledB,c,yd,n,symLineNumber);
				  errorResidue[n][0]=errorOrbit(a,scaledB,c,yd,n,symLineNumber);					

				  //find up residue
				  yu = funcyb(sol[MAX_ORBIT-1][0],rootRadius,a,scaledB,c,n, &err);
          residueOrbit[n][1]=yu;				  
          if(err==1){ printf("problem in finding up residue\n"); break;}
				  residueSol[n][1] = characterize(a,scaledB,c,yu,n,symLineNumber);
				  errorResidue[n][1]=errorOrbit(a,scaledB,c,yu,n,symLineNumber);
          rootRadius=(yu-yd)/1.5;
        }else{
				  //find down residue
				  yd = funcyb(sol[MAX_ORBIT-1][0],-rootRadius,a,scaledB,c,n, &err);
          residueOrbit[n][0]=yd;			  
          if(err==1){ printf("problem in finding down residue\n"); break;}
				  residueSol[n][0] = characterize(a,scaledB,c,yd,n,symLineNumber);
				  errorResidue[n][0]=errorOrbit(a,scaledB,c,yd,n,symLineNumber);					

				  //find up residue
				  yu = funcyb(sol[MAX_ORBIT-1][0],rootRadius,a,scaledB,c,n, &err);
          residueOrbit[n][1]=yu;				  
          if(err==1){ printf("problem in finding up residue\n"); break;}
				  residueSol[n][1] = characterize(a,scaledB,c,yu,n,symLineNumber);
				  errorResidue[n][1]=errorOrbit(a,scaledB,c,yu,n,symLineNumber);
          rootRadius=(yu-yd)/1.5;          
        }
					
				if(((fabsl(residueSol[n][0])>.3 && fabsl(residueSolOld[n][0])>.005) 
					&& fabsl((residueSol[n][0] - residueSolOld[n][0])/ residueSolOld[n][0]) > PER_RESIDUE ) ||
					((fabsl(residueSol[n][1])>.3 && fabsl(residueSolOld[n][1])>.005) 
					&& fabsl((residueSol[n][1] - residueSolOld[n][1])/ residueSolOld[n][1]) > PER_RESIDUE)){
							printf("\nUsing Percentage Criterion\n\n");
							printf("%d %.12Lf %.12Lf\n\n",n, residueSol[n][0], residueSolOld[n][0]);
							printf("%d %.12Lf %.12Lf\n\n",n, residueSol[n][1], residueSolOld[n][1]);
							goOn =1;
							//break;
				}else if(n>20 && fabsl(residueSol[n][0])>.25 && fabsl(residueSol[n][1])>.25){
          if(fabsl(residueSol[n][0])>fabsl(2.0*residueSol[n-12][0]) || fabsl(residueSol[n][1])>fabsl(2.0*residueSol[n-12][1])){
							printf("\nUsing Large Residue\n\n");
							printf("%d down %.12Lf %.12Lf\n\n",n+1, residueSol[n][0], residueSol[n-12][0]);
							printf("%d up %.12Lf %.12Lf\n\n",n+1, residueSol[n][1], residueSol[n-12][1]);
							goOn =1;
							//break;
				  }
        }
			}
			//fprintf(resData,"\n");
		}

		//print data

		if(err==0){
			if(goOn == 0){
				fprintf(data, "%.12Lf ", a);
				fprintf(resDataUp, "# a = %.12Lf b = %.12Lf y = %.12Lf scaledB = %.12Lf\n", a, sol[MAX_ORBIT-1][1], sol[MAX_ORBIT-1][0], scaledB);
				fprintf(resDataDown, "# a = %.12Lf b = %.12Lf y = %.12Lf scaledB = %.12Lf\n", a, sol[MAX_ORBIT-1][1], sol[MAX_ORBIT-1][0], scaledB);
				fprintf(resDataY, "# a = %.12Lf b = %.12Lf y = %.12Lf scaledB = %.12Lf\n", a, sol[MAX_ORBIT-1][1], sol[MAX_ORBIT-1][0], scaledB);
				for(n = 2; n<MAX_ORBIT; n++){
					oldSol[n][0]=sol[n][0];
					oldSol[n][1]=sol[n][1];
				}
				for(n=0; n<MAX_ORBIT_RES; n++){
					residueSolOld[n][0] = residueSol[n][0];
					residueSolOld[n][1] = residueSol[n][1];
				}
			}else{
				printf("goOn = 1 changed ainc\n");
				fprintf(data, "#%.11Lf ", a);
				fprintf(resDataUp, "## a = %.12Lf b = %.12Lf y = %.12Lf scaledB = %.12Lf\n", a, sol[MAX_ORBIT-1][1], sol[MAX_ORBIT-1][0], scaledB);
				fprintf(resDataDown, "## a = %.12Lf b = %.12Lf y = %.12Lf scaledB = %.12Lf\n", a, sol[MAX_ORBIT-1][1], sol[MAX_ORBIT-1][0], scaledB);
				fprintf(resDataY, "## a = %.12Lf b = %.12Lf y = %.12Lf scaledB = %.12Lf\n", a, sol[MAX_ORBIT-1][1], sol[MAX_ORBIT-1][0], scaledB);
				a-=ainc;
				ainc/=10.0;
				goOn = 0;
			}
			for(n=2; n<MAX_ORBIT; n++) fprintf(data, "%.12Lf %.12Lf ", sol[n][0], sol[n][1]);
			fprintf(data,"\n");

			for(n=0; n<MAX_ORBIT_RES; n=n+2){
//         fprintf(resDataUp, "%3d %16.12Lf\n", n+1, residueSol[n][1]);
//         fprintf(resDataDown, "%3d %16.12Lf\n", n+1, residueSol[n][0]);
         fprintf(resDataUp, "%3d %16.12Lf %16.12Lf\n", n+1, residueSol[n][1], errorResidue[n][1]);
         fprintf(resDataDown, "%3d %16.12Lf %16.12Lf\n", n+1, residueSol[n][0], errorResidue[n][0]);
         fprintf(resDataY, "%3d %16.12Lf %16.12Lf\n", n+1, residueOrbit[n][0], residueOrbit[n][1]);
      }
		}else{
			printf("error level set\n");
			a-=ainc;
			ainc/=10.0;
			err=0;
			goOn=0;
		}
		if(ainc<1.0e-9) break;


	}
	printf("error is set at %d\n", err);
	printf("exiting normally\n");
	printf("ainc is %.12Lf a is %.12Lf\n", ainc, a);
	fprintf(breakUp,"%15.12Lf %15.12Lf %15.12Lf %15.12Lf\n", a, oldSol[MAX_ORBIT-1][1], c, oldSol[MAX_ORBIT-1][0]);
	endclock = clock ();
    printf("time taken %9.2f seconds\n",
            (float)(endclock-startclock)/(float)CLOCKS_PER_SEC );
	fprintf(breakUp, "#time taken %9.2f seconds\n",
            (float)(endclock-startclock)/(float)CLOCKS_PER_SEC );
	fclose(breakUp);
	fclose(data);
	fclose(resDataUp);
  fclose(resDataDown);
	exit(EXIT_SUCCESS);
}
