#include <math.h>

main()
{
  long double pi=4.0L*atanl(1.0L);
  double shortpi=4.0*atan(1.0);
  long double pi2=4.0L*atanl(1.0L)+1.0e-15L;
  double shortpi2=4.0*atan(1.0)+1.0e-15L;
  long double pi3=4.0L*atanl(1.0L)+1.5e-19L;
  double shortpi3=4.0*atan(1.0)+1.5e-19L;
  
  printf("%Lf %f %Le \n",pi, shortpi, pi-shortpi);
  printf("%Lf %f %Le \n",pi2, shortpi2, pi2-shortpi2);
  printf("%Le %e %Le \n",pi-pi2, shortpi-shortpi2, 
	(pi-shortpi)-(pi2-shortpi2));
  printf("%Lf %f %Le \n",pi3, shortpi3, pi3-shortpi3);
  printf("%Le %e %Le \n",pi-pi3, shortpi-shortpi3, 
	(pi-shortpi)-(pi3-shortpi3));
}
