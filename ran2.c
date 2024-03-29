#include <stdio.h>
#include <stdlib.h>

#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

/*Long period (>2e10^18) random number generator of L'Ecuyer with Bays-Durham shuffle
and added safeguards. Returns a random deviate between 0.0 and 1.0 (exlclusive of the endpoint values).
Call with idum  a negative integer to initialize;
thereafter do not alter idum between successive deviates in a sequence.
RNMX should approxiamte the largest floating value that is less than 1. */

float ran2(long *idum);

/*this bit i entered i doubt it is even close to being correct
i always get the same number, which isnt even in the range 0.0 to 1.0!*/

/*int main()
{

   float answer;
   long *pointer;
   long initial = -0.2;        it says im supposed to initialise with a call to idum, am i doing that here?
   pointer=&initial;


   answer=ran2(pointer);

   printf ("%f", answer);

   printf ("\nPress ENTER to continue.\n");

   fflush (stdin);
   (void) getchar();
   return (0);

}*/
/*end of bit i entered!*/

/*this is the ran2 function listed in the book, straight from the cd*/

float ran2(long *idum)
{
	int j;
	long k;
	static long idum2=123456789;
	static long iy=0;
	static long iv[NTAB];
	float temp;

	if (*idum <= 0)
   {
		if (-(*idum) < 1) *idum=1;
		else *idum = -(*idum);
		idum2=(*idum);
		for (j=NTAB+7;j>=0;j--)
      {
			k=(*idum)/IQ1;
			*idum=IA1*(*idum-k*IQ1)-k*IR1;
			if (*idum < 0) *idum += IM1;
			if (j < NTAB) iv[j] = *idum;
		}
		iy=iv[0];
	}
	k=(*idum)/IQ1;
	*idum=IA1*(*idum-k*IQ1)-k*IR1;
	if (*idum < 0) *idum += IM1;
	k=idum2/IQ2;
	idum2=IA2*(idum2-k*IQ2)-k*IR2;
	if (idum2 < 0) idum2 += IM2;
	j=iy/NDIV;
	iy=iv[j]-idum2;
	iv[j] = *idum;
	if (iy < 1) iy += IMM1;
	if ((temp=AM*iy) > RNMX) return RNMX;
	else return temp;
}

#undef IM1
#undef IM2
#undef AM
#undef IMM1
#undef IA1
#undef IA2
#undef IQ1
#undef IQ2
#undef IR1
#undef IR2
#undef NTAB
#undef NDIV
#undef EPS
#undef RNMX
