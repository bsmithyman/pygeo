// Written / Provided by Lawrence Jones on site:
// http://www.velocityreviews.com/forums/t438469-c-code-for-converting-ibm-370-floating-point-to-ieee-754-a.html
// Modified slightly to use Py_ssize_t for array sizing.

/* ibm2ieee - Converts a number from IBM 370 single precision floating
point format to IEEE 754 single precision format. For normalized
numbers, the IBM format has greater range but less precision than the
IEEE format. Numbers within the overlapping range are converted
exactly. Numbers which are too large are converted to IEEE Infinity
with the correct sign. Numbers which are too small are converted to
IEEE denormalized numbers with a potential loss of precision (including
complete loss of precision which results in zero with the correct
sign). When precision is lost, rounding is toward zero (because it's
fast and easy -- if someone really wants round to nearest it shouldn't
be TOO difficult). */

#include "fpconvert.h"


void ibm2ieee(void *to, const void *from, Py_ssize_t len)
{
register unsigned fr; /* fraction */
register int exp; /* exponent */
register int sgn; /* sign */

for (; len-- > 0; to = (char *)to + 4, from = (char *)from + 4) {
/* split into sign, exponent, and fraction */
fr = ntohl(*(long *)from); /* pick up value */
sgn = fr >> 31; /* save sign */
fr <<= 1; /* shift sign out */
exp = fr >> 25; /* save exponent */
fr <<= 7; /* shift exponent out */

if (fr == 0) { /* short-circuit for zero */
exp = 0;
goto done;
}

/* adjust exponent from base 16 offset 64 radix point before first digit
to base 2 offset 127 radix point after first digit */
/* (exp - 64) * 4 + 127 - 1 == exp * 4 - 256 + 126 == (exp << 2) - 130 */
exp = (exp << 2) - 130;

/* (re)normalize */
while (fr < 0x80000000) { /* 3 times max for normalized input */
--exp;
fr <<= 1;
}

if (exp <= 0) { /* underflow */
if (exp < -24) /* complete underflow - return properly signed zero */
fr = 0;
else /* partial underflow - return denormalized number */
fr >>= -exp;
exp = 0;
} else if (exp >= 255) { /* overflow - return infinity */
fr = 0;
exp = 255;
} else { /* just a plain old number - remove the assumed high bit */
fr <<= 1;
}

done:
/* put the pieces back together and return it */
*(unsigned *)to = (fr >> 9) | (exp << 23) | (sgn << 31);
}
}

/* ieee2ibm - Converts a number from IEEE 754 single precision floating
point format to IBM 370 single precision format. For normalized
numbers, the IBM format has greater range but less precision than the
IEEE format. IEEE Infinity is mapped to the largest representable
IBM 370 number. When precision is lost, rounding is toward zero
(because it's fast and easy -- if someone really wants round to nearest
it shouldn't be TOO difficult). */

void ieee2ibm(void *to, const void *from, Py_ssize_t len)
{
register unsigned fr; /* fraction */
register int exp; /* exponent */
register int sgn; /* sign */

for (; len-- > 0; to = (char *)to + 4, from = (char *)from + 4) {
/* split into sign, exponent, and fraction */
fr = *(unsigned *)from; /* pick up value */
sgn = fr >> 31; /* save sign */
fr <<= 1; /* shift sign out */
exp = fr >> 24; /* save exponent */
fr <<= 8; /* shift exponent out */

if (exp == 255) { /* infinity (or NAN) - map to largest */
fr = 0xffffff00;
exp = 0x7f;
goto done;
}
else if (exp > 0) /* add assumed digit */
fr = (fr >> 1) | 0x80000000;
else if (fr == 0) /* short-circuit for zero */
goto done;

/* adjust exponent from base 2 offset 127 radix point after first digit
to base 16 offset 64 radix point before first digit */
exp += 130;
fr >>= -exp & 3;
exp = (exp + 3) >> 2;

/* (re)normalize */
while (fr < 0x10000000) { /* never executed for normalized input */
--exp;
fr <<= 4;
}

done:
/* put the pieces back together and return it */
fr = (fr >> 8) | (exp << 24) | (sgn << 31);
*(unsigned *)to = htonl(fr);
}
}

/* Test harness for IEEE systems */
#ifdef TEST
#define MAX	1000000	 /* number of iterations */
#define IBM_EPS	4.7683738e-7	/* worst case error */

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

double epsm;

void check(float f1)
{
int exp;
float f2;
double eps;
unsigned ibm1, ibm2;

frexp(f1, &exp);
ieee2ibm(&ibm1, &f1, 1);
ibm2ieee(&f2, &ibm1, 1);
ieee2ibm(&ibm2, &f2, 1);
if (memcmp(&ibm1, &ibm2, sizeof ibm1) != 0)
printf("Error: %08x <=> %08x\n", *(unsigned *)&ibm1, *(unsigned *)&ibm2);
eps = ldexp(fabs(f1 - f2), -exp);
if (eps > epsm) epsm = eps;
if (eps > IBM_EPS)
printf("Error: %.8g != %.8g\n", f1, f2);
}

int main()
{
int i;
float f1;

epsm = 0.0;
for (i = 0; i < MAX; i++) {
f1 = drand48();
check(f1);
check(-f1);
}
printf("Max eps: %g\n", epsm);
return 0;
}

#endif
