
#include <stdint.h>
#include <math.h>
#include "random_numrec.h"

uint64_t u, v, w;
double norm_num;
int norm_saved;

/***************************************************
	int64_t
***************************************************/
inline uint64_t int64()
{
   uint64_t x;
   u = u * 2862933555777941757LL + 7046029254386353087LL;
   v ^= v >> 17; v ^= v << 31; v ^= v >> 8;
   w = 4294957665U*(w & 0xffffffff) + (w >> 32);
   x = u ^ (u << 21); x ^= x >> 35; x ^= x << 4;
   return (x + v) ^ w;
}

/***************************************************
	setSeed
***************************************************/
inline void setseed(uint64_t seedmodifiervalue)
{
   v = 4101842887655102017LL;
   w = 1;
   u = seedmodifiervalue ^ v; int64();
   v = u; int64();
   w = v; int64();
   norm_saved = 0;
}

/***************************************************
	uniform distribution
***************************************************/\
inline double Rand()
{
   return 5.42101086242752217E-20 * int64();
}

/***************************************************
	gaussian distribution
***************************************************/
inline double gasdev()
{
   double x1, x2, w;
   if(!norm_saved)
   {
      do
      {
         x1 = 2.0 * Rand() - 1.0;
         x2 = 2.0 * Rand() - 1.0;
         w = x1 * x1 + x2 * x2;
      } while(w >= 1.0 || w == 0);
      w = sqrt((-2.0 * log(w)) / w);
      norm_num = x1 * w;
      norm_saved = 1;
      return (x2 * w);
   }
   else
   {
      norm_saved = 0;
      return norm_num;
   }
}

