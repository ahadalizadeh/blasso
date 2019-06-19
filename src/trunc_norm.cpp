#include <R.h>
#include <Rmath.h>
#include <assert.h>

// ======================================================================
// norm_rs(a, b)
// generates a sample from a N(0,1) RV restricted to be in the interval
// (a,b) via rejection sampling.
// ======================================================================
double
norm_rs(double a, double b)
{
   double	x;
   GetRNGstate();
   x = norm_rand();
   while( (x < a) || (x > b) ){
      x = norm_rand();
   }
   PutRNGstate();
   return x;
}

// ======================================================================
// half_norm_rs(a, b)
// generates a sample from a N(0,1) RV restricted to the interval
// (a,b) (with a > 0) using half normal rejection sampling.
// ======================================================================
double
half_norm_rs(double a, double b)
{
   double 	x;

   assert(a >= 0); // check it

   GetRNGstate();
   x = fabs(norm_rand());
   while( (x<a) || (x>b) ) x = fabs(norm_rand());
   PutRNGstate();
   return x;
}

// ======================================================================
// unif_rs(a, b)
// generates a sample from a N(0,1) RV restricted to the interval
// (a,b) using uniform rejection sampling. 
// ======================================================================
double
unif_rs(double a, double b)
{
   double xstar, logphixstar, x, logu;

   // Find the argmax (b is always >= 0)
   // This works because we want to sample from N(0,1)
   GetRNGstate();
   if(a <= 0.0) xstar = 0.0;
   else xstar = a;
   logphixstar = dnorm(xstar, 0.0, 1.0, 1);

   x = unif_rand()*(b-a) + a;
   logu = log(unif_rand());
   while( logu > (dnorm(x, 0.0, 1.0, 1) - logphixstar))
   {
      x = unif_rand()*(b-a) + a;
      logu = log(unif_rand());
   }
   PutRNGstate();
   return x;
}

// ======================================================================
// exp_rs(a, b)
// generates a sample from a N(0,1) RV restricted to the interval
// (a,b) using exponential rejection sampling.
// This function should be called by rnorm_truncated (where Get/PutRNGstate
// are invoked)
// ======================================================================
double
exp_rs(double a, double b)
{
   double	z;

   // Generate a proposal on (0, b-a)
   GetRNGstate();
   z = -log(1 - unif_rand()*(1 - exp(-(b-a)))) / a;
   while( log(unif_rand()) > (-0.5*z*z))
   {
      z = -log(1 - unif_rand()*(1 - exp(-(b-a)))) / a;
   }
   PutRNGstate();
   return(z+a);
}

// ======================================================================
// rnorm_trunc(sample, n, mu, sigma, lower, upper)
//
// generates 'n' random normal RVs with mean 'mu' and standard
// deviation 'sigma', truncated to the interval (lower,upper), where
// lower can be -Inf and upper can be Inf.
// Stores the result in 'sample'
// mu, sigma, lower and upper must be scalars -- same values for each 
// sample
// ======================================================================

void rnorm_truncated (double *sample,  int *n, double *mu, 
				  double *sigma, double *lower, double *upper)
{
  int		k;
  int		change = 0;
  double	a, b;
  double	logt1 = log(0.150), logt2 = log(2.18), t3 = 0.725, t4 = 0.45;
  double	z, tmp, lograt;

  a = (*lower - *mu)/(*sigma);
  b = (*upper - *mu)/(*sigma);

  if(a==b) Rprintf("Warning!! a=%f, b=%f\n",a,b);
  if(a>b) Rprintf("Warning!! a = %f > b = %f\n",a,b);
  
  for (k=0; k<(*n); k++)
  {
     change=0;
     // First scenario
     if( (a == R_NegInf) || (b == R_PosInf))
     {
        if(a == R_NegInf)
	{
           change = 1;
	   a = -b;
	   b = R_PosInf;
	}

	// The two possibilities for this scenario
     if(a <= t4) z = norm_rs(a, b);
	else z = exp_rs(a, b);
	if(change) z = -z;
     }
     // Second scenario
     else if((a * b) <= 0.0)
     {
        // The two possibilities for this scenario
        if((dnorm(a, 0.0, 1.0, 1) <= logt1) || (dnorm(b, 0.0, 1.0, 1) <= logt1))
	   {
           z = norm_rs(a, b);
	   }
	   else z = unif_rs(a,b);
     }
     // Third scenario
     else
     {
        if(b < 0)
	{
	   tmp = b; b = -a; a = -tmp; change = 1;
	}

	lograt = dnorm(a, 0.0, 1.0, 1) - dnorm(b, 0.0, 1.0, 1);
	if(lograt <= logt2) z = unif_rs(a,b);
	else if((lograt > logt1) && (a < t3)) z = half_norm_rs(a,b);
	else z = exp_rs(a,b);
	if(change) z = -z;
     }

     sample[k] = *sigma*z + *mu;
  }
}

