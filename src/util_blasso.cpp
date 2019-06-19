#include <R.h>
#include <Rmath.h>
#include <stdlib.h>
#include <stdio.h>

extern void rnorm_truncated(double*, int*, double*, double*, double*, double*);

extern "C"
{

#include "util_blasso.h"

double	easy_samp(int p,
			     double* h, 
				double* w,
				double* C,
				double* mu,
				double  s2)
{
   int	k;
   double wmax = R_NegInf;
   double	*wstar, *cumsum;
   double u, ans, tmpsd = sqrt(s2);
   int	one=1;

   wstar = new double[p+1];
   cumsum = new double[p+1];

   // First add the last piece of the constant
   for(k=0; k<(p+1); k++)
   {
      wstar[k] = w[k] + log(C[k]);
	 if(wmax < wstar[k]) wmax = wstar[k];
   }

   // Rescale based on max weight and find cumsum
   cumsum[0] = wstar[0] = exp(wstar[0] - wmax);
   for(k=1; k<(p+1); k++)
   {
      wstar[k] = exp(wstar[k] - wmax);
      cumsum[k] = cumsum[k-1] + wstar[k];
   }

   // Sample the component
   GetRNGstate();
   u = unif_rand();
   PutRNGstate();

   for(k=0; k<(p+1); k++) if(u <= cumsum[k]/cumsum[p]) break;

   rnorm_truncated(&ans, &one, mu+k, &tmpsd, h+k, h+k+1);

   // clean up
   delete[] wstar; wstar = NULL;
   delete[] cumsum; cumsum = NULL;

   return ans;
}

double	rej_samp(int p,
			    double* h,
			    double* w,
			    double* mu,
			    double  s2,
			    double  wmax,
			    int print)
{
   int 		l, L;
   double 	*ci, *mi, *mustar, *wstar, *hstar;
   double		*logomega, lomax=R_NegInf, maxtmp, *cumsum, u;
   double		prop, cl;	 // cl is to avoid numerical trouble

   // Allocate memory
   ci = new double[2*p];
   mi = new double[2*p];
   mustar = new double[2*p];
   wstar = new double[2*p];
   hstar = new double[2*p+1];
   logomega = new double[2*p];
   cumsum = new double[2*p];

   // First put the weights on a difference scale
   for(l=0; l<(p+1); l++) w[l] -= wmax;

   // First get the slopes and intercepts for the pieces
   GetSlopeInt(p, ci, mi, h, mu, s2, w, print);

   // Next extend the weights, means and knots, and compute weights
   mustar[0] = mu[0]; wstar[0] = w[0]; hstar[0] = h[0];

   // notice i use h[1] here rather than hstar[1] because the latter
   // hasn't been set yet (they are the same)
   logomega[0] = ci[0] - log(mi[0]) + mi[0]*h[1];
   if(lomax < logomega[0]) lomax = logomega[0];

   for(l=1; l<p; l++)
   {
      mustar[2*l-1] = mustar[2*l] = mu[l];
	 wstar[2*l-1] = wstar[2*l] = w[l];
	 hstar[2*l-1] = h[l];

	 // Fill in the knots where the tangents intersect
	 // (it's just the midpoint...)
	 hstar[2*l] = 0.5*(h[l] + h[l+1]);

	 // Compute the two weights
	 cl = fmax2(mi[2*l-1]*hstar[2*l], mi[2*l-1]*hstar[2*l-1]);
	 logomega[2*l-1] = ci[2*l-1] + cl + log(exp(mi[2*l-1]*hstar[2*l] - cl)/mi[2*l-1] - exp(mi[2*l-1]*hstar[2*l-1] - cl)/mi[2*l-1]);

      // note i use h[l+1] here rather than hstar[2*l+1] because the latter
	 // hasn't been assigned yet (they end up being equivalent)
	 cl = fmax2(mi[2*l]*h[l+1], mi[2*l]*hstar[2*l]);
	 logomega[2*l] = ci[2*l] + cl + log(exp(mi[2*l]*h[l+1] - cl)/mi[2*l] - exp(mi[2*l]*hstar[2*l] - cl)/mi[2*l]);

      maxtmp = fmax(logomega[2*l-1], logomega[2*l]);
	 if(lomax < maxtmp) lomax = maxtmp;
   }

   mustar[2*p-1] = mu[p];
   wstar[2*p-1] = w[p];
   hstar[2*p-1] = h[p];
   logomega[2*p-1] = ci[2*p-1] - log(-mi[2*p-1]) + mi[2*p-1]*hstar[2*p-1];

   hstar[2*p] = h[p+1];

   if(print) for(l=0;l<(2*p);l++)Rprintf("logomega[%d] = %f\n",l,logomega[l]);
   if(print) {for(l=0; l<=(2*p); l++) Rprintf("hstar[%d] = %.10f\n",l,hstar[l]);}

   if(lomax < logomega[2*p-1]) lomax = logomega[2*p-1];

   // normalize the weights
   cumsum[0] = exp(logomega[0] - lomax);
   for(l=1; l < (2*p); l++) cumsum[l] = cumsum[l-1] + exp(logomega[l] - lomax);

   if(print) for(l=0; l < (2*p); l++) Rprintf("cumsum[%d] = %f\n", l,cumsum[l]/cumsum[2*p-1]);

   // Now do the rejection sampling
   GetRNGstate();
   do
   {
      // Sample a component
	 u = unif_rand();
	 for(L=0; L<(2*p); L++) if(u <= cumsum[L]/cumsum[2*p-1]) break;

	 if(print) Rprintf("\nsampled %d\n", L);
	 
	 // Sample from the truncated exponential
	 // Need to be careful about the left tail
	 if(L!=0) prop = hstar[L] + log( (exp(mi[L]*(hstar[L+1]-hstar[L])) - 1.0)*unif_rand() + 1.0)/mi[L];
	 else prop = hstar[1] + log(unif_rand())/mi[0]; // L==0 case

   } while(log(unif_rand()) > (logtarget(prop, mustar[L], s2, wstar[L]) - ci[L] - prop*mi[L])); // check to see if we have to reject

   PutRNGstate();

   // clean up
   delete[] ci; ci = NULL;
   delete[] mi; mi = NULL;
   delete[] mustar; mustar = NULL;
   delete[] wstar; wstar = NULL;
   delete[] hstar; hstar = NULL;
   delete[] logomega; logomega= NULL;
   delete[] cumsum ; cumsum = NULL;

   return(prop);
}

void		GetSlopeInt(int	p,
				  double* ci,	// intercepts
				  double* mi,	// slopes
				  double* h,	// bounds
				  double* mu,	// component means
				  double  s2,	// component variance
				  double* w,
				  int	print)	// weights
{
   int	l;

   // the first one is for the left tail
   mi[0] = dlogtarget(h[1], mu[0], s2);
   ci[0] = w[0] - 0.5*(mu[0]*mu[0] - h[1]*h[1])/s2;

   // Now do the interior
   for(l=1; l<p; l++)
   {
      mi[2*l - 1] = dlogtarget(h[l], mu[l], s2);
	 mi[2*l] = dlogtarget(h[l+1], mu[l], s2);

	 ci[2*l - 1] = w[l] - 0.5*(mu[l]*mu[l] - h[l]*h[l])/s2;
	 ci[2*l] = w[l] - 0.5*(mu[l]*mu[l] - h[l+1]*h[l+1])/s2;
   }

   // the final one is for the right tail
   mi[2*p-1] = dlogtarget(h[p], mu[p], s2);

   ci[2*p-1] = w[p] - 0.5*(mu[p]*mu[p] - h[p]*h[p])/s2;
}

double	dlogtarget(double x, double m, double s2)
{
   return(-(x-m)/s2);
}

double 	logtarget(double x, double m, double s2, double W)
{
   return(W - 0.5*(x-m)*(x-m)/s2);
}

// This function performs rejection sampling to sample v = sigma^{-1}
// We transform before saving v^{-2} = sigma^2
// The function returns 0 is it is taking too long
int		sig2_rej_samp(double	*ans, // sig2 sample
				    double	a,	// "shape"
				    double	b,	// "rate"
				    double	tau,	// penalty
				    double	L1,  // l1 norm of betas
				    int		outer, // how far out should we look?
				    				  // "outer" sd's on each side
				    int		*count) // how many proposals?
{
   int			k,l;		// loop
   double			*knots, *iknots; // knots and interior knots
   double			*mi, *ci; // slopes and intercepts
   double			*logw;	// log weights for the chunks
   double			*cumsum;	// cumulative weights
   double			cl;		// used for numerical purposes
   double			maxw;	// keep track of the largest weight
   double			prop=0.0;	// proposal
   int			K;		// number of knots
   double			vhat;	// mode of distribution
   double			sd;		// normal approx at mode -- used
   						// to figure out spacing for knots
   double			u;		// for sampling
   int			L;		// ibid.
   int			cnt=-1;	// how many proposals must we take?
   int			bad=0;	// the value we return

   // Find the mode
   vhat = (-tau*L1 + sqrt(tau*tau*L1*L1 + 8.0*b*(2.0*a - 1.0)))/(4.0*b);

   // Find the second derivative at the mode
   sd = 1.0/sqrt(fabs((1.0-2.0*a)/(vhat*vhat) - 2.0*b));

   // make sure we go out at least 2 sds
   outer = (outer >= 2) ? outer : 2; 

   ///////////////////////////////////////
   // Set the knots to use for sampling //
   ///////////////////////////////////////
   
   // Find the smallest possible knot location
   for(k=outer; k>=1; k--) if( (vhat - (double)k*sd) > 0.0 ) break;

   // if vhat - k*sd were negative for all of them, check 0.5
   if(k == 0)
   {
      K = outer + 2;
	 knots = new double[K];

      // if 0.5 works, then we need -0.5, +0.5, +1, ..., +outer
      if( (vhat - 0.5*sd) > 0.0 ) knots[0] = vhat - 0.5*sd;

      // just use midpoint between zero and vhat
	 else knots[0] = 0.5*vhat;
   }
   else // other scenario is that we stopped earlier
   {
      K = outer + 2 + k;
	 knots = new double[K];
	 for(l=0; l<k; l++) knots[l] = vhat - (double)(k - l)*sd;
	 knots[l] = vhat - 0.5*sd;
   }
   // Now fill in the positive side
   knots[K-outer-1] = vhat + 0.5*sd;
   for(l=1; l<=outer; l++) knots[K-outer+l-1] = vhat + (double)l*sd;

   /////////////////////////////////
   // Done setting knot locations //
   /////////////////////////////////

   // Now compute the the slopes and intercepts through these knots
   // Also find the interior knots
   // AND find the weights
   ci = new double[K]; mi = new double[K]; iknots = new double[K-1];
   logw = new double[K]; cumsum = new double[K];
   maxw = R_NegInf;

   for(k=0; k<K; k++)
   {
      ci[k] = (2.0*a - 1.0)*log(knots[k]) - 2.0*a + 1.0 + b*knots[k]*knots[k];
	 mi[k] = (2.0*a - 1.0)/knots[k] - 2.0*b*knots[k] - tau*L1;

	 if(k < (K-1))
	 {
	    iknots[k] = ((2.0*a-1.0)*(log(knots[k+1])-log(knots[k])) + b*(knots[k+1] + knots[k])*(knots[k+1] - knots[k]))/( (2.0*a-1.0)*(1.0/knots[k] - 1.0/knots[k+1]) + 2.0*b*(knots[k+1] - knots[k]));

	    // Find the weight for this chunk
	    if(k==0) logw[k] = -log(mi[0]) + ci[0] + mi[0]*iknots[0] + log(1.0 - exp(-mi[0]*iknots[0]));
	    else 
	    {
            cl = fmax2(mi[k]*iknots[k], mi[k]*iknots[k-1]);
            logw[k] = ci[k] + cl + log(exp(mi[k]*iknots[k] - cl)/mi[k] - exp(mi[k]*iknots[k-1] - cl)/mi[k]);
	    }
	    
	    maxw = (logw[k] > maxw) ? logw[k] : maxw;
	 }
   }
   logw[K-1] = ci[K-1] + mi[K-1]*iknots[K-2]  - log(-mi[K-1]);
   maxw = (logw[K-1] > maxw) ? logw[K-1] : maxw;

   // Compute the cumulative weights
   cumsum[0] = exp(logw[0] - maxw);
   for(l=1; l<K; l++) cumsum[l] = cumsum[l-1] + exp(logw[l] - maxw);

   // Now do the rejection sampling
   GetRNGstate();
   do
   {
      if(++cnt == 500) {bad = 1; break;}
      // Sample a component
	 u = unif_rand();
	 for(L=0; L<K; L++) if(u <= cumsum[L]/cumsum[K-1]) break;

	 // Sample from the truncated exponential
	 // Need to be careful about the tails
	 if(L==0) prop = log(1.0-unif_rand()*(1.0-exp(mi[0]*iknots[0])))/mi[0];

	 else if(L<(K-1)) prop = iknots[L-1] + log(1.0 - unif_rand()*(1.0 - exp(mi[L]*(iknots[L]-iknots[L-1]))))/mi[L];

	 else prop = iknots[K-2] + log(unif_rand())/mi[K-1];

   } while(log(unif_rand()) > ((2.0*a - 1.0)*log(prop) - b*prop*prop - tau*L1*prop - ci[L] - mi[L]*prop)); // check to see if we have to reject

   PutRNGstate();

   *ans = exp(-2.0*log(prop));
   *count = ++cnt;

   // clean up memory
   delete[] knots; knots = NULL;
   delete[] iknots; iknots = NULL;
   delete[] ci; ci = NULL;
   delete[] mi; mi = NULL;
   delete[] logw; logw = NULL;
   delete[] cumsum; cumsum = NULL;

   return bad;
}


} // extern "C"


