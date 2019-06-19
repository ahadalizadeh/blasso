#include <stdio.h>
#include <stdlib.h>
#include <R.h>
#include <Rmath.h>

extern void rnorm_truncated(double*, int*, double*, double*, double*, double*);

extern "C"
{

#include "util_blasso.h"

// A Gibbs sampler for the Bayesian Lasso
void blassoGibbs(double* X, int* nn, int* pp, int* TT, int* BB, int* tthin,
   	            double* ttau, double* ssig2, double* pphi, 
			  double* sig2prior, int* fits2,
			  double* tauprior, int* fittau, int* modunc,
			  double* phiprior, int* fitphi, double* start, 
			  double* bdraws, double* sig2draws, 
			  double* taudraws, double* phidraws, double* marginc,
			  int* rrb, double* RB, double* YtY, double* YtX, 
			  int* NOISY, int* bprior, int* count)
{
   int		n = *nn; 			// no. observations
   int		p = *pp;			// no. variables
   int		T = *TT;			// no. Gibbs draws
   int		rb = *rrb;		// should we rao-blackwellize? 1(0)
   int		thin = *tthin;		// thinning the chain
   int		B = *BB;			// length of Burn in
   int		noisy = *NOISY;	// print progress
   double		tau = *ttau;		// penalty term
   double		sig2 = *ssig2;		// fixed error variance (or start value)
   double		phi = *pphi;		// prior variable incl. prob.
   double		post_shape=0.0, post_rate=0.0;	// posterior for sig2
   double		tau_rate=0.0, tau_shape=0.0;		// posterior for tau
   int		kvar=0;			// nvar in current model
   double		pscale;			// used to accomodate different
   							// priors on beta
   double		L1;				// l1 norm of betas
   int         i, j, k, tt, perc=2; // looping variables

   int		iter;
   int		nonzero;			// for var selection
   int		nvar;			// number of active predictors
   double		w;				// prob >< 0
   double		*beta;			// vector of current beta values
   double		mucon_pos, mucon_neg;  // conditional ``means'' for positive
   							// and negative components
   double		ldnpos, ldnneg, lpnpos, lpnneg; // normal dist stuff
   double		u;				// uniform RV
   int		tmpcnt;			// how many proposals for sig2 sampling
   int		badsamp;			// did sig2 sampler work?
   int		range;			// how far out should we go for sig2
   							// sampler (ie, how many cut points)?

   // for sampling from truncated normal
   int		one=1;
   double		sd, zero = 0.0;
   double		PosInf = R_PosInf, NegInf = R_NegInf;	// constants

   if(*fits2)
   {
	 post_shape = sig2prior[0] + 0.5*(double)n;

	 // add the extra bit for the different prior
	 if((*bprior) == 1) post_shape += 0.5*(double)p;
   }

   ////////////////////////////////////////////////////////////////
   ////////////////////////////////////////////////////////////////
   // Find the XtX matrix

   double* XtX = new double[p*p];

   for(i=0; i<p; i++)
   {
      for(j=0; j<p; j++)
      {
	    XtX[i*p + j] = 0.0;
	    for(k=0; k<n; k++) XtX[i*p + j] += X[i*n + k]*X[j*n + k];
      }

      RB[i] = 0.0; // the RB means...
   }

   ////////////////////////////////////////////////////////////////
   ////////////////////////////////////////////////////////////////

   // Set the starting values
   beta = new double[p];
   for(j=0; j<p; j++)
   {
      beta[j] = start[j];
	 if(*modunc && (beta[j]!=0)) kvar++;
   }

   ////////////////////////////////////////////////////////////////
   ////////////////////////////////////////////////////////////////
   //			START THE SAMPLER!!!			 //
   ////////////////////////////////////////////////////////////////
   ////////////////////////////////////////////////////////////////

   int progout = (thin*(T+B)*p) >= 40000;

   if(noisy && progout) // percent done
   {
      Rprintf("\nSampler Progress...\n| ");
	 for(i=0; i<100; i+=20) Rprintf("%d%%      ",i);
	 Rprintf("%d%% |\n|",i);
   }

   for(iter=0; iter < (T+B); iter++)
   {
      // Possibly thin the chain
      for(tt=0; tt<thin; tt++)
      {
         // allow for various priors on beta
	    if((*bprior) == 0) pscale = 1.0;
	    else pscale = sqrt(sig2);

	    L1 = 0.0;
	    nvar = 0; // reset at each iter

         // Loop over the full conditionals for beta
         for(j=0; j < p; j++)
         {
            // Conditional sd
		  sd = sqrt(sig2/XtX[j*p + j]);

            // Compute the conditional ``means''
		  mucon_pos = YtX[j]/XtX[j*p + j] - tau*sd*sd/pscale;
		  mucon_neg = YtX[j]/XtX[j*p + j] + tau*sd*sd/pscale;

            // Add the remaining piece
	       for(i=0; i<p; i++)
	       {
               if(i!=j)
	          {
	             mucon_pos -= beta[i]*XtX[i*p + j]/XtX[j*p + j];
	             mucon_neg -= beta[i]*XtX[i*p + j]/XtX[j*p + j];
	          }
	       }

		  ldnpos = dnorm(0.0, mucon_pos, sd, 1);
		  ldnneg = dnorm(0.0, mucon_neg, sd, 1);
		  lpnpos = pnorm(mucon_pos/sd, 0.0, 1.0, 1, 1);
		  lpnneg = pnorm(-mucon_neg/sd, 0.0, 1.0, 1,1);

            // nonzero may get set to 0 below if we're doing
		  // variable selection
		  nonzero = 1;
		  if(*modunc)
		  {
		     // weight in favor of = 0

               // If using the scaled prior (*bprior==1), 
			// the term is tau/(2*sigma). Otherwise the
			// term is just tau/2

			if((*bprior) == 1) w = 1.0/(1.0 + exp(log(phi) - log(1.0-phi) + log(tau) - log(2.0) - 0.5*log(sig2) + log(exp(lpnpos - ldnpos) + exp(lpnneg - ldnneg))));
			else w = 1.0/(1.0 + exp(log(phi) - log(1.0-phi) + log(tau) - log(2.0) + log(exp(lpnpos - ldnpos) + exp(lpnneg - ldnneg))));

               // Marginal inclusion probability
			//if((iter>=B) && (tt==(thin-1))) marginc[(iter-B)*p+j] = 1.0-w;
			if((iter>=B) && (tt==(thin-1))) marginc[j] += (1.0 - w) / (double)T;

			// CONTROL MODEL SIZE HERE!
			// If there are currently n variables in the model
			// and variable j is NOT already in the model,
			// force it to be zero
			//if((kvar==n) && (beta[j]==0.0)) w = 1.0;

			GetRNGstate();
			u = runif(0.0,1.0);
			PutRNGstate();

			if(u < w)
			{
			   // If beta[j] was currently nonzero, decrease
			   // the number of variables in the model
			   if(beta[j]!=0.0) kvar--;

                  beta[j] = 0.0;

			   nonzero = 0;
			}
			else
			{
			   // If beta[j] currently zero, increase
			   // the number of variables in the model
			   if(beta[j]==0.0) kvar++;
			}
		  }

            // If beta is nonzero for this draw
		  if(nonzero)
		  {
		     // Increment the number of active variables.
			// If fitting the full model, this will be
			// p at the end of the scan through beta
               nvar++; 

               // Compute the weight in favor of the nonnegative component
	          w = 1.0/(1.0 + exp(ldnpos - ldnneg + lpnneg - lpnpos));

               GetRNGstate();
	          u = runif(0.0,1.0);
	          PutRNGstate();

	          if(u < w)
		        rnorm_truncated(beta+j,&one,&mucon_pos,&sd,&zero,&PosInf);
	          else
		        rnorm_truncated(beta+j,&one,&mucon_neg,&sd,&NegInf,&zero);
		  }

	       // Update L1 norm
		  L1 += fabs(beta[j]);

            if((iter>=B) && (tt==(thin-1))) 
	       {
	          bdraws[(iter-B)*p + j] = beta[j];

	          // Compute rao-blackwellized mean if requested
	          if(rb==1) RB[j] += (w*(mucon_pos + sd*exp(dnorm(0.0,mucon_pos, sd, 1)-pnorm(mucon_pos/sd,0.0,1.0,1,1))) + (1-w)*(mucon_neg - sd*exp(dnorm(0.0, mucon_neg, sd, 1)-pnorm(-mucon_neg/sd, 0.0, 1.0, 1, 1))))/(double)T;
	       }
         }

	    // Update phi if so desired...
	    if(*modunc && *fitphi)
	    {
		  GetRNGstate();
		  phi = rbeta(phiprior[0] + kvar, phiprior[1] + p - kvar);
		  PutRNGstate();

            if((iter>=B) && (tt==(thin-1))) phidraws[iter-B] = phi;
	    }

	    // Sample sig2 if so desired...
	    if(*fits2)
	    {
            // This can all be made more efficient if we keep
		  // track of which variables are nonzero when doing
		  // variable selection...

            // compute the ``shape'' parameter
	       post_shape = sig2prior[0] + 0.5*(double)n;
            if((*bprior) == 1) post_shape += 0.5*(double)nvar;

            // compute the ``rate'' parameter
		  double tmp = *YtY;
		  if(nvar > 0)
		  {
		     for(j=0; j<p; j++) 
		     {
                  tmp += -2.0*YtX[j]*beta[j];
			   for(k=0; k<p; k++) tmp += beta[j]*XtX[j*p + k]*beta[k];
	          }
		  }
		  post_rate = sig2prior[1] + 0.5*tmp;

		  // If we are using the classic prior, just sample it...
		  if((*bprior)==0)
		  {
		     GetRNGstate();
		     sig2 = 1.0/rgamma(post_shape, 1.0/post_rate);
			PutRNGstate();
	       }
		  else // if using the rescaled prior
		  {
		     // rejection sampling here!
			// start with 3 sd's, increment if necessary
			range = 3;
			do
			{
			   badsamp = sig2_rej_samp(&sig2, post_shape, post_rate, 
			   					  tau, L1, range++, &tmpcnt);
			} while(badsamp);

			//if((iter>=B) && (tt==(thin-1))) count[iter-B] = tmpcnt;
		  }

		  if((iter>=B) && (tt==(thin-1))) sig2draws[iter-B] = sig2;
	    }

         // Sample tau if so desired
	    if(*fittau)
	    {
	       if((*bprior)==0) tau_rate = L1 + tauprior[1];
		  else tau_rate = L1/sqrt(sig2) + tauprior[1];

		  tau_shape = tauprior[0] + (double)nvar;

		  GetRNGstate();
		  tau = rgamma(tau_shape, 1.0/tau_rate);
		  PutRNGstate();

		  if((iter>=B) && (tt==(thin-1))) taudraws[iter-B] = tau;
	    }
      }

	 if(noisy && progout && ( (thin*(iter+1))/(double)(thin*(T+B)) >= perc/100.0))
	 {
         Rprintf("*");
	    perc += 2;
	 }
   }

   if(noisy && progout) Rprintf("|\n\n");

   delete[] beta; beta = NULL;
   delete[] XtX; XtX = NULL;

   return;
}


// A Gibbs sampler for the Bayesian Lasso -- sampling is done on a 
// transformation of the regression parameters
void	blassoGibbsTrans(int* nn, int* pp, int* TT, int* BB, int* tthin,
	                 double* ttau, double* ssig2, double* sig2prior,
				  int* fits2, double* tauprior, int* fittau, 
				  double* start, double* bdraws, double* sig2draws, 
				  double* taudraws, double* bhat, double* HH, 
				  double* Lambda, double* YtY, double*HtXtY,
				  int* NOISY, int* bprior, int* count)
{

   int		p = *pp;			// no. variables
   int		T = *TT;			// no. Gibbs draws
   int		thin = *tthin;		// thinning the chain
   int		B = *BB;			// length of Burn in
   double		tau = *ttau;		// penalty term
   							// on the original parameterization
   int		noisy = *NOISY;	// print progress
   double		**H;				// copy the HH matrix here
   double		**Z;				// matrix of 1's and -1's
   double		*Htbhat;			// H^T beta_hat
   double		sig2 = *ssig2;		// starting value for sig2
   double		post_shape=0.0, post_rate=0.0; // posterior pars for sig2
   double		tau_shape=0.0, tau_rate=0.0;	// posterior for tau
   double		pscale;			// used to accomodate different priors
   							// on beta
   double		L1=0.0;			// L1 norm of beta
   int		i, j, k, tt, perc=2;	// looping variables

   // Used in sampling
   int		iter;
   double		*eta;			// the transformed samples
   double		*bounds;			// the truncation points
   int		*ii;				// indices for sorting
   int		*tmp_index;
   double		*w;				// weights for the p+1 components
   double		wmax;			// max of the weights
   double		*etaH;			// used to compute weights
   double		*cumsum;
   double		*mu;				// means for the p+1 components
   double		*C;				// part of the weights
   int		easy;			// should we use rejection sampling?
   int		tmpcnt;			// how many proposals for sig2?
   int		badsamp;			// did sig2 sampler work?
   int		range;		     // how far out should we go for sig2
   							// sampler (ie, how many cut points)?

   if(*fits2) 
   {
      post_shape = sig2prior[0] + (double)(*nn)/2.0;
	 if(*bprior == 1) post_shape += 0.5*(double)p;
   }

   if(*fittau) tau_shape = p + tauprior[0];

   // If p==1, use the other function
   if(p==1)
   {
      Rprintf("\nThe option sampler = \"basic\" should be used when p = 1.\n\n");
      return;
   }


   // Precompute a few quantities
   C = new double[p+1];
   H = new double*[p];
   Z = new double*[p];
   eta = new double[p];
   Htbhat = new double[p];

   for(i=0; i<p; i++)
   {
      eta[i] = 0.0;

      H[i] = new double[p];
      for(j=0; j<p; j++) H[i][j] = HH[j*p+i];

      Z[i] = new double[p+1];
      for(j=0; j<(p+1); j++)
      {
         if(j<=i) Z[i][j] = -1.0;
	    else Z[i][j] = 1.0;
      }
   }

   for(j=0; j<p; j++)
   {
      eta[j] = Htbhat[j] = 0.0;
      for(i=0; i<p; i++)
      {
         // set the starting values
	    eta[j] += H[i][j]*start[i];
	    Htbhat[j] += H[i][j]*bhat[i];
      }
   }

   // allocate other memory
   bounds = new double[p+2];
   etaH = new double[p];
   tmp_index = new int[p];
   ii = new int[p];
   w = new double[p+1];
   mu = new double[p+1];
   cumsum = new double[p+1];

   ////////////////////////////////////////////////////////////////
   ////////////////////////////////////////////////////////////////
   //                   START THE SAMPLER!!!                     //
   ////////////////////////////////////////////////////////////////
   ////////////////////////////////////////////////////////////////

   int progout = (thin*(T+B)) >= 20000;

   if(noisy && progout) // percent done
   {
      Rprintf("\nSampler Progress...\n| ");
      for(i=0; i<100; i+=20) Rprintf("%d%%      ",i);
      Rprintf("%d%% |\n|",i);
   }

   for(iter=0; iter<(T+B); iter++)
   {
      // Possibly thin the chain
	 for(tt=0; tt<thin; tt++)
	 {
         if(*bprior==0) pscale = 1.0;
         else pscale = sqrt(sig2);

         for(j=0; j<p; j++)
         {
            // Compute the bounds
            for(k=0; k<p; k++)
	       {
               bounds[k+1] = 0.0;
	          etaH[k] = 0.0;
	          for(i=0; i<p; i++)
	          {
                  if(i!=j)
	             {
	                bounds[k+1] -= eta[i]*H[k][i];
			      etaH[k] += eta[i]*H[k][i];
		        }
	          }
	          bounds[k+1] = bounds[k+1]/H[k][j];
	          tmp_index[k] = ii[k] = k;
	       }
	       // Set the extremes
	       bounds[0] = R_NegInf; bounds[p+1] = R_PosInf;

            // Sort the bounds, get their ordering
	       // Don't include -Inf and Inf in the sorting (they are already sorted)
	       // ...want to get the index right for use later
	       rsort_with_index(bounds+1, tmp_index, p);

            // Now get the inverse indexing...this will be ii
            R_qsort_int_I(tmp_index, ii, 1, p);

            // Now loop over each of the p+1 possible configurations,
	       // compute the weights and store the means
	       easy = 0; wmax = R_NegInf;
	       for(k=0; k<(p+1); k++)
	       {
	          mu[k] = Htbhat[j];
	          w[k] = 0.0;
               for(i=0; i<p; i++)
  	          {
	             mu[k] -= tau*(sig2/pscale)*Lambda[j]*H[i][j]*Z[ii[i]][k]*sign(H[i][j]);
			   // Note: etaH is the correct sum -- it has removed i==j above
	             w[k] -= (tau/pscale)*Z[ii[i]][k]*sign(H[i][j])*etaH[i];
	          }

               if((bounds[k] <= mu[k]) && (mu[k] <= bounds[k+1])) easy = 1;

               C[k] = pnorm((bounds[k+1] - mu[k])/sqrt(sig2*Lambda[j]), 0.0, 1.0, 1, 0) - pnorm((bounds[k] - mu[k])/sqrt(sig2*Lambda[j]), 0.0, 1.0, 1, 0);

	          // compute the log unnormalized weight (still need to add C below)
	          w[k] += 0.5*mu[k]*mu[k]/(sig2*Lambda[j]);

		     if(w[k] > wmax) wmax = w[k];
	       }

            int PRINT=0;
	       if(easy) eta[j] = easy_samp(p, bounds, w, C, mu, sig2*Lambda[j]);
	       else eta[j] = rej_samp(p, bounds, w, mu, sig2*Lambda[j], wmax, PRINT);
         }

         // Compute the betas, L1 norm and save...
	    if( ((iter >= B) && (tt==(thin-1))) || ((*bprior == 1) || (*fittau == 1)))
	    {
	       L1 = 0.0;
            for(j=0; j<p; j++)
	       {
	          double tmp = 0.0;
	          for(k=0; k<p; k++) tmp += H[j][k]*eta[k];
	          L1 += fabs(tmp);
	          if((iter>=B) && (tt==(thin-1))) bdraws[(iter-B)*p + j] = tmp;
	       }
	    }

	    // Sample sig2 if so desired...
	    if(*fits2)
	    {
            post_rate = *YtY;
	       for(j=0; j<p; j++) post_rate += -2.0*eta[j]*HtXtY[j] + eta[j]*eta[j]/Lambda[j];
	       post_rate = sig2prior[1] + 0.5*post_rate;

	       // Sample...
            if(*bprior == 0)
	       {
	          GetRNGstate();
	          sig2 = 1.0 / rgamma(post_shape, 1.0/post_rate);
		     PutRNGstate();
	       }
	       else
	       {
	          // rejection sampling here!
		     // start with 3 sd's, increment if necessary
		     range = 3;
		     do
		     {
                  badsamp = sig2_rej_samp(&sig2, post_shape, post_rate,
							       tau, L1, range++, &tmpcnt);
		     } while(badsamp);

		     //if(iter>=B) count[iter-B] = tmpcnt;
	       }

	       if((iter>=B) && (tt==(thin-1))) sig2draws[(iter-B)] = sig2;
	    }

	    if(*fittau)
	    {
            if(*bprior==0) tau_rate = L1 + tauprior[1];
	       else tau_rate = L1/sqrt(sig2) + tauprior[1];

	       GetRNGstate();
	       tau = rgamma(tau_shape, 1.0/tau_rate);
	       PutRNGstate();

	       if((iter>=B) && (tt==(thin-1))) taudraws[iter-B] = tau;
	    }
	 } // end thinning

      if(noisy && progout && ( (thin*(iter+1))/(double)(thin*(T+B)) >= perc/100.0))
      {
         Rprintf("*");
         perc += 2;
      }
   }

   if(noisy && progout) Rprintf("|\n\n");

   // clean up memory!
   for(j=0; j<p; j++)
   {
      delete[] H[j]; H[j] = NULL;
      delete[] Z[j]; Z[j] = NULL;
   }
   delete[] H; H = NULL;
   delete[] Z; Z = NULL;
   delete[] eta; eta = NULL;
   delete[] Htbhat; Htbhat = NULL;
   delete[] bounds; bounds = NULL;
   delete[] tmp_index; tmp_index = NULL;
   delete[] ii; ii = NULL;
   delete[] mu; mu = NULL;
   delete[] w; w = NULL;
   delete[] cumsum; cumsum = NULL;
   delete[] etaH; etaH = NULL;
   delete[] C; C = NULL;

}


} // end extern "C"

