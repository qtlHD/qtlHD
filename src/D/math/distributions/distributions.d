/**
 * distributions: densities, tail probabilities, and quantiles
 * tests in test/simulate/distributions
 */

module math.distributions.distributions;

import std.stdio;
import std.math;
import std.mathspecial;
import std.algorithm;
import std.string;
import std.conv;

// normal density
double dnorm(double x, double mu, double sigma, bool give_log = false)
in {
  // sigma should be > 0 not just >= 0
  assert(sigma >= 0, "sigma must be >= 0");
}
body {
  double result = -(x-mu)*(x-mu)/sigma/sigma/2.0 - log(sigma) - log(2.0*PI)/2.0;

  if(give_log) return(result);
  else return(exp(result));
}

unittest {
  writeln("Unit test " ~ __FILE__);
  
  assert(abs(dnorm(0.0, 0.0, 1.0) -  1.0/sqrt(2.0*PI)) < 1e-12);
}

// tail probabilities for normal distribution
double pnorm(in double x, in double mu, in double sigma, in bool lower_tail = true)
in {
  assert(sigma >= 0, "sigma must be >= 0");
}
body {
  double z = (x-mu)/sigma;

  if(z > 0) {
    if(lower_tail) return(1.0 - normalDistribution(-z));
    else return(normalDistribution(-z));
  }
  else {
    if(lower_tail) return(normalDistribution(z));
    else return(1.0 - normalDistribution(z));
  }
}

unittest {
  assert(abs(pnorm(0.0, 0.0, 1.0) - 0.5) < 1e-12);
}

// quantiles of normal distribution
double qnorm(in double p, in double mu, in double sigma, in bool lower_tail = true)
in {
  assert(sigma >= 0, "sigma must be >= 0");
  assert(p >= 0 && p <= 1, "p must be in [0,1]");
}
body {
  double z = normalDistributionInverse(p);
  if(!lower_tail) z = -z;

  return(z*sigma + mu);
}

unittest {
  assert(abs(qnorm(0.5, 0.0, 1.0)) < 1e-12);

  assert(qnorm(1.0, 0.0, 1.0) == real.infinity);
  assert(qnorm(0.0, 0.0, 1.0) == -real.infinity);
  assert(qnorm(0.0, 0.0, 1.0) != real.infinity);
}



// gamma density
double dgamma(in double x, in double shape, in double scale, in bool give_log = false)
in {
  assert(x >= 0, "x should be >= 0");
  assert(shape > 0 && scale > 0, "shape and scale must be >0");
}
body {
  double result;

  if(x==0.0) { // deal with the 0 endpoint
    if(shape > 1.0) result = 0;
    else if(shape==1.0) result = 1.0;
    else result = real.infinity;

    if(give_log) return(log(result));
    else return(result);
  }

  double z = x / scale;
  result = -z + (shape-1.0) * log(z) - logGamma(shape) - log(scale);

  if(give_log) return(result);
  else return(exp(result));
}

unittest {
  assert(abs(dgamma(0.0, 1.0, 1.0) - 1.0) < 1e-12);
  assert(abs(dgamma(0.0, 2.0, 1.0)) < 1e-12);
  assert(dgamma(0.0, 0.1, 1.0) == real.infinity);

  assert(abs(dgamma(0.0, 1.0, 1.0, true)) < 1e-12);
  assert(dgamma(0.0, 2.0, 1.0, true) == -real.infinity);
  assert(dgamma(0.0, 0.1, 1.0, true) == real.infinity);
}


// tail probabilities of gamma distribution
double pgamma(in double x, in double shape, in double scale, in bool lower_tail = true)
in {
  assert(x >= 0, "x should be >= 0");
  assert(shape > 0 && scale > 0, "shape and scale must be >0");
}
body {
  if(lower_tail) return(gammaIncomplete(shape, x/scale));
  else return(gammaIncompleteCompl(shape, x/scale));
}

unittest {
  assert(abs(pgamma(0.0, 1.0, 1.0)) < 1e-12);
  assert(abs(pgamma(0.0, 0.1, 1.0)) < 1e-12);
  assert(abs(pgamma(0.0, 2.0, 1.0)) < 1e-12);
}


// quantiles of gamma distribution
double qgamma(in double p, in double shape, in double scale, in bool lower_tail = true)
in {
  assert(shape > 0 && scale > 0, "shape and scale  must be >0");
  assert(p >= 0 && p <= 1, "p must be in [0,1]");
}
body {
  if(lower_tail) return(gammaIncompleteComplInverse(shape, 1.0-p)*scale);
  else return(gammaIncompleteComplInverse(shape, p)*scale);
}

unittest {
  assert(abs(qgamma(0.0, 1.0, 1.0)) < 1e-12);
}

// poisson probabilities
double dpois(in uint x, in double lambda, in bool give_log = false)
in {
  // lambda = 0 should be okay
  assert(x >= 0, "x must be >= 0");
  assert(lambda > 0, "lambda must be > 0");
}
body {
  return dgamma(lambda, cast(double)(x+1), 1.0, give_log);
}

unittest {
  assert(abs(dpois(0, 5.0) - exp(-5.0)) < 1e-12);
  assert(abs(dpois(0, 5.0, true) - (-5.0)) < 1e-12);
}


// tail probabilities of poisson
double ppois(in uint x, in double lambda, in bool lower_tail = true)
in {
  assert(x >= 0, "x must be >= 0 [" ~ to!string(x) ~ "]");
  assert(lambda > 0, "lambda must be > 0");
}
body {
  return pgamma(lambda, cast(double)(x+1), 1.0, !lower_tail);
}

unittest {
  assert(abs(ppois(0, 1.0) - dpois(0, 1.0)) < 1e-12);
  assert(abs(ppois(1, 1.0) - 2*exp(-1.0)) < 1e-12);
}

// quantiles of poisson
uint qpois(in double p, in double lambda, in bool lower_tail = true)
in {
  assert(p >= 0 && p <= 1, "p must be in [0,1]");
  assert(lambda > 0, "lambda must be >0");
}
body {

  double pp; // can't change p because of the "in"
  if(lower_tail) pp = p;
  else pp = 1.0 - p;

  if(pp == 0.0) return(0);
  if(pp + 1.01*double.epsilon >= 1.0) return int.max; // close to 1

  // slight jitter to ensure left continuity
  pp *= 1.0 - 64*double.epsilon;

  // initial guess
  // following R ver 2.15.0 (src/nmath/qpois.c)
  double mu = lambda;
  double sigma = sqrt(lambda);
  double z = qnorm(pp, 0.0, 1.0, true);
  int guess = cast(int)(mu + sigma*z + (z*z-1.0)/6.0 + 0.5);
  if(guess < 0) guess = 0;

  double current_prob = ppois(guess, lambda, true);

  // Want the smallest integer with ppois >= pp
  // This is a very crude method: step by 1 from the initial guess
  if(pp <= current_prob) { // may need to move down
    while(pp <= current_prob) {
      current_prob -= dpois(guess, lambda, false);
      guess--;
    }
    return(guess+1);
  }
  else { // need to move up
    while(pp > current_prob) {
      guess++;
      current_prob += dpois(guess, lambda, false);
    }
    return(guess);
  }
}


unittest {
  assert( qpois( ppois(4, 3.0)+dpois(5, 3.0)/2.0, 3.0) == 5);
  assert( qpois( ppois(0, 3.0)+dpois(1, 3.0)/2.0, 3.0) == 1);
  assert( qpois( dpois(0, 3.0)/2.0, 3.0) == 0);
}
