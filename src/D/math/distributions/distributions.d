/*
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
  // x == 0 should be valid
  assert(x > 0 &&  shape > 0 && scale > 0, "Arguments must be >0");
}
body {
  double z = x / scale;
  double result = -z + (shape-1.0) * log(z) - logGamma(shape) - log(scale);

  if(give_log) return(result);
  else return(exp(result));
}

unittest {
}


// tail probabilities of gamma distribution
double pgamma(in double x, in double shape, in double scale, in bool lower_tail = true)
in {
  assert(x > 0 && shape > 0 && scale > 0, "Arguments must be >0");
}
body {
  if(lower_tail) return(gammaIncomplete(shape, x/scale));
  else return(gammaIncompleteCompl(shape, x/scale));
}

unittest {
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
}

// poisson probabilities
double dpois(in int x, in double lambda, in bool give_log = false)
in {
  // lambda = 0 should be okay
  assert(x >= 0, "x must be >= 0");
  assert(lambda > 0, "lambda must be > 0");
}
body {
  return dgamma(lambda, cast(double)(x+1), 1.0, give_log);
}

unittest {
}


// tail probabilities of poisson
double ppois(in int x, in double lambda, in bool lower_tail = true)
in {
  assert(x >= 0, "x must be >= 0 [" ~ to!string(x) ~ "]");
  assert(lambda > 0, "lambda must be > 0");
}
body {
  return pgamma(lambda, cast(double)(x+1), 1.0, !lower_tail);
}

unittest {
}

// quantiles of poisson
int qpois(in double p, in double lambda, in bool lower_tail = true)
in {
  assert(p >= 0 && p <= 1, "p must be in [0,1]");
  assert(lambda > 0, "lambda must be >0");
}
body {
  // initial part follows R ver 2.15.0 (src/nmath/qpois.c)
  double mu = lambda;
  double sigma = sqrt(lambda);
  double gamma = 1.0/sigma;

  double pp; // can't change p because of the "in"
  if(lower_tail) pp = p;
  else pp = 1.0 - p;

  if(pp == 0.0) return(0);
  if(pp == 1.0) return(int.max);

  // p too close to 1
  if(pp + 1.01*double.epsilon >= 1.0) return int.max;

  // initial guess
  double z = qnorm(pp, 0.0, 1.0, true);
  int q = cast(int)floor(mu + sigma*(z + gamma * (z*z-1.0)/6.0) + 0.5);
  if(q < 0) q = 0;

  double current_prob = ppois(q, lambda, true);

  // slight jitter to ensure left continuity
  pp *= 1.0 - 64*double.epsilon;

  // Want the smallest integer with ppois >= pp
  // This is a very crude method:
  //     step by 1 from the initial guess
  if(pp <= current_prob) { // may need to move down
    while(pp <= current_prob) {
      current_prob -= dpois(q, lambda, false);
      q--;
    }
    return(q+1);
  }
  else { // need to move up
    while(pp > current_prob) {
      q++;
      current_prob += dpois(q, lambda, false);
    }
    return(q);
  }
}


unittest {
}
