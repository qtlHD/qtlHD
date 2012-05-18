/*
 * distributions: densities, tail probabilities, and quantiles
 */

module math.distributions.distributions;

import std.algorithm;
import std.math;
import std.mathspecial;
import std.stdio;
import std.string;
import std.conv;


// normal density
double dnorm(double x, double mu, double sigma, bool give_log = false)
in {
  assert(sigma >= 0, "sigma must be >= 0");
}
body {
  double result = -(x-mu)*(x-mu)/sigma/sigma/2.0 - log(sigma) - log(2.0*PI)/2.0;

  if(give_log) return(result);
  else return(exp(result));
}

unittest {
  double[] x = [-2.0, 0.0, 3.0, 4.0];
  double[] logdnorm1_3 = [-2.5175508218727822296, -2.0731063774283380319, -2.2397730440950045505, -2.5175508218727822296];
  double[] logdnorm5_01 = [-2448.616353440210787085, -1248.616353440210559711, -198.616353440210644976, -48.616353440210630765];
  double[] dnorm1_3 = [0.080656908173047783817, 0.125794409230997716875, 0.106482668507450750628, 0.080656908173047783817];
  double[] dnorm0_2 = [0.120985362259571682664, 0.199471140200716351432, 0.064758797832945871886, 0.026995483256594031418];

  foreach(i, xv; x) {
    double y = dnorm(xv, 1.0, 3.0, true);
    assert(abs(y - logdnorm1_3[i]) < 1e-12);

    y = dnorm(xv, 5.0, 0.1, true);
    assert(abs(y - logdnorm5_01[i]) < 1e-12);

    y = dnorm(xv, 1.0, 3.0);
    assert(abs(y - dnorm1_3[i]) < 1e-12);

    y = dnorm(xv, 0.0, 2.0);
    assert(abs(y - dnorm0_2[i]) < 1e-12);
  }
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
  double[] x = [-2.0, 0.0, 3.0, 4.0];
  double[] pnorm1_3 = [0.15865525393145704647, 0.36944134018176366663, 0.74750746245307708726, 0.84134474606854292578];
  double[] pnorm0_2 = [0.15865525393145704647, 0.50000000000000000000, 0.93319279873114191481, 0.97724986805182079141];

  foreach(i, xv; x) {
    double y = pnorm(xv, 1.0, 3.0);
    assert(abs(y - pnorm1_3[i]) < 1e-12);

    y = pnorm(xv, 1.0, 3.0, false);
    assert(abs(y - (1.0 - pnorm1_3[i])) < 1e-12);

    y = pnorm(xv, 0.0, 2.0);
    assert(abs(y - pnorm0_2[i]) < 1e-12);

    y = pnorm(xv, 0.0, 2.0, false);
    assert(abs(y - (1.0 - pnorm0_2[i])) < 1e-12);
  }
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
  double[] p = [1e-5, 0.01, 0.25, 0.999];
  double[] qnorm1_3 = [-11.7946723817684766544, -5.9790436221225222724, -1.0234692505882452274, 10.2706969185034395764];
  double[] qnorm0_2 = [-8.5297815878456511030, -4.6526957480816815149, -1.3489795003921634109, 6.1804646123356263843];
  double[] uqnorm1_3 = [13.7946723817684766544, 7.9790436221225222724, 3.0234692505882452274, -8.2706969185034395764];
  double[] uqnorm0_2 = [8.5297815878456511030, 4.6526957480816815149, 1.3489795003921634109, -6.1804646123356263843];

  foreach(i, pv; p) {
    double y = qnorm(pv, 1.0, 3.0);
    assert(abs(y - qnorm1_3[i]) < 1e-12);

    y = qnorm(pv, 1.0, 3.0, false);
    assert(abs(y - uqnorm1_3[i]) < 1e-12);

    y = qnorm(pv, 0.0, 2.0);
    assert(abs(y - qnorm0_2[i]) < 1e-12);

    y = qnorm(pv, 0.0, 2.0, false);
    assert(abs(y - uqnorm0_2[i]) < 1e-12);
  }

  // check the 0, 1 cases
  assert(qnorm(1.0, 0.0, 1.0) == real.infinity);
  assert(qnorm(0.0, 0.0, 1.0) == -real.infinity);
  assert(qnorm(0.0, 0.0, 1.0) != real.infinity);
}



// gamma density
double dgamma(in double x, in double shape, in double scale, in bool give_log = false)
in {
  assert(x > 0 && shape > 0 && scale > 0, "Arguments must be >0");
}
body {
  double z = x / scale;
  double result = -z + (shape-1.0) * log(z) - logGamma(shape) - log(scale);

  if(give_log) return(result);
  else return(exp(result));
}

unittest {
  double[] x = [0.1, 2.0, 3.5, 7.0];
  double[] dgamma2_3 = [0.010746845560911176543, 0.114092693118353794013, 0.121101253744565762194, 0.075422641672315049455];
  double[] logdgamma2_3 = [-4.5331430036635982361, -2.1707440634429406856, -2.1111282755075180262, -2.5846477616142395917];
  double[] dgamma9_2 = [4.6078124249205551814e-16, 4.5619970383363459403e-06, 1.8955643271209153978e-04, 8.4326320176074404805e-03];
  double[] logdgamma9_2 = [-35.3136082717371166950, -12.2977500833051944795, -8.5708237798218132042, -4.7756463353422518026];

  foreach(i, xv; x) {
    double y = dgamma(xv, 2.0, 3.0);
    assert(abs(y - dgamma2_3[i]) < 1e-12);

    y = dgamma(xv, 2.0, 3.0, true);
    assert(abs(y - logdgamma2_3[i]) < 1e-12);

    y = dgamma(xv, 9.0, 2.0);
    assert(abs(y - dgamma9_2[i]) < 1e-12);

    y = dgamma(xv, 9.0, 2.0, true);
    assert(abs(y - logdgamma9_2[i]) < 1e-12);
  }
}


// tail probabilities of gamma distribution
/*
double pgamma(in double x, in double shape, in double scale, in bool lower_tail)
in {
  assert(x > 0 && shape > 0 && scale > 0, "Arguments must be >0");
}
body {


}
*/

// quantiles of gamma distribution

// poisson probabilities

// tail probabilities of poisson

// quantiles of poisson

