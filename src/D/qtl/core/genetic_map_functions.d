/**
 * Genetic map function module
 **/

module qtl.core.map_functions;

import std.stdio;
import std.math;
import std.numeric; // contains findRoot()
import std.conv;

enum GeneticMapFunc { Haldane, Kosambi, Morgan, Carter_Falconer };

// convert cM distance to recombination fraction
// Note: one chromosome at a time
double[] dist_to_recfrac(double[] dist_cM, GeneticMapFunc which_dist_to_recfrac = GeneticMapFunc.Haldane)
in {
  foreach(dist; dist_cM) {
    assert(dist >= 0, "dist_cM must be >= 0, was " ~ to!string(dist));
  }
}
body {
  auto rec_frac = new double[dist_cM.length];
  enum TOL = 1e-14;

  double mfcfsub_dist;
  double mfcfsub(double rf)
  {
    return recfrac_to_dist([rf], GeneticMapFunc.Carter_Falconer)[0]-mfcfsub_dist;
  }

  foreach(i, dist; dist_cM) {
    switch(which_dist_to_recfrac) {
    case GeneticMapFunc.Haldane:
      rec_frac[i] = 0.5*(1-exp(-dist/50));
      break;
    case GeneticMapFunc.Kosambi:
      rec_frac[i] = 0.5*tanh(dist/50);
      break;
    case GeneticMapFunc.Morgan:
      if(dist >= 50) {
        rec_frac[i] = 0.5;
      }
      else {
        rec_frac[i] = dist/100.0;
      }
      break;
    case GeneticMapFunc.Carter_Falconer:
      if(dist >= 300) { dist = 300; }
      if(dist == 0) {
        rec_frac[i] = 0; 
      }
      else {
        mfcfsub_dist = dist;
        rec_frac[i] = findRoot(&mfcfsub, 0.0, 0.5-TOL);
      }
      break;
    default: break;
    }
  }
  return rec_frac;
}

unittest {
  writeln("Unit test " ~ __FILE__, " : dist_to_recfrac");

  auto dist = [0, 0.5, 2, 4, 8, 10, 66, 1000];
  auto rec_frac_h = [0,
                     0.004975083125415946661008,
                     0.019605280423838411518744,
                     0.038441826806682122263936,
                     0.073928105516894326854072,
                     0.090634623461009089506746,
                     0.366432349017074832087104,
                     0.499999998969423209427276];
  auto rec_frac_k = [0,
                     0.004999833339999730161263,
                     0.019989340155581784841399,
                     0.039914884555565681434341,
                     0.079324252148749455071375,
                     0.098687660112452002536543,
                     0.433391964424909348352344,
                     0.5];
  auto rec_frac_m = [0, 0.005, 0.02, 0.04, 0.08, 0.1, 0.5, 0.5];
  auto rec_frac_cf = [0,
                      0.0049999999899611822368,
                      0.0199999897600116466334,
                      0.0399996723259650951987,
                      0.0799895172931675141337,
                      0.0999680227368065638105,
                      0.4771558601757918682829,
                      0.4999999998184003757729];

  auto dist_orig = dist.dup;

  auto result_h = dist_to_recfrac(dist, GeneticMapFunc.Haldane);
  auto result_k = dist_to_recfrac(dist, GeneticMapFunc.Kosambi);
  auto result_m = dist_to_recfrac(dist, GeneticMapFunc.Morgan);
  auto result_cf = dist_to_recfrac(dist, GeneticMapFunc.Carter_Falconer);

  assert(dist == dist_orig);

  foreach(i; 0..dist.length) {
    assert(abs(result_h[i] - rec_frac_h[i]) < 1e-14);
    assert(abs(result_k[i] - rec_frac_k[i]) < 1e-14);
    assert(abs(result_m[i] - rec_frac_m[i]) < 1e-14);
    assert(abs(result_cf[i] - rec_frac_cf[i]) < 1e-13);
  }
}

// convert recombination fraction to cM distance
// Note: one chromosome at a time
double[] recfrac_to_dist(double[] rec_frac, GeneticMapFunc which_dist_to_recfrac = GeneticMapFunc.Haldane)
in {
  foreach(rf; rec_frac) {
    assert(rf >= 0.0 && rf <= 0.5, "rec_frac must be >= 0 and <= 0.5");
  }
}
body {
  auto dist_cM = new double[rec_frac.length];
  enum TOL = 1e-14;

  foreach(i, rf; rec_frac) {
    if(rf >= 0.5-TOL) { rf = 0.5 - TOL; }
    switch(which_dist_to_recfrac) {
    case GeneticMapFunc.Haldane:
      dist_cM[i] = -50*log(1-2*rf);
      break;
    case GeneticMapFunc.Kosambi:
      dist_cM[i] = 50*atanh(2*rf);
      break;
    case GeneticMapFunc.Morgan:
      dist_cM[i] = rf*100;
      break;
    case GeneticMapFunc.Carter_Falconer:
      dist_cM[i] = 12.5*(log(1+2*rf)-log(1-2*rf))+25*atan(2*rf);
    default: break;
    }
  }

  return dist_cM;
}

unittest {
  writeln("Unit test " ~ __FILE__, " : recfrac_to_dist");

  auto rec_frac = [0, 0.005, 0.02, 0.04, 0.08, 0.1, 0.4, 0.5];
  auto dist_h = [0,
                 0.5025167926750725433394,
                 2.0410997260127583530220,
                 4.1690804469525506448235,
                 8.7176693572388899156067,
                 11.1571775657104854673207,
                 80.4718956217050305212979,
                 1577.1921859393446538888384];
  auto dist_k = [0,
                 0.500016667666738134912,
                 2.001067691838410933514,
                 4.008566251879484454435,
                 8.069334806576277330237,
                 10.136627702704108955345,
                 54.930614433405487773143,
                 805.924772483670722067473];
  auto dist_m = [0, 0.5, 2, 4, 8, 10, 40, 50];
  auto dist_cf = [0,
                  0.5000000010000001937627,
                  2.0000010240014569617983,
                  4.0000327687456751490913,
                  8.0010489579481731681199,
                  10.0032028475990735216783,
                  44.3338307722915629938143,
                  422.5973403267713024433760];
  
  auto rec_frac_orig = rec_frac.dup;

  auto result_h = recfrac_to_dist(rec_frac, GeneticMapFunc.Haldane);
  auto result_k = recfrac_to_dist(rec_frac, GeneticMapFunc.Kosambi);
  auto result_m = recfrac_to_dist(rec_frac, GeneticMapFunc.Morgan);
  auto result_cf = recfrac_to_dist(rec_frac, GeneticMapFunc.Carter_Falconer);

  assert(rec_frac == rec_frac_orig);

  foreach(i; 0..rec_frac.length) {
    assert(abs(result_h[i] - dist_h[i]) < 1e-11);
    assert(abs(result_k[i] - dist_k[i]) < 1e-11);
    assert(abs(result_m[i] - dist_m[i]) < 1e-11);
    assert(abs(result_cf[i] - dist_cf[i]) < 1e-11);
  }
}
