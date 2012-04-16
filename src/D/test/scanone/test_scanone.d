/**
 * Test scanone routines
 */

module test.scanone.test_scanone;

import std.stdio;
import std.conv;
import std.string;
import std.path;
import std.algorithm;
import std.math;
alias std.algorithm.find find;

import qtl.core.primitives;
import qtl.core.chromosome;
import qtl.core.phenotype;
import qtl.core.marker;
import qtl.core.genotype;
import qtl.core.phenotype;
import qtl.plugins.qtab.read_qtab;
import qtl.core.map.map;
import qtl.core.map.make_map;
import qtl.core.scanone.scanone_hk;

import qtl.core.map.genetic_map_functions;
import qtl.core.hmm.bc;
import qtl.core.scanone.scanone_hk;

unittest {
  writeln("Unit test " ~ __FILE__);
}

unittest {
  alias std.path.buildPath buildPath;
  auto dir = to!string(dirName(__FILE__) ~ sep ~
                       buildPath("..","..","..","..","test","data", "input", "hyper_noX_qtab"));

  // Read founder info
  auto founder_fn = to!string(buildPath(dir, "hyper_noX_founder.qtab"));
  writeln("reading ", founder_fn);
  auto info = read_founder_settings_qtab(founder_fn);
  assert(info["Cross"] == "BC");

  // First read symbol information (the GenotypeCombinators)
  auto symbol_fn = to!string(buildPath(dir,"hyper_noX_symbol.qtab"));
  writeln("reading ",symbol_fn);
  auto fs = File(symbol_fn,"r");
  auto symbols = read_genotype_symbol_qtab(fs);

  // Read genotype matrix
  auto genotype_fn = to!string(buildPath(dir,"hyper_noX_genotype.qtab"));
  writeln("reading ",genotype_fn);
  auto fg = File(genotype_fn,"r");
  auto ret = read_genotype_qtab(fg, symbols);
  auto individuals = ret[0];
  auto genotype_matrix = ret[1];

  // reading phenotypes
  auto pheno_fn = to!string(buildPath(dir,"hyper_noX_phenotype.qtab"));
  writeln("reading ",pheno_fn);
  auto p_res = read_phenotype_qtab!(Phenotype!double)(pheno_fn);
  Phenotype!double[][] pheno = p_res[0];

  // Marker map reader
  auto marker_map_fn = to!string(buildPath(dir,"hyper_noX_marker_map.qtab"));
  writeln("reading ",marker_map_fn);
  auto markers = read_marker_map_qtab!(Marker)(marker_map_fn);

  // split up markers into chromosomes
  auto markers_by_chr = get_markers_by_chromosome(markers);

  // sorted chromosomes
  auto markers_by_chr_sorted = sort_chromosomes_by_marker_id(markers_by_chr);

  // chr 4, no pseudomarkers
  auto chr4_map = markers_by_chr_sorted[3][1];
  sort(chr4_map); // sort in place

  writeln(" --Scanone for hyper chr 4, no pseudomarkers");
  // calc_geno_prob
  auto rec_frac = recombination_fractions(chr4_map, GeneticMapFunc.Carter_Falconer);
  auto chr4probs = calc_geno_prob_BC(genotype_matrix, chr4_map, rec_frac, 0.001);

  // empty covariate matrices
  auto addcovar = new double[][](0, 0);
  auto intcovar = new double[][](0, 0);
  auto weights = new double[](0);

  // pull out the first phenotype; also do it as sqrt
  auto pheno_rev = pheno.dup;
  foreach(i; 0..pheno_rev.length) {
    pheno_rev[i][0].value = pheno[i][0].value;
    pheno_rev[i][1].value = sqrt(pheno[i][0].value);
  }

  // run scanone and calculate LOD scores
  auto rss = scanone_hk(chr4probs, pheno_rev, addcovar, intcovar, weights);
  auto rss0 = scanone_hk_null(pheno_rev, addcovar, weights);
  auto lod = rss_to_lod(rss, rss0, pheno_rev.length);

  /*******************************
   * in R:

   data(hyper)
   hyper <- hyper["-X",]
   hyper$pheno[,2] <- sqrt(hyper$pheno[,1])
   hyper <- calc.genoprob(hyper, err=0.001, map="c-f")
   out <- scanone(hyper, phe=1:2, method="hk", chr=4)

   *
   ******************************/

  auto Rlod = [[2.6254865829665732235, 2.7020254097664917481],
               [5.4468751777026227501, 5.5347680212903469510],
               [5.4855367913328336726, 5.5335861088395006391],
               [6.5521802135129973976, 6.5926615345241259547],
               [6.5145853225407108766, 6.5538831293538066802],
               [6.8321199282964926169, 6.8838582280987736794],
               [5.8424139084056605498, 5.8504519713298464012],
               [5.8424139084054331761, 5.8504519713297042927],
               [6.2766093214407874257, 6.2887786121572162301],
               [6.3019283482516357253, 6.3261703036766050445],
               [8.0936828368751321250, 8.0813636210099275559],
               [6.3856324540660125422, 6.3630175147217187259],
               [5.1433416237845221985, 5.1180677048473910418],
               [5.1433416237839537644, 5.1180677048468510293],
               [4.8915420586123445901, 4.8718203736263774317],
               [4.7662276124896152396, 4.7568692716767770889],
               [3.7425634657936370786, 3.7194855889463269705],
               [2.7920054024624505473, 2.7570138296044888193],
               [2.4484441531932361613, 2.4159175588472407981],
               [2.9348031147566189247, 2.9540248545710312555]];

  foreach(i; 0..lod.length)
    foreach(j; 0..lod[0].length)
      assert(abs(lod[i][j] - Rlod[i][j]) < 1e-8);

  writeln(" --Scanone for hyper chr 4, with pseudomarkers");
  auto pmap_stepped_chr4 = add_stepped_markers_autosome(chr4_map, 2.4, 0.0);
}