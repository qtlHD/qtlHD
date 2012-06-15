/**
 * Test scanone routines: intercross (listeria data)
 */

module test.scanone.test_scanone_f2;

import std.stdio;
import std.conv;
import std.string;
import std.path;
import std.algorithm;
import std.math;
import std.algorithm;

import qtl.core.primitives;
import qtl.core.chromosome;
import qtl.core.phenotype;
import qtl.core.marker;
import qtl.core.genotype;
import qtl.core.phenotype;
import qtl.plugins.qtab.read_qtab;
import qtl.core.map.map;
import qtl.core.map.make_map;

import qtl.core.map.genetic_map_functions;
import qtl.core.hmm.cross;
import qtl.core.hmm.calcgenoprob;
import qtl.core.scanone.scanone_hk;
import qtl.core.scanone.util;
import qtl.core.util.data_manip;


unittest {
  writeln("Unit test " ~ __FILE__);

  alias std.path.buildPath buildPath;
  auto dir = to!string(dirName(__FILE__) ~ dirSeparator ~
                       buildPath("..","..","..","..","test","data", "input", "listeria_qtab"));

  // load founder info
  auto founder_fn = to!string(buildPath(dir, "listeria_founder.qtab"));
  writeln("reading ", founder_fn);
  auto info = read_founder_settings_qtab(founder_fn);

  // load symbols
  auto symbol_fn = to!string(buildPath(dir,"listeria_symbol.qtab"));
  // First read symbol information (the GenotypeCombinators)
  writeln("reading ",symbol_fn);
  auto fs = File(symbol_fn,"r");
  auto symbols = read_genotype_symbol_qtab(fs);

  // Read genotype matrix
  auto genotype_fn = to!string(buildPath(dir,"listeria_genotype.qtab"));
  writeln("reading ",genotype_fn);
  auto fg = File(genotype_fn,"r");
  auto ret = read_genotype_qtab(fg, symbols);
  auto individuals = ret[0];
  auto genotype_matrix = ret[1];

  // reading phenotypes
  auto pheno_fn = to!string(buildPath(dir,"listeria_phenotype.qtab"));
  writeln("reading ",pheno_fn);
  auto p_res = read_phenotype_qtab!(Phenotype!double)(pheno_fn);
  Phenotype!double[][] pheno = p_res[0];

  auto ind_to_omit = is_any_phenotype_missing(pheno);
  auto n_to_omit = count(ind_to_omit, true);
  writeln("Omitting ", n_to_omit, " individuals with missing phenotype");

  genotype_matrix = omit_ind_from_genotypes(genotype_matrix, ind_to_omit);
  pheno = omit_ind_from_phenotypes(pheno, ind_to_omit);

  // Marker map reader
  auto marker_map_fn = to!string(buildPath(dir,"listeria_marker_map.qtab"));
  writeln("reading ",marker_map_fn);
  auto markers = read_marker_map_qtab!(Marker)(marker_map_fn);

  // split up markers into chromosomes
  auto markers_by_chr = get_markers_by_chromosome(markers);

  // sorted chromosomes
  auto markers_by_chr_sorted = sort_chromosomes_by_marker_id(markers_by_chr);

  // chr 5
  auto chr5_map = markers_by_chr_sorted[4][1];
  sort(chr5_map); // sort in place

  // add pseudomarkers
  auto chr5_map_wpmark = add_minimal_markers_autosome(chr5_map, 2.0);

  // form cross
  auto f2 = form_cross("F2");
  
  // calc_genoprob
  auto rec_frac = recombination_fractions(chr5_map_wpmark, GeneticMapFunc.Kosambi);
  auto chr5probs = calc_geno_prob(f2, genotype_matrix, chr5_map_wpmark, rec_frac, 0.01);

  // empty covariate matrices
  auto addcovar = new double[][](0, 0);
  auto intcovar = new double[][](0, 0);
  auto weights = new double[](0);

  // run scanone and calculate LOD scores
  auto rss = scanone_hk(chr5probs, pheno, addcovar, intcovar, weights);
  auto rss0 = scanone_hk_null(pheno, addcovar, weights);
  auto lod = rss_to_lod(rss, rss0, pheno.length);

  auto peak = get_peak_scanone(lod, chr5_map_wpmark);
  foreach(i; 0..peak.length) {
    writefln("Peak for phenotype %d: max lod = %7.2f at pos = %7.2f", i,
             peak[i][0], peak[i][1].get_position);
  }

  /*******************************
   * in R:

   data(listeria)
   for(i in seq(along=listeria$geno)) {
     m <- listeria$geno[[i]]$map
     listeria$geno[[i]]$map <- round(m, ifelse(m < 10, 5, 4))
   }
   listeria$geno[[5]]$map[c(3,8)] <- listeria$geno[[5]]$map[c(3,8)] - 0.0001
   listeria <- calc.genoprob(listeria, step=2, stepwidth="max", err=0.01, map="kosambi")
   out <- scanone(listeria, method="hk", chr=5)
   paste0("auto Rlod = [", paste(sprintf("%.20f", out$lod), collapse=", "), "];")

   *
   ******************************/

  auto Rlod = [1.75082274582229047155, 1.90066877180601068176, 2.01704759176658399156, 2.08712206839663849678, 2.10492801579255228717, 2.70630302141955780826, 3.37761157217386198681, 4.05941906864683232925, 4.67519200471798512808, 5.16234843579849211892, 5.49575099785062093360, 5.68654862303930030976, 6.03045265834794008697, 6.13583424663443111058, 5.99580013107305376252, 5.64987507857529180910, 6.37031259593959475751, 6.67654731409692203670, 6.55134209099810504995, 6.06133741741672338321, 6.06133792672636673160, 6.03243655551284518879, 5.84745205087784825082, 5.78282981073931523497, 5.41716092472790933243, 4.81528990129646672358, 4.52366484133540325274, 4.00127114746788947741, 3.35266025459043248702, 3.40387261361178161678, 3.38560952847876706073, 3.29543985511463688454, 3.14208080912311515931, 3.07376939057337494887, 2.99184705489955149460, 2.90165064079496914928, 2.80647666224700742532, 2.70835009464775566812, 2.60898442489252602172];

  assert(lod.length == Rlod.length);
  foreach(i; 0..lod.length) {
    assert(abs(lod[i][0] - Rlod[i]) < 1e-10,
           to!string(i) ~ "  " ~ to!string(lod[i][0]) ~ "  " ~ to!string(Rlod[i]) ~
           "  " ~ to!string(log10(abs(lod[i][0] - Rlod[i]))));
  }

  // chr 12
  auto chr12_map = markers_by_chr_sorted[11][1];
  sort(chr12_map); // sort in place

  // add pseudomarkers
  auto chr12_map_wpmark = add_minimal_markers_autosome(chr12_map, 2.0);

  // calc_genoprob
  rec_frac = recombination_fractions(chr12_map_wpmark, GeneticMapFunc.Haldane);
  auto chr12probs = calc_geno_prob(f2, genotype_matrix, chr12_map_wpmark, rec_frac, 0.02);

  // run scanone and calculate LOD scores
  rss = scanone_hk(chr12probs, pheno, addcovar, intcovar, weights);
  rss0 = scanone_hk_null(pheno, addcovar, weights);
  lod = rss_to_lod(rss, rss0, pheno.length);

  peak = get_peak_scanone(lod, chr12_map_wpmark);
  foreach(i; 0..peak.length) {
    writefln("Peak for phenotype %d: max lod = %7.2f at pos = %7.2f", i,
             peak[i][0], peak[i][1].get_position);
  }

  /*******************************
   * in R:

   data(listeria)
   for(i in seq(along=listeria$geno)) {
     m <- listeria$geno[[i]]$map
     listeria$geno[[i]]$map <- round(m, ifelse(m < 10, 5, 4))
   }
   listeria <- calc.genoprob(listeria, step=2, stepwidth="max", err=0.02, map="haldane")
   out <- scanone(listeria, method="hk", chr=12)
   paste0("Rlod = [", paste(sprintf("%.20f", out$lod), collapse=", "), "];")

   *
   ******************************/

  Rlod = [0.44686557371397839233, 0.31719830021842199130, 0.19908135983155261783, 0.10854117488901238175, 0.05291343553381011588, 0.04992887275039947781, 0.05695208741678925435, 0.07481227234819698424, 0.10117561165805000201, 0.13111996071154408128, 0.15946237918950600942, 0.18290757393202738967, 0.20054151877616277488, 0.34883936750361499435, 0.52911249767839763081, 0.71472671490977290887, 0.87743476263176489738, 0.99597119444501913677, 1.15220600484849455825, 1.34424240068182143659, 1.56255509077720944333, 1.78976639152438110614, 2.00406055549115080794, 2.18537665668486624782, 2.20631684624248691762, 2.16134232404624526680, 2.04393453744398811978, 1.86435123038620531588, 1.64683971340781454273, 1.41958383376515939744, 1.20499752515956970456];

  assert(lod.length == Rlod.length);
  foreach(i; 0..lod.length)
    assert(abs(lod[i][0] - Rlod[i]) < 1e-10,
           to!string(i) ~ "  " ~ to!string(lod[i][0]) ~ "  " ~ to!string(Rlod[i]) ~
           "  " ~ to!string(log10(abs(lod[i][0] - Rlod[i]))));

  // chr 13
  auto chr13_map = markers_by_chr_sorted[12][1];
  sort(chr13_map); // sort in place

  // add pseudomarkers
  auto chr13_map_wpmark = add_minimal_markers_autosome(chr13_map, 2.5);

  // calc_genoprob
  rec_frac = recombination_fractions(chr13_map_wpmark, GeneticMapFunc.Haldane);
  auto chr13probs = calc_geno_prob(f2, genotype_matrix, chr13_map_wpmark, rec_frac, 0.001);

  // run scanone and calculate LOD scores
  rss = scanone_hk(chr13probs, pheno, addcovar, intcovar, weights);
  rss0 = scanone_hk_null(pheno, addcovar, weights);
  lod = rss_to_lod(rss, rss0, pheno.length);

  peak = get_peak_scanone(lod, chr13_map_wpmark);
  foreach(i; 0..peak.length) {
    writefln("Peak for phenotype %d: max lod = %7.2f at pos = %7.2f", i,
             peak[i][0], peak[i][1].get_position);
  }

  /*******************************
   * in R:

   data(listeria)
   for(i in seq(along=listeria$geno)) {
     m <- listeria$geno[[i]]$map
     listeria$geno[[i]]$map <- round(m, ifelse(m < 10, 5, 4))
   }
   listeria <- calc.genoprob(listeria, step=2.5, stepwidth="max", err=0.001, map="haldane")
   out <- scanone(listeria, method="hk", chr=13)
   paste0("Rlod = [", paste(sprintf("%.20f", out$lod), collapse=", "), "];")

   *
   ******************************/

  Rlod = [1.30859906993367758332, 1.63906677087584284891, 2.08903219504446724386, 2.52643940519760690222, 2.89361004827401302464, 3.14916355799124403347, 3.28383514339594739795, 3.52057799439478458225, 3.67436228566020872677, 3.67449073914536938901, 4.12740908043326726329, 4.38903247282928532513, 4.40806342317449662005, 3.88850412574578285785, 4.55939781766801388585, 4.62626215214106650819, 5.87536984413850404962, 4.56342340093880238783, 4.56365900796430423725, 3.97415853535017049580, 3.14762891112053466713, 2.26658588253167181392, 1.50094396044403310952];

  assert(lod.length == Rlod.length);
  foreach(i; 0..lod.length)
    assert(abs(lod[i][0] - Rlod[i]) < 1e-10,
           to!string(i) ~ "  " ~ to!string(lod[i][0]) ~ "  " ~ to!string(Rlod[i]) ~
           "  " ~ to!string(log10(abs(lod[i][0] - Rlod[i]))));


  // Additive allele model

  // collapse genotype probabilities to allele probabilities
  auto chr5_allele_probs = collapse_geno_prob_to_allele_prob(chr5probs, f2.all_true_geno);

  // run scanone and calculate LOD scores
  rss = scanone_hk(chr5_allele_probs, pheno, addcovar, intcovar, weights);
  rss0 = scanone_hk_null(pheno, addcovar, weights);
  lod = rss_to_lod(rss, rss0, pheno.length);

  peak = get_peak_scanone(lod, chr5_map_wpmark);
  foreach(i; 0..peak.length) {
    writefln("Additive allele model: Peak for phenotype %d: max lod = %7.2f at pos = %7.2f", i,
             peak[i][0], peak[i][1].get_position);
  }

  // collapse genotype probabilities to allele probabilities
  auto chr13_allele_probs = collapse_geno_prob_to_allele_prob(chr13probs, f2.all_true_geno);

  // run scanone and calculate LOD scores
  rss = scanone_hk(chr13_allele_probs, pheno, addcovar, intcovar, weights);
  rss0 = scanone_hk_null(pheno, addcovar, weights);
  lod = rss_to_lod(rss, rss0, pheno.length);

  peak = get_peak_scanone(lod, chr13_map_wpmark);
  foreach(i; 0..peak.length) {
    writefln("Additive allele model: Peak for phenotype %d: max lod = %7.2f at pos = %7.2f", i,
             peak[i][0], peak[i][1].get_position);
  }
}
