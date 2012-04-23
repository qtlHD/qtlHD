/**
 * Test scanone routines
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
import qtl.core.hmm.f2;
import qtl.core.scanone.scanone_hk;
import qtl.core.util.data_manip;


unittest {
  writeln("Unit test " ~ __FILE__);
}

unittest {
  alias std.path.buildPath buildPath;
  auto dir = to!string(dirName(__FILE__) ~ sep ~
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

  // calc_genoprob
  auto rec_frac = recombination_fractions(chr5_map_wpmark, GeneticMapFunc.Kosambi);
  auto chr5probs = calc_geno_prob_F2(genotype_matrix, chr5_map_wpmark, rec_frac, 0.01);

  // empty covariate matrices
  auto addcovar = new double[][](0, 0);
  auto intcovar = new double[][](0, 0);
  auto weights = new double[](0);

  // run scanone and calculate LOD scores
  auto rss = scanone_hk(chr5probs, pheno, addcovar, intcovar, weights);
  auto rss0 = scanone_hk_null(pheno, addcovar, weights);
  auto lod = rss_to_lod(rss, rss0, pheno.length);

  /*******************************
   * in R:

   data(listeria)
   listeria <- listeria["-X",]
   listeria <- calc.genoprob(listeria, step=2, stepwidth="max", err=0.01, map="kosambi")
   out <- scanone(listeria, method="hk", chr=5)
   paste0("auto Rlod = [", paste(sprintf("%.20f", out$lod), collapse=", "), "];")

   *
   ******************************/

  auto Rlod = [1.75082264702388101796, 1.90066849718965613647, 2.01704712283992648736, 2.08712140364627884992, 2.10492717358022218832, 2.70630273019253309030, 3.37761230506606580093, 4.05942123556485512381, 4.67519573823500422804, 5.16235351860865421258, 5.49575699915243376381, 5.68655510311691614334, 6.03045154338030897634, 6.13583408757335746486, 5.99580084184214001652, 5.64987645100131885556, 6.37031242544765063940, 6.67654669747520301826, 6.55134165341996776988, 6.06133769955556545028, 6.06133826161379829500, 6.03243273069080032656, 5.84745187330491944522, 5.78282966210315407807, 5.41716088137815177106, 4.81528999604944374369, 4.52366494344465763788, 4.00127126667382526648, 3.35266038068522220783, 3.40387254237720071615, 3.38560912296014748790, 3.29543897580413158721, 3.14207936272970300706, 3.07376784424855031830, 2.99184546370133830351, 2.90164906063142780113, 2.80647516232494353972, 2.70834876275660008105, 2.60898335637313039115];

  assert(lod.length == Rlod.length);
  foreach(i; 0..lod.length)
    assert(abs(lod[i][0] - Rlod[i]) < 1e-5,
           to!string(i) ~ "  " ~ to!string(lod[i][0]) ~ "  " ~ to!string(Rlod[i]) ~
           "  " ~ to!string(abs(lod[i][0] - Rlod[i])));

  // chr 12
  auto chr12_map = markers_by_chr_sorted[11][1];
  sort(chr12_map); // sort in place

  // add pseudomarkers
  auto chr12_map_wpmark = add_minimal_markers_autosome(chr12_map, 2.0);

  // calc_genoprob
  rec_frac = recombination_fractions(chr12_map_wpmark, GeneticMapFunc.Haldane);
  auto chr12probs = calc_geno_prob_F2(genotype_matrix, chr12_map_wpmark, rec_frac, 0.02);

  // run scanone and calculate LOD scores
  rss = scanone_hk(chr12probs, pheno, addcovar, intcovar, weights);
  rss0 = scanone_hk_null(pheno, addcovar, weights);
  lod = rss_to_lod(rss, rss0, pheno.length);

  /*******************************
   * in R:

   data(listeria)
   listeria <- listeria["-X",]
   listeria <- calc.genoprob(listeria, step=2, stepwidth="max", err=0.02, map="haldane")
   out <- scanone(listeria, method="hk", chr=12)
   paste0("Rlod = [", paste(sprintf("%.20f", out$lod), collapse=", "), "];")

   *
   ******************************/

  Rlod = [0.44686560296548805127, 0.31719831335118442439, 0.19908135479340671736, 0.10854115584191958988, 0.05291341020324580313, 0.04992884967498412152, 0.05695205886627263681, 0.07481223183197016624, 0.10117555491427765446, 0.13111988581692912703, 0.15946228610266643955, 0.18290746381467215542, 0.20054139326998665638, 0.34883922930487187841, 0.52911240293053651840, 0.71472669973269375987, 0.87743482697680974525, 0.99597121996833948288, 1.15220595382913870708, 1.34424225094528537738, 1.56255486107983188049, 1.78976615496708291175, 2.00406042107897519600, 2.18537672208657340889, 2.20631717307816188622, 2.16134284799142051270, 2.04393512304562818827, 1.86435173175937052292, 1.64684004291768815165, 1.41958398569079236040, 1.20499755036638589445];

  assert(lod.length == Rlod.length);
  foreach(i; 0..lod.length)
    assert(abs(lod[i][0] - Rlod[i]) < 1e-5,
           to!string(i) ~ "  " ~ to!string(lod[i][0]) ~ "  " ~ to!string(Rlod[i]) ~
           "  " ~ to!string(abs(lod[i][0] - Rlod[i])));

  // chr 13
  auto chr13_map = markers_by_chr_sorted[12][1];
  sort(chr13_map); // sort in place

  // add pseudomarkers
  auto chr13_map_wpmark = add_minimal_markers_autosome(chr13_map, 2.5);

  // calc_genoprob
  rec_frac = recombination_fractions(chr13_map_wpmark, GeneticMapFunc.Haldane);
  auto chr13probs = calc_geno_prob_F2(genotype_matrix, chr13_map_wpmark, rec_frac, 0.001);

  // run scanone and calculate LOD scores
  rss = scanone_hk(chr13probs, pheno, addcovar, intcovar, weights);
  rss0 = scanone_hk_null(pheno, addcovar, weights);
  lod = rss_to_lod(rss, rss0, pheno.length);

  /*******************************
   * in R:

   data(listeria)
   listeria <- listeria["-X",]
   listeria <- calc.genoprob(listeria, step=2.5, stepwidth="max", err=0.001, map="haldane")
   out <- scanone(listeria, method="hk", chr=13)
   paste0("Rlod = [", paste(sprintf("%.20f", out$lod), collapse=", "), "];")

   *
   ******************************/

  Rlod = [1.30859909425367959557, 1.63906716218406245389, 2.08903273767634800606, 2.52644014150325801893, 2.89361109626952384133, 3.14916508988642362965, 3.28383730510029181460, 3.52058129824172283406, 3.67436603879303902431, 3.67449449111808235102, 4.12741119832122649314, 4.38903322260358663698, 4.40806340509595884214, 3.88850409211187297842, 4.55939775508318234643, 4.62626207478717788035, 5.87536962625478054179, 4.56342478195790590689, 4.56366037963141479850, 3.97415979519860229630, 3.14762991240513656521, 2.26658650640547421062, 1.50094419835136250185];

  assert(lod.length == Rlod.length);
  foreach(i; 0..lod.length)
    assert(abs(lod[i][0] - Rlod[i]) < 1e-5,
           to!string(i) ~ "  " ~ to!string(lod[i][0]) ~ "  " ~ to!string(Rlod[i]) ~
           "  " ~ to!string(abs(lod[i][0] - Rlod[i])));

}
