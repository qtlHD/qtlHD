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
   paste0("Rlod = [", paste(print(out$lod, digits=20), collapse=", "), "];")

   *
   ******************************/

  auto Rlod = [1.75082264702388, 1.90066849718966, 2.01704712283993, 2.08712140364628, 2.10492717358022, 2.70630273019253, 3.37761230506607, 4.05942123556486, 4.675195738235, 5.16235351860865, 5.49575699915243, 5.68655510311692, 6.03045154338031, 6.13583408757336, 5.99580084184214, 5.64987645100132, 6.37031242544765, 6.6765466974752, 6.55134165341997, 6.06133769955557, 6.0613382616138, 6.0324327306908, 5.84745187330492, 5.78282966210315, 5.41716088137815, 4.81528999604944, 4.52366494344466, 4.00127126667383, 3.35266038068522, 3.4038725423772, 3.38560912296015, 3.29543897580413, 3.1420793627297, 3.07376784424855, 2.99184546370134, 2.90164906063143, 2.80647516232494, 2.7083487627566, 2.60898335637313];

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
   paste0("Rlod = [", paste(print(out$lod, digits=20), collapse=", "), "];")

   *
   ******************************/

  Rlod = [0.446865602965488, 0.317198313351184, 0.199081354793407, 0.10854115584192, 0.0529134102032458, 0.0499288496749841, 0.0569520588662726, 0.0748122318319702, 0.101175554914278, 0.131119885816929, 0.159462286102666, 0.182907463814672, 0.200541393269987, 0.348839229304872, 0.529112402930537, 0.714726699732694, 0.87743482697681, 0.995971219968339, 1.15220595382914, 1.34424225094529, 1.56255486107983, 1.78976615496708, 2.00406042107898, 2.18537672208657, 2.20631717307816, 2.16134284799142, 2.04393512304563, 1.86435173175937, 1.64684004291769, 1.41958398569079, 1.20499755036639];

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
   paste0("Rlod = [", paste(print(out$lod, digits=20), collapse=", "), "];")

   *
   ******************************/

  Rlod = [1.30859909425368, 1.63906716218406, 2.08903273767635, 2.52644014150326, 2.89361109626952, 3.14916508988642, 3.28383730510029, 3.52058129824172, 3.67436603879304, 3.67449449111808, 4.12741119832123, 4.38903322260359, 4.40806340509596, 3.88850409211187, 4.55939775508318, 4.62626207478718, 5.87536962625478, 4.56342478195791, 4.56366037963141, 3.9741597951986, 3.14762991240514, 2.26658650640547, 1.50094419835136];

  assert(lod.length == Rlod.length);
  foreach(i; 0..lod.length)
    assert(abs(lod[i][0] - Rlod[i]) < 1e-5,
           to!string(i) ~ "  " ~ to!string(lod[i][0]) ~ "  " ~ to!string(Rlod[i]) ~
           "  " ~ to!string(abs(lod[i][0] - Rlod[i])));

}
