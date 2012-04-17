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

   *
   ******************************/

    auto Rlod = [1.75082264702388, 1.90066849718966, 2.01704712283993, 2.08712140364622, 2.10492717358017, 2.70630273019253,
                 3.37761230506607, 4.05942123556486, 4.67519573823500, 5.16235351860865, 5.49575699915243, 5.68655510311692,
                 6.03045154338031, 6.13583408757336, 5.99580084184214, 5.64987645100132, 6.37031242544765, 6.67654669747520,
                 6.55134165341997, 6.06133769955562, 6.06133826161380, 6.03243273069080, 5.84745187330492, 5.78282966210315,
                 5.41716088137809, 4.81528999604944, 4.52366494344466, 4.00127126667383, 3.35266038068522, 3.40387254237714,
                 3.38560912296015, 3.29543897580413, 3.14207936272970, 3.07376784424855, 2.99184546370128, 2.90164906063137,
                 2.80647516232494, 2.70834876275660, 2.60898335637313];

  assert(lod.length == Rlod.length);
  foreach(i; 0..lod.length)
    assert(abs(lod[i][0] - Rlod[i]) < 1e-5);


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

   *
   ******************************/

  Rlod = [0.4468656029654881, 0.3171983133512413, 0.1990813547934067, 0.1085411558419196, 0.0529134102031890, 0.0499288496749841,
          0.0569520588662726, 0.0748122318320270, 0.1011755549142777, 0.1311198858168723, 0.1594622861026096, 0.1829074638146722,
          0.2005413932699867, 0.3488392293048719, 0.5291124029305365, 0.7147266997326369, 0.8774348269768666, 0.9959712199683395,
          1.1522059538291387, 1.3442422509452854, 1.5625548610798319, 1.7897661549670261, 2.0040604210789752, 2.1853767220865734,
          2.2063171730782187, 2.1613428479913637, 2.0439351230456282, 1.8643517317593705, 1.6468400429176882, 1.4195839856907924,
          1.2049975503663291];

  assert(lod.length == Rlod.length);
  foreach(i; 0..lod.length)
    assert(abs(lod[i][0] - Rlod[i]) < 1e-6);

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
   for(i in seq(along=listeria$geno))
     listeria$geno[[i]]$map <- round(listeria$geno[[i]]$map, 4)
   listeria <- calc.genoprob(listeria, step=2.5, stepwidth="max", err=0.001, map="haldane")
   out <- scanone(listeria, method="hk", chr=13)

   *
   ******************************/

  Rlod = [1.30859426838913, 1.63907750026016, 2.08904295511854, 2.52644926134673, 2.89361799872603, 3.14916884205684,
          3.28383743561648, 3.52057919680880, 3.67436254196735, 3.67449099540994, 4.12740932094681, 4.38903263328922,
          4.40806342365391, 3.88850412605399, 4.55939781788635, 4.62626215214118, 5.87536984413850, 4.56342340093875,
          4.56365900796430, 3.97415853535017, 3.14762891112053, 2.26658588253167, 1.50094396044409];

  assert(lod.length == Rlod.length);
  foreach(i; 0..lod.length)
    assert(abs(lod[i][0] - Rlod[i]) < 5e-4);

}
