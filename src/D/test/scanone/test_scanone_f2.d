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
   for(i in 1:nchr(listeria)) { # round map to be like in qtab file
     listeria$geno[[i]]$map <- as.numeric(substr(listeria$geno[[i]]$map, 1, 7))
     names(listeria$geno[[i]]$map) <- colnames(listeria$geno[[i]]$data)
   }
   listeria <- calc.genoprob(listeria, step=2, stepwidth="max", err=0.01, map="kosambi")
   out <- scanone(listeria, method="hk", chr=5)
   paste0("auto Rlod = [", paste(sprintf("%.20f", out$lod), collapse=", "), "];")

   *
   ******************************/

  auto Rlod = [1.75082274582462105172, 1.90066877182243842981, 2.01704759179830261928, 2.08712206844444381204, 2.10492801585661482022, 2.70630302156502011712, 3.37761157243420484519, 4.05941906905133009786, 4.67519200527806333412, 5.16234843650192942732, 5.49575099866888194811, 5.68654862394015481186, 6.03045266166685678400, 6.13583525190000500515, 5.99580199981136274801, 5.64987753525338121108, 6.37031211391189344795, 6.67654657431808118417, 6.55134185597233908993, 6.06133825259968261889, 6.06133881500704774226, 6.03243308028930869114, 5.84745174506542753079, 5.78282955114690366827, 5.41716089002414946663, 4.81529017274021953199, 4.52366508445351200862, 4.00127134254677230274, 3.35266038968075008597, 3.40387211705490244640, 3.38560793530706405363, 3.29543670063321769703, 3.14207577905830248710, 3.07376440222537894442, 2.99184242318608539790, 2.90164662696315645007, 2.80647348734129309378, 2.70834795579105502839, 2.60898348839060645332];

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
   for(i in 1:nchr(listeria)) { # round map to be like in qtab file
     listeria$geno[[i]]$map <- as.numeric(substr(listeria$geno[[i]]$map, 1, 7))
     names(listeria$geno[[i]]$map) <- colnames(listeria$geno[[i]]$data)
   }
   listeria <- calc.genoprob(listeria, step=2, stepwidth="max", err=0.02, map="haldane")
   out <- scanone(listeria, method="hk", chr=12)
   paste0("Rlod = [", paste(sprintf("%.20f", out$lod), collapse=", "), "];")

   *
   ******************************/

  Rlod = [0.44686557370641821763, 0.31719830021603456771, 0.19908135984303498844, 0.10854117492033310555, 0.05291343558741345987, 0.04992887271254176085, 0.05695208748323921100, 0.07481227281692781617, 0.10117561287131593417, 0.13111996295117478439, 0.15946238259550682415, 0.18290757849905503463, 0.20054152440951611425, 0.34883933227411034750, 0.52911241086792415445, 0.71472657919201765253, 0.87743459204699547627, 0.99597084099229959975, 1.15220539884620620796, 1.34424151723368368039, 1.56255400476413797151, 1.78976530327673799547, 2.00405974435153666491, 2.18537637130151551901, 2.20631717509814961886, 2.16134314105408975593, 2.04393554274116695524, 1.86435209499614984452, 1.64684024437053722068, 1.41958402988689158519, 1.20499751437273516785];

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
   for(i in 1:nchr(listeria)) { # round map to be like in qtab file
     listeria$geno[[i]]$map <- as.numeric(substr(listeria$geno[[i]]$map, 1, 7))
     names(listeria$geno[[i]]$map) <- colnames(listeria$geno[[i]]$data)
   }
   listeria <- calc.genoprob(listeria, step=2.5, stepwidth="max", err=0.001, map="haldane")
   out <- scanone(listeria, method="hk", chr=13)
   paste0("Rlod = [", paste(sprintf("%.20f", out$lod), collapse=", "), "];")

   *
   ******************************/

  Rlod = [1.30859914180359737657, 1.63906665800357131957, 2.08903247621617538243, 2.52644026687539735576, 2.89361190865412254425, 3.14916697347706531218, 3.28384055386112549968, 3.52058635649603957063, 3.67437183797466104807, 3.67450028908729109389, 4.12741448044937442319, 4.38903439423523877849, 4.40806339313496664545, 3.88850410486634245899, 4.55939782295905615683, 4.62626217515548887604, 5.87536984215734037207, 4.56342204443120635915, 4.56365765088105490577, 3.97415778431093258405, 3.14762897276216335740, 2.26658651753189133160, 1.50094471689578767837];

  assert(lod.length == Rlod.length);
  foreach(i; 0..lod.length)
    assert(abs(lod[i][0] - Rlod[i]) < 1e-5,
           to!string(i) ~ "  " ~ to!string(lod[i][0]) ~ "  " ~ to!string(Rlod[i]) ~
           "  " ~ to!string(abs(lod[i][0] - Rlod[i])));

}
