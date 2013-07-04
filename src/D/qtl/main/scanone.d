// Entry point for command line (CLI) scanone tool

import std.getopt;
import std.stdio;
import std.conv;
import std.exception;
import std.file;
import std.string;
import std.path;
import std.algorithm;
import std.math;
import std.container;
import std.typecons;
import std.variant;
import std.array;
import std.range;

import qtl.plugins.csv.read_csv;
import qtl.plugins.qtab.read_qtab;
import qtl.core.chromosome;
import qtl.core.primitives;
import qtl.core.marker;
import qtl.core.genotype;
import qtl.core.phenotype;
import qtl.plugins.qtab.read_qtab;
import qtl.core.map.map;
import qtl.core.map.make_map;
import qtl.core.data.matrix;

import qtl.core.map.genetic_map_functions;
import qtl.core.hmm.cross;
import qtl.core.hmm.calcgenoprob;
import qtl.core.scanone.scanone_hk;
import qtl.core.scanone.util;

static string ver = import("VERSION");

string credits = "Karl W. Broman, Pjotr Prins and Danny Arends";
string copyright = "; qtlHD project (c) 2012-2013";
string usage = "
  usage: scanone [options] inputfile(s)

  options:

    --format          qtab|csv (default qtab)
    --phenocol        vector of numeric indices, of phenotypes to scan

  options for CSV files 

    --cross           F2|BC|RISELF|RISIB (default F2)
    --genotypes       genotype codes (default for BC is 'A H')
    --na              missing data identifiers (default '- NA')
    --sex             name of sex phenotype (default 'sex')
    --crossdir        name of cross direction phenotype (default 'pgm')

  other options:

    -v --verbosity    Set verbosity level (default 1)
    -d --debug        Set debug message level (default 0)
    --credits         Show list of contributors

  examples:

    Execute scanone with the listeria dataset, the csv version

      ./scanone --format csv ../../test/data/input/listeria.csv

    the qtab version

      ./scanone --format qtab ../../test/data/input/listeria_qtab/listeria_symbol.qtab ../../test/data/input/listeria_qtab/listeria_founder.qtab ../../test/data/input/listeria_qtab/listeria_marker_map.qtab ../../test/data/input/listeria_qtab/listeria_genotype.qtab ../../test/data/input/listeria_qtab/listeria_phenotype.qtab

    should display:

       --Peaks with LOD > 2:
       ----Chr 1  : peak for phenotype 0: max lod =    2.10 at pos =   81.40
       ----Chr 5  : peak for phenotype 0: max lod =    6.68 at pos =   27.30
       ----Chr 6  : peak for phenotype 0: max lod =    3.33 at pos =   59.37
       ----Chr 12 : peak for phenotype 0: max lod =    2.16 at pos =   43.60
       ----Chr 13 : peak for phenotype 0: max lod =    5.90 at pos =   26.16
       ----Chr 15 : peak for phenotype 0: max lod =    3.18 at pos =   22.40

";

int main(string[] args) {
  writeln("scanone ",strip(ver)," ",copyright);
  if (args.length == 1) {
    writeln(usage);
    return 0;
  }
  bool show_help = false;
  uint verbosity = 1;
  uint debug_level = 0;
  bool contributors = false;
  string format = "qtab";
  string cross = "F2";
  string genotype_ids = "A H B D C";
  string na_ids = "- NA";
  string sexcol = "sex";
  string crossdircol = "pgm";
  string phenocol = "";

  // ---- parse arguments to grab cross type
  getopt(args,
         std.getopt.config.passThrough,
         "cross", (string o, string s) { cross = s.toUpper; },
         );

  // ---- different default genotypes for other crosses
  if(cross == "BC") {
    genotype_ids = "A H B";
  }
  else if(cross == "RISIB" || cross == "RISELF") {
    genotype_ids = "A B";
  }

  // ---- parse arguments again; 'cross' has been removed
  getopt(args, "v|verbose", (string o, string v) { verbosity = to!int(v); },
               "d|debug", (string o, string d) { debug_level = to!int(d); },
               "h|help", (string o) { show_help = true; },
               "format", (string o, string s) { format = s; },
               "na", (string o, string s) { na_ids = s; },
               "genotypes", (string o, string s) { genotype_ids = s; },
               "credits", (string o) { contributors = true; },
               "phenocol", (string o, string s) { phenocol = s; },
               "sex", (string o, string s) { sexcol = s; },
               "crossdir", (string o, string s) { crossdircol = s; }
         );

  if (show_help) {
    writeln(usage);
    return 0;
  }
  if (debug_level > 0) writeln(args);
  if (contributors) {
    writeln("  by ",credits);
    return 1;
  }
  writeln("Verbosity: ",verbosity);
  writeln("Debug level: ",debug_level);

  SymbolSettings s;
  Founders f;
  Marker[] ms;
  Inds i;
  PhenotypeMatrix p; 
  // ObservedGenotypes observed;  // unused
  GenotypeMatrix g;
 
  switch(format) {
    case "qtab" :
      auto res = load_qtab(args[1..$]);
      s  = res[0];  // symbols
      f  = res[1];  // founder format (contains Cross information)
      ms = res[2];  // markers
      i  = res[3];  // individuals
      p  = res[4];  // phenotype matrix
      // observed  = res[5];  // observed genotypes
      g  = res[6];  // observed genotype matrix
      break;
    case "csv" : 
      writeln("genotype_ids: ", genotype_ids);
      auto observed_genotypes = parse_genotype_ids(cross,genotype_ids,na_ids);
      writeln("Observed ", observed_genotypes.toEncodingString(), " ", observed_genotypes);
      auto res = load_csv(args[1], observed_genotypes);
      f["Cross"] = cross;
      ms = res[0];  // markers
      i  = res[1];  // individuals
      p  = res[2];  // phenotype matrix
      g  = res[4];
      // observed = observed_genotypes;  // unused
      break;
    default :
      throw new Exception("Unknown format "~format);
  }
  

  if (debug_level > 2) {
    writeln("* Format");
    writeln(format);
    writeln("* Symbol data");
    writeln(s);
    writeln("* Individuals");
    writeln(i);
    // writeln("* Observed genotypes");
    // writeln(observed);
    writeln("* Genotype data (partial)");
    writeln(g[0..3]);
    writeln("* Phenotype data (partial)");
    writeln(p[0..3]);
    writeln("* Marker data (partial)");
    writeln(ms[0..3]);
  }

  // ---- Find individuals missing phenotype
  auto ind_missing_a_phenotype = 
    test_matrix_by_row_element!Phenotype(p, element => isNA(element)).array();
  writeln("Omitting ", count!"a[0]==true"(ind_missing_a_phenotype), " individuals with missing phenotype");
  if (g.length != p.length) 
    throw new Exception("Genotype individuals does not match phenotype individuals");

  // ---- drop individuals missing phenotype from phenotype list
  auto ind_to_include = map!"a[1]"( filter!"a[0]==false"(ind_missing_a_phenotype) ).array();
  auto pheno = indexed(p,ind_to_include).array();
  writeln("done omitting from phenotypes");

  //      and do the same for genotype list
  auto genotype_matrix = indexed(g,ind_to_include).array();
  assert(genotype_matrix.length == pheno.length);
  writeln("done omitting from genotypes");

  // ---- set the cross type
  auto cross_class = form_cross(f["Cross"]);
  writeln("formed cross class, ", f["Cross"]);

  auto markers_by_chr = sort_chromosomes_by_marker_id(get_markers_by_chromosome(ms));

  // ---- add pseudomarkers at 2.0 cM spacing
  auto pmar_by_chr = add_minimal_markers(markers_by_chr, 2.0);

  // ---- inter-marker recombination fractions
  auto rec_frac = recombination_fractions(pmar_by_chr, GeneticMapFunc.Haldane);

  // ---- empty covariate matrices
  auto addcovar = new double[][](0, 0);
  auto intcovar = new double[][](0, 0);
  auto weights = new double[](0);

  // ---- null model
  auto rss0 = scanone_hk_null(pheno, addcovar, weights);

  // ---- storage for LOD curves and peaks for all chromosomes
  double[][] lod;
  Tuple!(double, Marker)[][] peaks;

  // ---- calc genoprob for each chromosome, followed by scanone
  foreach(j, chr; pmar_by_chr) {
    auto genoprobs = calc_geno_prob(cross_class, genotype_matrix, chr[1], rec_frac[j][0], 0.002);
    auto rss = scanone_hk(genoprobs, pheno, addcovar, intcovar, weights);
    auto lod_this_chr = rss_to_lod(rss, rss0, pheno.length);
    lod ~= lod_this_chr;

    auto peak_this_chr = get_peak_scanone(lod_this_chr, chr[1]);
    peaks ~= peak_this_chr;
  }

  // ---- print peaks
  double threshold = 2;
  writeln(" --Peaks with LOD > ", threshold, ":");
  foreach(peak; peaks) {
    foreach(j; 0..peak.length) {
      if(peak[j][0] > threshold)
        writefln(" ----Chr %-2s : peak for phenotype %d: max lod = %7.2f at pos = %7.2f", peak[j][1].chromosome.name, j,
                 peak[j][0], peak[j][1].get_position);
    }
  }
  return 0;
}
