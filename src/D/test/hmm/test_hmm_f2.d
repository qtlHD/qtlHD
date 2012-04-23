/**
 * test reading data and then apply HMM
 * (intercross example: listeria)
 */

module test.hmm.test_hmm_f2;

import std.math, std.stdio, std.path;
import std.exception, std.conv;
import std.random;
import std.algorithm;

import qtl.core.primitives;
import qtl.core.phenotype, qtl.core.chromosome, qtl.core.genotype;
import qtl.plugins.qtab.read_qtab;
import qtl.core.map.make_map, qtl.core.map.map;
import qtl.core.map.genetic_map_functions;
import qtl.core.hmm.f2;


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
  assert(info["Cross"] == "F2");
  assert(info["Phase"] == "unknown");

  // load symbols
  auto symbol_fn = to!string(buildPath(dir,"listeria_symbol.qtab"));
  // First read symbol information (the GenotypeCombinators)
  writeln("reading ",symbol_fn);
  auto fs = File(symbol_fn,"r");
  auto symbols = read_genotype_symbol_qtab(fs);

  // Test working of symbols
  assert(symbols.decode("A") == symbols.decode("AA"));
  assert(to!string(symbols.decode("NA")) == "[NA]");
  assert(to!string(symbols.decode("A")) == "[(0,0)]");
  assert(to!string(symbols.decode("H")) == "[(0,1), (1,0)]");
  assert(to!string(symbols.decode("B")) == "[(1,1)]");
  assert(to!string(symbols.decode("HorA")) == "[(0,0), (0,1), (1,0)]");
  assert(to!string(symbols.decode("HorB")) == "[(0,1), (1,0), (1,1)]");

  // Read genotype matrix
  auto genotype_fn = to!string(buildPath(dir,"listeria_genotype.qtab"));
  writeln("reading ",genotype_fn);
  auto fg = File(genotype_fn,"r");
  auto ret = read_genotype_qtab(fg, symbols);
  auto individuals = ret[0];
  auto genotype_matrix = ret[1];

  // Show the first individual and genotypes
  assert(individuals.list.length == 120);
  assert(individuals.list[15].name == "16");
  assert(genotype_matrix[119].length == 133);

  // by symbol
  assert(genotype_matrix[0][0] == symbols.decode("B"));
  assert(genotype_matrix[0][3] == symbols.decode("H"));

  assert(genotype_matrix[15][0] == symbols.decode("H"));
  assert(genotype_matrix[15][130] == symbols.decode("HorB"));

  assert(genotype_matrix[18][129] == symbols.decode("B"));
  assert(genotype_matrix[18][130] == symbols.decode("NA"));
  assert(genotype_matrix[18][131] == symbols.decode("A"));

  // by founders
  assert(genotype_matrix[0][0].list[0].homozygous == true);
  assert(genotype_matrix[0][0].list[0].heterozygous == false);
  assert(genotype_matrix[0][0].list[0].founders[0] == 1);
  assert(genotype_matrix[0][0].list[0].founders[1] == 1);

  assert(genotype_matrix[0][3].list[0].homozygous == false);
  assert(genotype_matrix[0][3].list[0].heterozygous == true);
  assert(genotype_matrix[0][3].list[0].founders[0] == 0);
  assert(genotype_matrix[0][3].list[0].founders[1] == 1);

  assert(genotype_matrix[15][0].list[0].homozygous == false);
  assert(genotype_matrix[15][0].list[0].heterozygous == true);
  assert(genotype_matrix[15][0].list[0].founders[0] == 0);
  assert(genotype_matrix[15][0].list[0].founders[1] == 1);

  assert(genotype_matrix[18][129].list[0].homozygous == true);
  assert(genotype_matrix[18][129].list[0].heterozygous == false);
  assert(genotype_matrix[18][129].list[0].founders[0] == 1);
  assert(genotype_matrix[18][129].list[0].founders[1] == 1);

  assert(genotype_matrix[18][131].list[0].homozygous == true);
  assert(genotype_matrix[18][131].list[0].heterozygous == false);
  assert(genotype_matrix[18][131].list[0].founders[0] == 0);
  assert(genotype_matrix[18][131].list[0].founders[1] == 0);

  // reading phenotypes
  auto pheno_fn = to!string(buildPath(dir,"listeria_phenotype.qtab"));
  writeln("reading ",pheno_fn);
  auto p_res = read_phenotype_qtab!(Phenotype!double)(pheno_fn);
  Phenotype!double[][] pheno = p_res[0];

  assert(pheno.length == 120);
  foreach(p; pheno) assert(p.length == 1);

  // 1st ind, 1st phenotype
  assert(pheno[0][0].value == 118.317);
  // 3rd ind, 1st phenotype
  assert(pheno[2][0].value == 194.917);
  assert(to!string(pheno[29][0]) == "NA"); // missing value

  // Marker map reader
  auto marker_map_fn = to!string(buildPath(dir,"listeria_marker_map.qtab"));
  writeln("reading ",marker_map_fn);
  auto markers = read_marker_map_qtab!(Marker)(marker_map_fn);
  assert(markers.length == 133);
  assert(markers[0].name == "D10M44");
  assert(markers[0].chromosome.name == "1");
  assert(markers[3].get_position == 40.4136);
  assert(markers[128].name == "D19M117");
  assert(markers[128].chromosome.name == "19");
  assert(markers[128].get_position == 16.364);

  // marker id == numeric index
  foreach(i, m; markers) assert(m.id == i);

  // test splitting up of markers into chromosomes
  //    note: chromosomes not necessarily ordered
  auto markers_by_chr = get_markers_by_chromosome(markers);
  foreach(chr; markers_by_chr) {
    // check that markers within chromosome are in order:
    //    contiguous ids; non-decreasing position
    assert(chr[1][0].chromosome.name == chr[0].name);
    for(auto i=1; i<chr[1].length; i++) {
      assert(chr[1][i].id == chr[1][i-1].id+1);
      assert(chr[1][i].get_position >= chr[1][i-1].get_position);
      assert(chr[1][i].chromosome.name == chr[0].name);
    }
  }

  auto markers_by_chr_sorted = sort_chromosomes_by_marker_id(markers_by_chr);
  Marker[] pmap_stepped, pmap_minimal;
  foreach(chr; markers_by_chr_sorted) {
    // check that markers within chromosome are in order:
    //    contiguous ids; non-decreasing position
    assert(chr[1][0].chromosome.name == chr[0].name);
    for(auto i=1; i<chr[1].length; i++) {
      assert(chr[1][i].id == chr[1][i-1].id+1);
      assert(chr[1][i].get_position >= chr[1][i-1].get_position);
      assert(chr[1][i].chromosome.name == chr[0].name);
    }

    pmap_stepped = add_stepped_markers_autosome(chr[1], 5.0, 0.0);
    pmap_minimal = add_minimal_markers_autosome(chr[1], 5.0, 0.0);
  }

  // test calc_geno_prob with listeria data
  writeln("Test calc_geno_prob with listeria data, chr 4");
  auto chr4_map = markers_by_chr_sorted[3][1];
  auto rec_frac = recombination_fractions(chr4_map, GeneticMapFunc.Haldane);

  auto genoprobs = calc_geno_prob_F2(genotype_matrix, chr4_map, rec_frac, 0.002);

  double[][] genoprobs_from_rqtl;

  /* probs from R/qtl for individual 1 */
  genoprobs_from_rqtl = [[0.99365258749878116, 0.005366774350785078, 0.00098063815043388934],
                         [0.01597337476839351, 0.984010456989931837, 0.00001616824167473149],
                         [0.98819056080900758, 0.010834274013193156, 0.00097516517779884938],
                         [0.00156449602454221, 0.998274388090304443, 0.00016111588515345984]];

  foreach(i; 0..genoprobs.length)
    foreach(j; 0..2)
      assert(abs(genoprobs[i][0][j] - genoprobs_from_rqtl[i][j]) < 1e-7);

  /* probs from R/qtl for individual 88 */
  genoprobs_from_rqtl = [[0.0001172399151822272,  0.000682931055218006, 0.999199829029599917],
                         [0.0009524442550889619,  0.058082586355804780, 0.940964969389106454],
                         [0.0001167398194463191,  0.001182526695390220, 0.998700733485163750],
                         [0.0001586426302537041,  0.998263173522426106, 0.001578183847321077]];

  foreach(i; 0..genoprobs.length)
    foreach(j; 0..2)
      assert(abs(genoprobs[i][87][j] - genoprobs_from_rqtl[i][j]) < 1e-6);


  /* probs from R/qtl for individual 103 */
  genoprobs_from_rqtl = [[0.3938894726346546, 0.467329331452057628, 0.1387811959132879691],
                         [0.4722757014939976, 0.429690530637569956, 0.0980337678684325559],
                         [0.5756148456641348, 0.365802567955339664, 0.0585825863805257627],
                         [0.9970029970029969, 0.001998001998001998, 0.0009990009990009992]];

  foreach(i; 0..genoprobs.length)
    foreach(j; 0..2)
      assert(abs(genoprobs[i][102][j] - genoprobs_from_rqtl[i][j]) < 1e-6);

  /* probs from R/qtl for individual 106 */
  genoprobs_from_rqtl = [[0.1387811959132879691, 0.467329331452057628, 0.3938894726346546],
                         [0.0980337678684325975, 0.429690530637569956, 0.4722757014939976],
                         [0.0585825863805258112, 0.365802567955339497, 0.5756148456641349],
                         [0.0009990009990009992, 0.001998001998001998, 0.9970029970029971]];

  foreach(i; 0..genoprobs.length)
    foreach(j; 0..2)
      assert(abs(genoprobs[i][105][j] - genoprobs_from_rqtl[i][j]) < 1e-6);

  /* probs from R/qtl for individual 107 */
  genoprobs_from_rqtl = [[0.2336319623541089074, 0.5327360752917820, 0.2336319623541089074],
                         [0.2147748854695731846, 0.5704502290608535, 0.2147748854695732956],
                         [0.1827669522138612168, 0.6344660955722773, 0.1827669522138613833],
                         [0.0005005005005005007, 0.9989989989989990, 0.0005005005005005007]];

  foreach(i; 0..genoprobs.length)
    foreach(j; 0..2)
      assert(abs(genoprobs[i][106][j] - genoprobs_from_rqtl[i][j]) < 1e-6);

  // test estmap with listeria data, chr 13
  writeln("Test estmap with listeria data, chr 13");

  auto chr13_map = markers_by_chr_sorted[12][1];
  rec_frac = recombination_fractions(chr13_map, GeneticMapFunc.Kosambi);
  auto rec_frac_rev = estmap_F2(genotype_matrix, chr13_map, rec_frac, 0.002, 100, 1e-6, false);

  /******************************
   * in R:

   data(listeria)
   listeria <- listeria["13",] # pull out chr 13
   listeria$geno[["13"]]$map <- as.numeric(substr(listeria$geno[["13"]]$map, 1, 7))
   names(listeria$geno[["13"]]$map) <- colnames(listeria$geno[["13"]]$data)
   rf <- mf.k(diff(unlist(pull.map(listeria))))
   map <- est.map(listeria, err=0.002, map.function="haldane", tol=1e-7)
   rfrev <- mf.h(diff(unlist(map)))
   paste0("auto rec_frac_rqtl = [", paste(sprintf("%.20f", rf), collapse=", "), "];")
   paste0("auto rec_frac_rev_rqtl = [", paste(sprintf("%.20f", rfrev), collapse=", "), "];")

   *
   ******************************/

  auto rec_frac_rqtl = [0.00286746856284019434, 0.09944712696578776601, 0.02681424948073934250, 0.00000999999999867889, 0.05831343140397290958, 0.02102459540041785868, 0.03855133275858499409, 0.01283917692671171160, 0.02231716200601452718, 0.00000999999999864336, 0.07535458078316857600];
  auto rec_frac_rev_rqtl = [0.01183048207281678366, 0.10662947366989428133, 0.03230106006565025556, 0.00000000010000000827, 0.05675828936511045919, 0.02110557519841721463, 0.03833376626430756717, 0.01427505015825597523, 0.02370366835248577386, 0.00000000009999995276, 0.08059502076315699926];

  assert(rec_frac_rqtl.length == rec_frac.length);
  foreach(i; 0..rec_frac.length) {
    assert(abs(rec_frac[i] - rec_frac_rqtl[i]) < 1e-6);
  }

  assert(rec_frac_rev_rqtl.length == rec_frac_rev.length);
  foreach(i; 0..rec_frac_rev.length) {
    assert(abs(rec_frac_rev[i] - rec_frac_rev_rqtl[i]) < 1e-7);
  }



  // test estmap with listeria data, chr 7
  writeln("Test estmap with listeria data, chr 7");

  auto chr7_map = markers_by_chr_sorted[6][1];
  rec_frac = recombination_fractions(chr7_map, GeneticMapFunc.Haldane);
  rec_frac_rev = estmap_F2(genotype_matrix, chr7_map, rec_frac, 0.01, 100, 1e-6, false);

  /******************************
   * in R:

   data(listeria)
   listeria <- listeria["7",] # pull out chr 7
   listeria$geno[["7"]]$map <- as.numeric(substr(listeria$geno[["7"]]$map, 1, 7))
   names(listeria$geno[["7"]]$map) <- colnames(listeria$geno[["7"]]$data)
   rf <- mf.h(diff(unlist(pull.map(listeria))))
   map <- est.map(listeria, err=0.01, map.function="kosambi", tol=1e-7)
   rfrev <- mf.k(diff(unlist(map)))
   paste0("rec_frac_rqtl = [", paste(sprintf("%.20f", rf), collapse=", "), "];")
   paste0("rec_frac_rev_rqtl = [", paste(sprintf("%.20f", rfrev), collapse=", "), "];")

   *
   ******************************/

  rec_frac_rqtl = [0.15661986512953757211, 0.13781102774534148558, 0.05760104042887342901, 0.15864052585018517672, 0.10645158449454483751];
  rec_frac_rev_rqtl = [0.17951889350312594251, 0.15585635633447825210, 0.06089437508606505844, 0.18207986332242809269, 0.11747108991514787490];

  assert(rec_frac_rqtl.length == rec_frac.length);
  foreach(i; 0..rec_frac.length) {
    assert(abs(rec_frac[i] - rec_frac_rqtl[i]) < 1e-6);
  }

  assert(rec_frac_rev_rqtl.length == rec_frac_rev.length);
  foreach(i; 0..rec_frac_rev.length) {
    assert(abs(rec_frac_rev[i] - rec_frac_rev_rqtl[i]) < 1e-7);
  }
}
