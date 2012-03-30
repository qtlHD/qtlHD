/**
 * test reading data and then apply HMM
 * (backcross example: hyper with no X chr)
 */

module test.hmm.test_hmm_bc;

import std.math, std.stdio, std.path;
import std.exception, std.conv;
import std.random;
import std.algorithm;

import qtl.core.primitives;
import qtl.core.phenotype, qtl.core.chromosome, qtl.core.genotype;
import qtl.plugins.qtab.read_qtab;
import qtl.core.map.make_map, qtl.core.map.map;
import qtl.core.map.genetic_map_functions;
import qtl.core.hmm.hmm_bc;

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

  // Test working of symbols
  assert(to!string(symbols.decode("NA")) == "[NA]");
  assert(to!string(symbols.decode("A")) == "[(0,0)]");
  assert(to!string(symbols.decode("H")) == "[(1,0)]");

  // Read genotype matrix
  auto genotype_fn = to!string(buildPath(dir,"hyper_noX_genotype.qtab"));
  writeln("reading ",genotype_fn);
  auto fg = File(genotype_fn,"r");
  auto ret = read_genotype_qtab(fg, symbols);
  auto individuals = ret[0];
  auto genotype_matrix = ret[1];

  // Show the first individual and genotypes
  assert(individuals.list.length == 250);
  assert(individuals.list[244].name == "245");
  assert(genotype_matrix[244].length == 170);

  // by symbol
  assert(genotype_matrix[239][0] == symbols.decode("NA"));
  assert(genotype_matrix[239][2] == symbols.decode("A"));
  assert(genotype_matrix[239][37] == symbols.decode("H"));

  assert(genotype_matrix[244][0] == symbols.decode("NA"));
  assert(genotype_matrix[244][2] == symbols.decode("H"));
  assert(genotype_matrix[244][169] == symbols.decode("NA"));

  writeln("genotype_matrix[244][0]: ", genotype_matrix[244][0]);
  writeln("genotype_matrix[244][0].length: ", genotype_matrix[244][0].length);
  writeln("genotype_matrix[244][0].list.length: ", genotype_matrix[244][0].list.length);
  foreach(g; genotype_matrix[244][0].list)
    writeln("genotype_matrix[244][0].list[i]: ", g, " ", typeid(g));

  assert(genotype_matrix[19][0] == symbols.decode("H"));
  assert(genotype_matrix[19][2] == symbols.decode("A"));
  assert(genotype_matrix[19][3] == symbols.decode("NA"));
  assert(genotype_matrix[19][169] == symbols.decode("A"));

  writeln("genotype_matrix[19][2]: ", genotype_matrix[19][2]);
  writeln("genotype_matrix[19][2].length: ", genotype_matrix[19][2].length);
  writeln("genotype_matrix[19][2].list.length: ", genotype_matrix[19][2].list.length);
  foreach(g; genotype_matrix[19][2].list)
    writeln("genotype_matrix[19][2].list[i]: ", g, " ", typeid(g));

  writeln("genotype_matrix[19][0]: ", genotype_matrix[19][0]);
  writeln("genotype_matrix[19][0].length: ", genotype_matrix[19][0].length);
  writeln("genotype_matrix[19][0].list.length: ", genotype_matrix[19][0].list.length);
  foreach(g; genotype_matrix[19][0].list)
    writeln("genotype_matrix[19][0].list[i]: ", g, " ", typeid(g));

  // by founders
  assert(genotype_matrix[239][2].list[0].homozygous == true);
  assert(genotype_matrix[239][2].list[0].heterozygous == false);
  assert(genotype_matrix[239][2].list[0].founders[0] == 0);
  assert(genotype_matrix[239][2].list[0].founders[1] == 0);

  assert(genotype_matrix[244][2].list[0].homozygous == false);
  assert(genotype_matrix[244][2].list[0].heterozygous == true);
  assert(genotype_matrix[244][2].list[0].founders[0] == 1);
  assert(genotype_matrix[244][2].list[0].founders[1] == 0);

  // reading phenotypes
  auto pheno_fn = to!string(buildPath(dir,"hyper_noX_phenotype.qtab"));
  writeln("reading ",pheno_fn);
  auto p_res = read_phenotype_qtab!(Phenotype!double)(pheno_fn);
  Phenotype!double[][] pheno = p_res[0];

  assert(pheno.length == 250);
  foreach(p; pheno) assert(p.length == 2);
  // 1st ind, 1st phenotype
  assert(pheno[0][0].value == 109.6);
  // 3rd ind, 1st phenotype
  assert(pheno[2][0].value == 110.1);
  // 133rd ind
  assert(pheno[132][0].value == 97);
  // all have 2nd phenotype == 1
  foreach(ph; pheno)
    assert(ph[1].value == 1);

  // Marker map reader
  auto marker_map_fn = to!string(buildPath(dir,"hyper_noX_marker_map.qtab"));
  writeln("reading ",marker_map_fn);
  auto markers = read_marker_map_qtab!(Marker)(marker_map_fn);
  assert(markers.length == 170);
  assert(markers[0].name == "D1Mit296");
  assert(markers[0].chromosome.name == "1");
  assert(markers[0].get_position == 3.3);
  assert(markers[3].name == "D1Mit178");
  assert(markers[3].chromosome.name == "1");
  assert(markers[3].get_position == 35);
  assert(markers[164].name == "D18Mit50");
  assert(markers[164].chromosome.name == "18");
  assert(markers[164].get_position == 26.2);
  assert(markers[169].name == "D19Mit137");
  assert(markers[169].chromosome.name == "19");
  assert(markers[169].get_position == 55.7);

  // marker id == numeric index
  foreach(i, m; markers) assert(m.id == i);

  // test splitting up of markers into chromosomes
  auto markers_by_chr = get_markers_by_chromosome(markers);
  writeln("markers_by_chr.length: ", markers_by_chr.length);
  writeln("markers_by_chr[0].length: ", markers_by_chr[0].length);
  writeln("markers_by_chr[0][0].name: ", markers_by_chr[0][0].name);
  writeln("markers_by_chr[0][1][0].name: ", markers_by_chr[0][1][0].name);
  foreach(chr; markers_by_chr) {
    writefln("%2s (%2d): %-9s (%3d)", chr[0].name, chr[1].length, chr[1][0].name, chr[1][0].id);
    // check that markers within chromosome are in order:
    //    contiguous ids; non-decreasing position
    assert(chr[1][0].chromosome.name == chr[0].name);
    for(auto i=1; i<chr[1].length; i++) {
      assert(chr[1][i].id == chr[1][i-1].id+1);
      assert(chr[1][i].get_position >= chr[1][i-1].get_position);
      assert(chr[1][i].chromosome.name == chr[0].name);
    }
  }

  // sorted chromosomes
  writeln("sorted chromosomes: ");
  auto markers_by_chr_sorted = sort_chromosomes_by_marker_id(markers_by_chr);
  foreach(chr; markers_by_chr_sorted) {
    writef("%2s (%2d): %-9s (%3d)", chr[0].name, chr[1].length, chr[1][0].name, chr[1][0].id);
    // check that markers within chromosome are in order:
    //    contiguous ids; non-decreasing position
    assert(chr[1][0].chromosome.name == chr[0].name);
    for(auto i=1; i<chr[1].length; i++) {
      assert(chr[1][i].id == chr[1][i-1].id+1);
      assert(chr[1][i].get_position >= chr[1][i-1].get_position);
      assert(chr[1][i].chromosome.name == chr[0].name);
    }

    auto pmap_stepped = add_stepped_markers_autosome(chr[1], 1, 0);
    auto pmap_minimal = add_minimal_markers_autosome(chr[1], 1, 0);
    writefln("\tmarkers: %3d\tpmar (stepped): %3d\tpmar (minimal): %3d", chr[1].length,
             pmap_stepped.length-chr[1].length, pmap_minimal.length-chr[1].length);
  }

  // test calc_geno_prob with hyper data
  writeln("Test calc_geno_prob with hyper data, chr 4");
  auto chr4_map = markers_by_chr_sorted[3][1];
  sort(chr4_map); // sort in place
  auto pmap_stepped_chr4 = add_stepped_markers_autosome(chr4_map, 1.0, 0.0);

  auto rec_frac = recombination_fractions(pmap_stepped_chr4, GeneticMapFunc.Carter_Falconer);
  writeln("rec_frac.length: ", rec_frac.length);
  writeln("pmap_stepped_chr4.length: ", pmap_stepped_chr4.length);
  writeln("running calc_geno_prob");
  auto chr4probs = calc_geno_prob_BC(genotype_matrix, pmap_stepped_chr4, rec_frac, 0.001);

  // unit test I'd had previously (calc genotype; no pseudomarkers)
  writeln("Test calc_geno_prob hyper hyper data, chr 5, no pseudomarkers");
  auto chr5_map = markers_by_chr_sorted[4][1];

  rec_frac = recombination_fractions(chr5_map, GeneticMapFunc.Carter_Falconer);

  double[] rec_frac_from_rqtl = [0.05499838959578960, 0.05399853076176147, 0.03299987476879847,
				 0.01099999948567994, 0.03299987476879847, 0.14181577930458458,
				 0.18529361912102324, 0.08698405702618413, 0.06599599307740488,
				 0.06499628754224140, 0.01099999948567990, 0.03299987476879847,
				 0.04399947228314834];

  foreach(i; 0..rec_frac.length)
    assert(abs(rec_frac[i] - rec_frac_from_rqtl[i]) < 1e-11);

  writeln("      - Run calcGenoprob for BC");
  auto genoprobs = calc_geno_prob_BC(genotype_matrix, chr5_map, rec_frac, 0.002);

  writeln("      - Compare results to R/qtl");
  /* probs from R/qtl for individual 1 */
  auto genoprobs_from_rqtl =  [[0.999769660491391820578, 0.0002303395086083074390],
                               [0.996571637787474373660, 0.0034283622125256579989],
                               [0.999992001840160571469, 0.0000079981598392023363],
                               [0.999999229385373444856, 0.0000007706146264155516],
                               [0.999999226845584399115, 0.0000007731544155547456],
                               [0.999988106195807358034, 0.0000118938041927743214],
                               [0.998546612868447480693, 0.0014533871315524147898],
                               [0.000838837664664903023, 0.9991611623353346960386],
                               [0.000014752882287688144, 0.9999852471177120838419],
                               [0.000009884987083730629, 0.9999901150129172355818],
                               [0.000006327226189669525, 0.9999936727738113484421],
                               [0.000385258672705765363, 0.9996147413272942205964],
                               [0.000004366470355155772, 0.9999956335296438236782],
                               [0.000092406805046532986, 0.9999075931949532591858]];

  foreach(i; 0..genoprobs[0].length)
    foreach(j; 0..1)
      assert(abs(genoprobs[0][i][j] - genoprobs_from_rqtl[i][j]) < 1e-6);


  /* probs from R/qtl for individual 88 */
  genoprobs_from_rqtl = [[0.9667756240797861222447, 0.033224375920213690405],
                         [0.0018989493588511250667, 0.998101050641148845877],
                         [0.0000061727288390725428, 0.999993827271161594972],
                         [0.0000007694148986088805, 0.999999230585101317459],
                         [0.0000007727830675603257, 0.999999227216932462525],
                         [0.0000113389799089367530, 0.999988661020090652265],
                         [0.0000754209410127622390, 0.999924579058987128555],
                         [0.0000562473769980956022, 0.999943752623000836621],
                         [0.0027022235245518625747, 0.997297776475448261024],
                         [0.9980350818191259243406, 0.001964918180874547071],
                         [0.9999911337652792608921, 0.000008866234720407318],
                         [0.9996128260665828602072, 0.000387173933416828464],
                         [0.9999956281497925925095, 0.000004371850207192675],
                         [0.9999075929709730914396, 0.000092407029026649721]];

  foreach(i; 0..genoprobs[87].length)
    foreach(j; 0..1)
      assert(abs(genoprobs[87][i][j] - genoprobs_from_rqtl[i][j]) < 1e-6);

  /* probs from R/qtl for individual 103 */
  genoprobs_from_rqtl =  [[0.0002366636922487103, 0.99976333630775027217],
                          [0.0036190176590203296, 0.99638098234097982608],
                          [0.0003819364312548427, 0.99961806356874582935],
                          [0.0549223423890561752, 0.94507765761094364443],
                          [0.0722315882263140224, 0.92776841177368585267],
                          [0.1240732453345061992, 0.87592675466549385632],
                          [0.3563901310091487917, 0.64360986899085126378],
                          [0.6573119740963657698, 0.34268802590363423022],
                          [0.7877017986121217508, 0.21229820138787824924],
                          [0.8907905734005259202, 0.10920942659947401043],
                          [0.9998456411901881502, 0.00015435880981190807],
                          [0.9994803764517402600, 0.00051962354825971387],
                          [0.9999028296560041884, 0.00009717034399581742],
                          [0.9559119082605704865, 0.04408809173942938170]];

  foreach(i; 0..genoprobs[102].length)
    foreach(j; 0..1)
      assert(abs(genoprobs[102][i][j] - genoprobs_from_rqtl[i][j]) < 1e-6);

  /* probs from R/qtl for individual 106 */
  genoprobs_from_rqtl =  [[0.00023263214083326644, 0.9997673678591675],
                          [0.00349747829854584087, 0.9965025217014539],
                          [0.00014355762173019074, 0.9998564423782705],
                          [0.01991084477676950681, 0.9800891552232307],
                          [0.02586339947739364070, 0.9741366005226065],
                          [0.04268032029532820015, 0.9573196797046718],
                          [0.09532197837358877268, 0.9046780216264112],
                          [0.09341117498753083448, 0.9065888250124690],
                          [0.06741531961276907292, 0.9325846803872311],
                          [0.03801966338345039859, 0.9619803366165495],
                          [0.00005802675751773932, 0.9999419732424817],
                          [0.00044600617494278827, 0.9995539938250571],
                          [0.00009309207758853352, 0.9999069079224117],
                          [0.04408437235616159688, 0.9559156276438385]];

  foreach(i; 0..genoprobs[105].length)
    foreach(j; 0..1)
      assert(abs(genoprobs[105][i][j] - genoprobs_from_rqtl[i][j]) < 1e-6);

  /* probs from R/qtl for individual 107 */
  genoprobs_from_rqtl =  [[0.00023263214083326644, 0.9997673678591675],
                          [0.00349747829854584087, 0.9965025217014539],
                          [0.00014355762173019074, 0.9998564423782705],
                          [0.01991084477676950681, 0.9800891552232307],
                          [0.02586339947739364070, 0.9741366005226065],
                          [0.04268032029532820015, 0.9573196797046718],
                          [0.09532197837358877268, 0.9046780216264112],
                          [0.09341117498753083448, 0.9065888250124690],
                          [0.06741531961276907292, 0.9325846803872311],
                          [0.03801966338345039859, 0.9619803366165495],
                          [0.00005802675751773932, 0.9999419732424817],
                          [0.00044600617494278827, 0.9995539938250571],
                          [0.00009309207758853352, 0.9999069079224117],
                          [0.04408437235616159688, 0.9559156276438385]];

  foreach(i; 0..genoprobs[106].length)
    foreach(j; 0..1)
      assert(abs(genoprobs[106][i][j] - genoprobs_from_rqtl[i][j]) < 1e-6);
}
