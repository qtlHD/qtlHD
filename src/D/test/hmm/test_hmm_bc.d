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
import qtl.core.hmm.bc;

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

  auto rec_frac = recombination_fractions(pmap_stepped_chr4, GeneticMapFunc.Kosambi);
  writeln("running calc_geno_prob");
  auto genoprobs = calc_geno_prob_BC(genotype_matrix, pmap_stepped_chr4, rec_frac, 0.001);

  /******************************
   * in R:

   data(hyper)
   hyper$geno[[4]]$map <- round(hyper$geno[[4]]$map,1) # round to 1 digit to match qtab file
   hyper <- calc.genoprob(hyper[4,], map.function="kosambi", step=1, err=0.001)
   rf <- mf.k(diff(attr(hyper$geno[["4"]]$prob, "map")))
   paste0("double[] rec_frac_stepped_rqtl = [", paste(sprintf("%.20f", rf), collapse=", "), "];")
   for(ind in c(1, 51, 128, 250))
     print(paste0("genoprobs_from_rqtl = [",
                  paste(apply(hyper$geno[["4"]]$prob[ind,,], 1, function(a) paste0("[", paste(sprintf("%.20f", a), collapse=", "), "]")), collapse=", "),
                  "];"))

   *
   ******************************/
  double[] rec_frac_stepped_rqtl = [0.00999866687996546662, 0.00999866687996546662, 0.00999866687996546662, 0.00999866687996546662, 0.00999866687996546662, 0.00999866687996546662, 0.00999866687996546662, 0.00999866687996546662, 0.00999866687996546662, 0.00999866687996546662, 0.00999866687996546662, 0.00999866687996546662, 0.00999866687996546662, 0.00999866687996546662, 0.00199998933340159255, 0.00799931740323116432, 0.00999866687996546662, 0.00399991466885112924, 0.00599971201658784734, 0.00499983333999973016, 0.00499983333999973016, 0.00599971201658784734, 0.00399991466885112924, 0.00999866687996546662, 0.00999866687996546662, 0.00899902812595466776, 0.00099999866666881411, 0.00999866687996546662, 0.00000000000000000000, 0.00999866687996546662, 0.00999866687996546662, 0.00099999866666881411, 0.00899902812595466776, 0.00999866687996546662, 0.00999866687996546662, 0.00399991466885112924, 0.00599971201658784734, 0.00499983333999973016, 0.00499983333999973016, 0.00599971201658784734, 0.00399991466885112924, 0.00699954270251874899, 0.00000000000000000000, 0.00299996400051839949, 0.00799931740323112962, 0.00199998933340162768, 0.00899902812595466776, 0.00099999866666881411, 0.00999866687996546662, 0.00999866687996546662, 0.00999866687996546662, 0.00999866687996546662, 0.00999866687996546662, 0.00999866687996546662, 0.00999866687996546662, 0.00999866687996546662, 0.00999866687996546662, 0.00999866687996546662, 0.00999866687996546662, 0.00999866687996546662, 0.00999866687996546662, 0.00999866687996546662, 0.00999866687996546662, 0.00999866687996546662, 0.00999866687996546662, 0.00999866687996546662, 0.00999866687996546662, 0.00999866687996546662, 0.00999866687996546662, 0.00999866687996546662, 0.00799931740323112962, 0.00199998933340162768, 0.00999866687996546662, 0.00999866687996546662, 0.00999866687996546662, 0.00999866687996546662, 0.00999866687996546662, 0.00999866687996546662, 0.00999866687996546662, 0.00999866687996546662, 0.00999866687996546662, 0.00999866687996546662, 0.00999866687996546662, 0.00999866687996546662, 0.00999866687996546662, 0.00999866687996546662, 0.00999866687996546662, 0.00999866687996546662, 0.00999866687996546662, 0.00299996400051836393];
  assert(rec_frac.length == rec_frac_stepped_rqtl.length);
  foreach(i; 0..rec_frac.length)
    assert(abs(rec_frac[i] - rec_frac_stepped_rqtl[i]) < 1e-12);

  auto genoprobs_from_rqtl = [[0.00697533053351557814, 0.99302466946648448864], [0.07724982761110457674, 0.92275017238889545101], [0.14735182007044453356, 0.85264817992955543868], [0.21730991326343712666, 0.78269008673656259578], [0.28715265382340843248, 0.71284734617659151201], [0.35690854131365717938, 0.64309145868634265408], [0.42660603985679840244, 0.57339396014320154205], [0.49627358974964180804, 0.50372641025035824747], [0.56593961906835033027, 0.43406038093164983627], [0.63563255526861117684, 0.36436744473138854561], [0.70538083678555585365, 0.29461916321444425737], [0.77521292463815483220, 0.22478707536184511229], [0.84515731404283567230, 0.15484268595716449424], [0.91524254604104682276, 0.08475745395895320500], [0.98549721914552856905, 0.01450278085447143789], [0.99945950995821275509, 0.00054049004178735886], [0.99895627438052747582, 0.00104372561947255302], [0.99850904692273645402, 0.00149095307726353253], [0.99838774865812407455, 0.00161225134187583871], [0.99826591547874155008, 0.00173408452125860453], [0.99821987860280791960, 0.00178012139719209840], [0.99822416359536858543, 0.00177583640463130138], [0.99829583871570393150, 0.00170416128429608542], [0.99838375778060706178, 0.00161624221939270642], [0.99874699806455147222, 0.00125300193544856074], [0.99931375372500708121, 0.00068624627499292562], [0.99999726983824788196, 0.00000273016175201058], [0.99996600006829816643, 0.00003399993170180258], [0.99976258634980175177, 0.00023741365019843994], [0.99976258634980130768, 0.00023741365019853990], [0.99976310242200450151, 0.00023689757799542162], [0.99996754849549240873, 0.00003245150450776578], [0.99999892060861028664, 0.00000107939138966911], [0.99977920589304214438, 0.00022079410695784849], [0.99972848283539106351, 0.00027151716460908860], [0.99988167565238650347, 0.00011832434761360928], [0.99999926086128398950, 0.00000073913871575836], [0.99990211788414073979, 0.00009788211585909965], [0.99987681745866296978, 0.00012318254133700696], [0.99990200625802716505, 0.00009799374197298069], [0.99999901513465161429, 0.00000098486534820288], [0.99983685782696318789, 0.00016314217303687207], [0.99963020160247595403, 0.00036979839752410647], [0.99963020160247550994, 0.00036979839752433047], [0.99957223076264778427, 0.00042776923735215799], [0.99950608408843755193, 0.00049391591156241034], [0.99950980254985832829, 0.00049019745014162650], [0.99962722074439913644, 0.00037277925560095230], [0.99965024738314223196, 0.00034975261685764916], [0.99999464527024717331, 0.00000535472975279562], [0.99887788517641473707, 0.00112211482358521479], [0.99796469386801645030, 0.00203530613198347641], [0.99725469871429250901, 0.00274530128570762500], [0.99674760999935119177, 0.00325239000064885811], [0.99644322080395109076, 0.00355677919604909486], [0.99634140692106498705, 0.00365859307893504071], [0.99644212680519872460, 0.00355787319480128233], [0.99674542155543688349, 0.00325457844456307054], [0.99725141493221425293, 0.00274858506778570934], [0.99796031340781643504, 0.00203968659218366471], [0.99887240625063100641, 0.00112759374936890447], [0.99998806564318609702, 0.00001193435681393077], [0.99909319503552018737, 0.00090680496447962007], [0.99840198107119537685, 0.00159801892880481614], [0.99791414169804615319, 0.00208585830195373101], [0.99762947785158995195, 0.00237052214841013869], [0.99754787337379358103, 0.00245212662620630014], [0.99766929496567824387, 0.00233070503432163210], [0.99799379217372941042, 0.00200620782627047725], [0.99852149741011553274, 0.00147850258988447810], [0.99925262600671882485, 0.00074737399328136088], [0.99998261405276656077, 0.00001738594723355960], [0.99963439705369561139, 0.00036560294630448261], [0.99800146335152850607, 0.00199853664847150043], [0.99657174080794586057, 0.00342825919205398070], [0.99534464601989447452, 0.00465535398010576574], [0.99431967826725098369, 0.00568032173274902324], [0.99349641930850529636, 0.00650358069149478864], [0.99287453321009222496, 0.00712546678990779066], [0.99245376620931491374, 0.00754623379068525106], [0.99223394661079433732, 0.00776605338920562886], [0.99221498471641067507, 0.00778501528358946544], [0.99239687278869903597, 0.00760312721130105857], [0.99277968504769442681, 0.00722031495230566513], [0.99336357770121608279, 0.00663642229878382527], [0.99414878900860981226, 0.00585121099139027621], [0.99513563937796856163, 0.00486436062203136205], [0.99632453149687838501, 0.00367546850312144933], [0.99771595049673478250, 0.00228404950326498374], [0.99931046415070301503, 0.00068953584929692393], [0.99982486353059174533, 0.00017513646940798419]];
  auto ind=0;
  foreach(i; 0..genoprobs.length)
    foreach(j; 0..1)
      assert(abs(genoprobs[i][ind][j] - genoprobs_from_rqtl[i][j]) < 1e-12,
             to!string(i) ~ "  " ~
             to!string(genoprobs[i][ind][j]) ~ "  " ~
             to!string(genoprobs_from_rqtl[i][j]) ~ "  " ~
             to!string(abs(genoprobs[i][ind][j] - genoprobs_from_rqtl[i][j])));

  genoprobs_from_rqtl = [[0.99302466946648404456, 0.00697533053351557207], [0.92275017238889578408, 0.07724982761110457674], [0.85264817992955577175, 0.14735182007044453356], [0.78269008673656337294, 0.21730991326343712666], [0.71284734617659184508, 0.28715265382340793288], [0.64309145868634287613, 0.35690854131365684632], [0.57339396014320154205, 0.42660603985679840244], [0.50372641025035824747, 0.49627358974964136396], [0.43406038093164983627, 0.56593961906834966413], [0.36436744473138887868, 0.63563255526861117684], [0.29461916321444442390, 0.70538083678555518752], [0.22478707536184561189, 0.77521292463815483220], [0.15484268595716482730, 0.84515731404283500616], [0.08475745395895328826, 0.91524254604104682276], [0.01450278085447148819, 0.98549721914552901314], [0.00054049004178735984, 0.99945950995821275509], [0.00104372561947255389, 0.99895627438052791991], [0.00149095307726353123, 0.99850904692273645402], [0.00161225134187584153, 0.99838774865812451864], [0.00173408452125860908, 0.99826591547874155008], [0.00178012139719210013, 0.99821987860280836369], [0.00177583640463130289, 0.99822416359536902952], [0.00170416128429608542, 0.99829583871570393150], [0.00161624221939270923, 0.99838375778060706178], [0.00125300193544855966, 0.99874699806455147222], [0.00068624627499292378, 0.99931375372500663712], [0.00000273016175201058, 0.99999726983824832605], [0.00003399993170180258, 0.99996600006829816643], [0.00023741365019843994, 0.99976258634980219586], [0.00023741365019854031, 0.99976258634980219586], [0.00023689757799542202, 0.99976310242200494560], [0.00003245150450776595, 0.99996754849549196464], [0.00000107939138966929, 0.99999892060861117482], [0.00022079410695789714, 0.99977920589304125620], [0.00027151716460919182, 0.99972848283539106351], [0.00011832434761376713, 0.99988167565238561529], [0.00000073913871593786, 0.99999926086128310132], [0.00009788211589176772, 0.99990211788410832128], [0.00012318254139672473, 0.99987681745860390592], [0.00009799374205975481, 0.99990200625794056766], [0.00000098486546748835, 0.99999901513453259838], [0.00016314219476064364, 0.99983685780523878783], [0.00036979845717293011, 0.99963020154282711260], [0.00036979845717315346, 0.99963020154282711260], [0.00042776931319318871, 0.99957223068680722822], [0.00049391603081305800, 0.99950608396918660947], [0.00049019758018445814, 0.99950980241981512986], [0.00037277943457660051, 0.99962722056542374638], [0.00034975280123049118, 0.99965024719876938164], [0.00000535496862781535, 0.99999464503137236449], [0.00112216945205138993, 0.99887783054794831994], [0.00203541517233213733, 0.99796458482766758902], [0.00274546478243304064, 0.99725453521756657338], [0.00325260802046641083, 0.99674739197953365899], [0.00355705182792250864, 0.99644294817207768045], [0.00365892043411265811, 0.99634107956588680022], [0.00355825540686147915, 0.99644174459313861192], [0.00325501566946875820, 0.99674498433053149160], [0.00274907748394795154, 0.99725092251605207405], [0.00204023440053452323, 0.99795976559946519835], [0.00112819717344341116, 0.99887180282655629870], [0.00001259364284129496, 0.99998740635715899217], [0.00096260219631365382, 0.99903739780368627788], [0.00170897687469211807, 0.99829102312530815233], [0.00225202223865431312, 0.99774797776134560578], [0.00259195987967273658, 0.99740804012032702186], [0.00272892851043592642, 0.99727107148956395388], [0.00266298402145063009, 0.99733701597854917953], [0.00239409950384811927, 0.99760590049615216479], [0.00192216523840389969, 0.99807783476159628311], [0.00124698865076647183, 0.99875301134923333279], [0.00056181885170347647, 0.99943818114829685584], [0.01194625408226506792, 0.98805374591773509341], [0.06920727549383844102, 0.93079272450616190593], [0.12629251049693926423, 0.87370748950306054148], [0.18322525291216171239, 0.81677474708783825985], [0.24002873433499127831, 0.75997126566500872169], [0.29672613361553373812, 0.70327386638446631739], [0.35334058631672299322, 0.64665941368327717331], [0.40989519415486364018, 0.59010480584513613778], [0.46641303442636644627, 0.53358696557363338719], [0.52291716942451194061, 0.47708283057548789285], [0.57943065585009845275, 0.42056934414990115867], [0.63597655421980658907, 0.36402344578019296684], [0.69257793827611757731, 0.30742206172388214513], [0.74925790440263306813, 0.25074209559736704289], [0.80603958104862916123, 0.19396041895137039468], [0.86294613816670162620, 0.13705386183329842931], [0.92000079666733680650, 0.07999920333266360983], [0.97722683789428510526, 0.02277316210571537700], [0.99431226749532686249, 0.00568773250467338471]];
  ind = 50;
  foreach(i; 0..genoprobs.length)
    foreach(j; 0..1)
      assert(abs(genoprobs[i][ind][j] - genoprobs_from_rqtl[i][j]) < 1e-12,
             to!string(i) ~ "  " ~
             to!string(genoprobs[i][ind][j]) ~ "  " ~
             to!string(genoprobs_from_rqtl[i][j]) ~ "  " ~
             to!string(abs(genoprobs[i][ind][j] - genoprobs_from_rqtl[i][j])));

  genoprobs_from_rqtl = [[0.12473022081297534258, 0.87526977918702453252], [0.11707269570316086726, 0.88292730429683941029], [0.10925891582928166923, 0.89074108417071828914], [0.10128569275239159009, 0.89871430724760836828], [0.09314977297221328778, 0.90685022702778694814], [0.08484783659953677726, 0.91515216340046312560], [0.07637649600152791873, 0.92362350399847203963], [0.06773229441939299100, 0.93226770558060712002], [0.05891170455783703575, 0.94108829544216288099], [0.04991112714573924664, 0.95008887285426091296], [0.04072688946745838673, 0.95927311053254160633], [0.03135524386416869286, 0.96864475613583156388], [0.02179236620461550561, 0.97820763379538477889], [0.01203435432466561789, 0.98796564567533429191], [0.00207722643501609969, 0.99792277356498426677], [0.00007754728911664333, 0.99992245271088187053], [0.00062919898963974589, 0.99937080101035991042], [0.00113691605524871334, 0.99886308394475109562], [0.00128222492314158205, 0.99871777507685854935], [0.00144010618860475380, 0.99855989381139542704], [0.00151612037371586785, 0.99848387962628415426], [0.00154178602534930714, 0.99845821397465073233], [0.00150605645325654035, 0.99849394354674347873], [0.00144203712976029947, 0.99855796287024001234], [0.00113885699603418389, 0.99886114300396600107], [0.00063211490857323501, 0.99936788509142704395], [0.00000253798221195826, 0.99999746201778749732], [0.00003381370934738953, 0.99996618629065181505], [0.00023728749755948180, 0.99976271250244053146], [0.00023728749755958211, 0.99976271250244053146], [0.00023683144359491371, 0.99976316855640490200], [0.00003244536135934612, 0.99996755463863984392], [0.00000107919337595059, 0.99999892080662344807], [0.00022079396296209721, 0.99977920603703818170], [0.00027151708063776321, 0.99972848291936144971], [0.00011832432363244565, 0.99988167567636709876], [0.00000073913858403861, 0.99999926086141655013], [0.00009788211576813127, 0.99990211788423266626], [0.00012318254127995350, 0.99987681745871936911], [0.00009799374194983661, 0.99990200625805059076], [0.00000098486536578770, 0.99999901513463407277], [0.00016314217628721450, 0.99983685782371345407], [0.00036979840644931004, 0.99963020159355087113], [0.00036979840644953469, 0.99963020159355087113], [0.00042776924870026068, 0.99957223075129952861], [0.00049391592940604989, 0.99950608407059404747], [0.00049019746960013932, 0.99950980253039944934], [0.00037277928238145958, 0.99962722071761855869], [0.00034975264444575809, 0.99965024735555430002], [0.00000535476549626150, 0.99999464523450320996], [0.00112212299777922288, 0.99887787700222041565], [0.00203532244796353476, 0.99796467755203666083], [0.00274532575013152142, 0.99725467424986824483], [0.00325242262349938888, 0.99674757737650043765], [0.00355681999063811377, 0.99644318000936171753], [0.00365864206190891954, 0.99634135793809086579], [0.00355793038614768842, 0.99644206961385239918], [0.00325464386761912432, 0.99674535613238046672], [0.00274865874924750146, 0.99725134125075287628], [0.00203976856211716094, 0.99796023143788248344], [0.00112768404122220826, 0.99887231595877812307], [0.00001203300743094004, 0.99998796699256797815], [0.00091515404381600914, 0.99908484595618429758], [0.00161462184373025938, 0.99838537815626937633], [0.00211072182734306111, 0.99788927817265715703], [0.00240365642990783644, 0.99759634357009230321], [0.00249354518436647457, 0.99750645481563304795], [0.00238042477012512851, 0.99761957522987454450], [0.00206424902802136495, 0.99793575097197884105], [0.00154488894148876489, 0.99845511105851092459], [0.00082213258391131085, 0.99917786741608871637], [0.00009885080424079000, 0.99990114919575823205], [0.00209844473653425859, 0.99790155526346591142], [0.01205514831672647828, 0.98794485168327339508], [0.02181274437227692653, 0.97818725562772324000], [0.03137521452280987722, 0.96862478547718999788], [0.04074646076617333829, 0.95925353923382672416], [0.04993030707066171969, 0.95006969292933829419], [0.05893050093539932421, 0.94106949906460057864], [0.06775071491951970004, 0.93224928508048021669], [0.07639454814076565570, 0.92360545185923459410], [0.08486552774412109423, 0.91513447225587896128], [0.09316711034107470757, 0.90683288965892538958], [0.10130268342010141647, 0.89869731657989848639], [0.10927556672893845580, 0.89072443327106165523], [0.11708901362921977918, 0.88291098637078035960], [0.12474621242402060550, 0.87525378757597938062], [0.13225028765885568638, 0.86774971234114417484], [0.13960430139666013538, 0.86039569860333986462], [0.14681125446727461004, 0.85318874553272538996], [0.14893036151124747213, 0.85106963848875250012]];
  ind = 127;
  foreach(i; 0..genoprobs.length)
    foreach(j; 0..1)
      assert(abs(genoprobs[i][ind][j] - genoprobs_from_rqtl[i][j]) < 1e-12,
             to!string(i) ~ "  " ~
             to!string(genoprobs[i][ind][j]) ~ "  " ~
             to!string(genoprobs_from_rqtl[i][j]) ~ "  " ~
             to!string(abs(genoprobs[i][ind][j] - genoprobs_from_rqtl[i][j])));

  genoprobs_from_rqtl = [[0.12479706555589180350, 0.87520293444410846018], [0.11714090443893196425, 0.88285909556106856311], [0.10932851639070967353, 0.89067148360929071504], [0.10135671354021769108, 0.89864328645978286403], [0.09322224296670689192, 0.90677775703329310808], [0.08492178537232104862, 0.91507821462767890974], [0.07645195372764611308, 0.92354804627235376202], [0.06780929188962203646, 0.93219070811037807456], [0.05899027319125160762, 0.94100972680874894749], [0.04999129900253231756, 0.95000870099746792530], [0.04080869726202351888, 0.95919130273797670316], [0.03143872097844949143, 0.96856127902155031428], [0.02187754670172688612, 0.97812245329827307572], [0.01212127296279367529, 0.98787872703720636114], [0.00216591868160103361, 0.99783408131839890221], [0.00016659572755563801, 0.99983340427244404847], [0.00143565910630060649, 0.99856434089369916496], [0.00284223784056989977, 0.99715776215942975025], [0.00334512857125040305, 0.99665487142874964466], [0.00404070865603974906, 0.99595929134396021798], [0.00456464556320258091, 0.99543535443679731500], [0.00503854184809374855, 0.99496145815190584205], [0.00554169514066817949, 0.99445830485933173204], [0.00583653784626064417, 0.99416346215373951889], [0.00643717849951201161, 0.99356282150048813584], [0.00683641919679464898, 0.99316358080320565893], [0.00702357041096137325, 0.99297642958903831190], [0.09656485115411657383, 0.90343514884588282943], [0.99999809700898423248, 0.00000190299101614499], [0.99999809700898512066, 0.00000190299101524978], [0.99988656679323273391, 0.00011343320676696375], [0.99997901695881419304, 0.00002098304118555143], [0.99999929027443601459, 0.00000070972556387527], [0.99977947464025462843, 0.00022052535974542448], [0.99972863944236811662, 0.00027136055763190572], [0.99988172018303234090, 0.00011827981696766794], [0.99999926083529067089, 0.00000073916470965824], [0.99990206855445629142, 0.00009793144554389656], [0.99987672707823060936, 0.00012327292176967207], [0.99990187481771741762, 0.00009812518228223880], [0.99999883435554781652, 0.00000116564445206409], [0.99980393498550768072, 0.00019606501449207062], [0.99953980252984875321, 0.00046019747015127839], [0.99953980252984875321, 0.00046019747015158170], [0.99945729205303945886, 0.00054270794696088516], [0.99932535717043491363, 0.00067464282956538962], [0.99931271984572522804, 0.00068728015427489957], [0.99935597930379982401, 0.00064402069620015680], [0.99937082637918372985, 0.00062917362081614016], [0.99963262504453931356, 0.00036737495546082034], [0.91608727072740803177, 0.08391272927259167680], [0.83271170220927426264, 0.16728829779072554307], [0.74947189781124134988, 0.25052810218875887216], [0.66633389125341380144, 0.33366610874658636510], [0.58326375779484518347, 0.41673624220515476102], [0.50022760039043490110, 0.49977239960956482134], [0.41719153585912377435, 0.58280846414087617013], [0.33412168105774886495, 0.66587831894225124607], [0.25098413905491540055, 0.74901586094508509905], [0.16774498529924031232, 0.83225501470075979871], [0.08437025377632698109, 0.91562974622367299116], [0.00082592314881693046, 0.99917407685118297955], [0.00164588014708512767, 0.99835411985291477865], [0.00226248208410684357, 0.99773751791589293525], [0.00267598056635962034, 0.99732401943364079600], [0.00288654432327868291, 0.99711345567672149404], [0.00289425927610753499, 0.99710574072389224209], [0.00269912857295842796, 0.99730087142704149095], [0.00230107259009692611, 0.99769892740990262503], [0.00169992889945120543, 0.99830007110054874708], [0.00089545220233272916, 0.99910454779766699751], [0.00010694667580652065, 0.99989305332419387007], [0.00210650822478643729, 0.99789349177521335932], [0.01206305055671282392, 0.98793694944328680485], [0.02182048858853275278, 0.97817951141146741723], [0.03138280387538858185, 0.96861719612461139040], [0.04075389835193536647, 0.95924610164806434209], [0.04993759592453880097, 0.95006240407546127535], [0.05893764403163266608, 0.94106235596836762536], [0.06775771517287348944, 0.93224228482712656607], [0.07640140840771668385, 0.92359859159228363534], [0.08487225082402412302, 0.91512774917597572433], [0.09317369897730497230, 0.90682630102269523587], [0.10130914030117378621, 0.89869085969882589460], [0.10928189448960487495, 0.89071810551039498627], [0.11709521485154399112, 0.88290478514845582847], [0.12475228963843221219, 0.87524771036156789883], [0.13225624334518221170, 0.86774375665481784381], [0.13961013798513935424, 0.86038986201486089556], [0.14681697433954588039, 0.85318302566045456370], [0.14893604706469706911, 0.85106395293530312518]];
  ind = 249;
  foreach(i; 0..genoprobs.length)
    foreach(j; 0..1)
      assert(abs(genoprobs[i][ind][j] - genoprobs_from_rqtl[i][j]) < 1e-12,
             to!string(i) ~ "  " ~
             to!string(genoprobs[i][ind][j]) ~ "  " ~
             to!string(genoprobs_from_rqtl[i][j]) ~ "  " ~
             to!string(abs(genoprobs[i][ind][j] - genoprobs_from_rqtl[i][j])));

  // unit test I'd had previously (calc genotype; no pseudomarkers)
  writeln("Test calc_geno_prob hyper hyper data, chr 5, no pseudomarkers");
  auto chr5_map = markers_by_chr_sorted[4][1];

  rec_frac = recombination_fractions(chr5_map, GeneticMapFunc.Carter_Falconer);
  /******************************
   * in R:

   data(hyper)
   m <- pull.map(hyper)[[5]]
   rf <- mf.cf(diff(m))
   paste0("double[] rec_frac_from_rqtl = [", paste(sprintf("%.20f", rf), collapse=", "), "];")

   *
   ******************************/

  double[] rec_frac_from_rqtl = [0.05499838959578959990, 0.05399853076176146932, 0.03299987476879847476, 0.01099999948567993539, 0.03299987476879847476, 0.14181577930458447168, 0.18529361912102329524, 0.08698405702618412783, 0.06599599307740487697, 0.06499628754224140437, 0.01099999948567990243, 0.03299987476879835679, 0.04399947228314848163];

  foreach(i; 0..rec_frac.length) {
    writefln("%.8f", log(abs(rec_frac[i] - rec_frac_from_rqtl[i])));
    assert(abs(rec_frac[i] - rec_frac_from_rqtl[i]) < 1e-11);
  }

  writeln("      - Run calcGenoprob for BC");
  genoprobs = calc_geno_prob_BC(genotype_matrix, chr5_map, rec_frac, 0.002);

  writeln("      - Compare results to R/qtl");
  /* probs from R/qtl for individual 1 */
  /******************************
   * in R:

   data(hyper)
   hyper$geno[["5"]]$map <- round(hyper$geno[["5"]]$map, 1)
   hyper <- calc.genoprob(hyper[5,], step=0, err=0.002, map.function="c-f")
   for(ind in c(1, 88, 103, 106, 107))
     print(paste0("genoprobs_from_rqtl = [",
                  paste(apply(hyper$geno[["5"]]$prob[ind,,], 1, function(a) paste0("[", paste(sprintf("%.20f", a), collapse=", "), "]")), collapse=", "),
                  "];"))

   *
   ******************************/
  genoprobs_from_rqtl = [[0.99976966049139626147, 0.00023033950860386717], [0.99657163778760493589, 0.00342836221239520116], [0.99999200184016101556, 0.00000799815983879991], [0.99999922938537344486, 0.00000077061462632092], [0.99999922684558439911, 0.00000077315441545992], [0.99998810619580780212, 0.00001189380419231287], [0.99854661286845014523, 0.00145338713154999572], [0.00083883766465980478, 0.99916116233534046920], [0.00001475288228728719, 0.99998524711771386020], [0.00000988498708340747, 0.99999011501291723558], [0.00000632722618927354, 0.99999367277381046026], [0.00038525867265866486, 0.99961474132734129405], [0.00000436647035485508, 0.99999563352964471186], [0.00009240680504433276, 0.99990759319495681190]];

  ind = 0;
  foreach(i; 0..genoprobs.length)
    foreach(j; 0..1)
      assert(abs(genoprobs[i][ind][j] - genoprobs_from_rqtl[i][j]) < 1e-12,
             to!string(i) ~ "  " ~
             to!string(genoprobs[i][ind][j]) ~ "  " ~
             to!string(genoprobs_from_rqtl[i][j]) ~ "  " ~
             to!string(abs(genoprobs[i][ind][j] - genoprobs_from_rqtl[i][j])));


  /* probs from R/qtl for individual 88 */
  genoprobs_from_rqtl = [[0.96677562407916872722, 0.03322437592083161279], [0.00189894935884924961, 0.99810105064115062223], [0.00000617272883884484, 0.99999382727116159497], [0.00000076941489851444, 0.99999923058510131746], [0.00000077278306746554, 0.99999922721693246253], [0.00001133897990848853, 0.99998866102009065226], [0.00007542094101165588, 0.99992457905898890491], [0.00005624737699739280, 0.99994375262300083662], [0.00270222352456174833, 0.99729777647543860208], [0.99803508181912592434, 0.00196491818087399239], [0.99999113376528014907, 0.00000886623471993617], [0.99961282606662993366, 0.00038717393336970112], [0.99999562814979348069, 0.00000437185020689194], [0.99990759297097575597, 0.00009240702902444948]];

  ind = 87;
  foreach(i; 0..genoprobs.length)
    foreach(j; 0..1)
      assert(abs(genoprobs[i][ind][j] - genoprobs_from_rqtl[i][j]) < 1e-12,
             to!string(i) ~ "  " ~
             to!string(genoprobs[i][ind][j]) ~ "  " ~
             to!string(genoprobs_from_rqtl[i][j]) ~ "  " ~
             to!string(abs(genoprobs[i][ind][j] - genoprobs_from_rqtl[i][j])));

  /* probs from R/qtl for individual 103 */
  genoprobs_from_rqtl = [[0.00023666369224434103, 0.99976333630775560124], [0.00361901765888844899, 0.99638098234111172058], [0.00038193643125135275, 0.99961806356874938206], [0.05492234238792238765, 0.94507765761207762623], [0.07223158822379814603, 0.92776841177620206214], [0.12407324533105812403, 0.87592675466894165393], [0.35639013100723504479, 0.64360986899276495521], [0.65731197409711217272, 0.34268802590288777177], [0.78770179861309896907, 0.21229820138690097542], [0.89079057340114586871, 0.10920942659885408965], [0.99984564119019325723, 0.00015435880980659578], [0.99948037645179133026, 0.00051962354820860393], [0.99990282965600840726, 0.00009717034399169741], [0.95591190826157412808, 0.04408809173842598989]];

  ind = 102;
  foreach(i; 0..genoprobs.length)
    foreach(j; 0..1)
      assert(abs(genoprobs[i][ind][j] - genoprobs_from_rqtl[i][j]) < 1e-12,
             to!string(i) ~ "  " ~
             to!string(genoprobs[i][ind][j]) ~ "  " ~
             to!string(genoprobs_from_rqtl[i][j]) ~ "  " ~
             to!string(abs(genoprobs[i][ind][j] - genoprobs_from_rqtl[i][j])));

  /* probs from R/qtl for individual 106 */
  genoprobs_from_rqtl = [[0.00023263214082880717, 0.99976736785917019734], [0.00349747829841354002, 0.99650252170158670761], [0.00014355762172606518, 0.99985644237827497882], [0.01991084477597598185, 0.98008915522402406673], [0.02586339947601953154, 0.97413660052398043376], [0.04268032029338705846, 0.95731967970661291378], [0.09532197837124513351, 0.90467802162875465832], [0.09341117498538110964, 0.90658882501461890424], [0.06741531961110046323, 0.93258468038889963392], [0.03801966338247741301, 0.96198033661752269108], [0.00005802675751460497, 0.99994197324248523806], [0.00044600617489225556, 0.99955399382510767836], [0.00009309207758432434, 0.99990690792241432305], [0.04408437235515812180, 0.95591562764484172554]];

  ind = 105;
  foreach(i; 0..genoprobs.length)
    foreach(j; 0..1)
      assert(abs(genoprobs[i][ind][j] - genoprobs_from_rqtl[i][j]) < 1e-12,
             to!string(i) ~ "  " ~
             to!string(genoprobs[i][ind][j]) ~ "  " ~
             to!string(genoprobs_from_rqtl[i][j]) ~ "  " ~
             to!string(abs(genoprobs[i][ind][j] - genoprobs_from_rqtl[i][j])));

  /* probs from R/qtl for individual 107 */
  genoprobs_from_rqtl = [[0.00023263214082880717, 0.99976736785917019734], [0.00349747829841354002, 0.99650252170158670761], [0.00014355762172606518, 0.99985644237827497882], [0.01991084477597598185, 0.98008915522402406673], [0.02586339947601953154, 0.97413660052398043376], [0.04268032029338705846, 0.95731967970661291378], [0.09532197837124513351, 0.90467802162875465832], [0.09341117498538110964, 0.90658882501461890424], [0.06741531961110046323, 0.93258468038889963392], [0.03801966338247741301, 0.96198033661752269108], [0.00005802675751460497, 0.99994197324248523806], [0.00044600617489225556, 0.99955399382510767836], [0.00009309207758432434, 0.99990690792241432305], [0.04408437235515812180, 0.95591562764484172554]];

  ind = 106;
  foreach(i; 0..genoprobs.length)
    foreach(j; 0..1)
      assert(abs(genoprobs[i][ind][j] - genoprobs_from_rqtl[i][j]) < 1e-12,
             to!string(i) ~ "  " ~
             to!string(genoprobs[i][ind][j]) ~ "  " ~
             to!string(genoprobs_from_rqtl[i][j]) ~ "  " ~
             to!string(abs(genoprobs[i][ind][j] - genoprobs_from_rqtl[i][j])));
}
