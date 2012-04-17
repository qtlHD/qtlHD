/**
 * Test scanone routines
 */

module test.scanone.test_scanone_bc;

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

  assert(lod.length == Rlod.length);
  foreach(i; 0..lod.length) {
    assert(lod[i].length == Rlod[i].length);
    foreach(j; 0..lod[i].length)
      assert(abs(lod[i][j] - Rlod[i][j]) < 1e-8);
  }

  writeln(" --Scanone for hyper chr 4, with pseudomarkers");
  /*******************************
   * Now, with pseudomarkers; need to round the map here to be like the qtab file

   data(hyper)
   hyper <- hyper["-X",]
   for(i in seq_along(hyper$geno))
     hyper$geno[[i]]$map <- round(hyper$geno[[i]]$map, 1)
   hyper$pheno[,2] <- sqrt(hyper$pheno[,1])
   hyper <- calc.genoprob(hyper, err=0.002, map="kosambi", step=1)
   out <- scanone(hyper, phe=1:2, method="hk", chr=4)

   *
   ******************************/

  // map with pseudomarkers
  auto pmap_stepped_chr4 = add_stepped_markers_autosome(chr4_map, 1.0, 0.0);

  // calc_geno_prob
  rec_frac = recombination_fractions(pmap_stepped_chr4, GeneticMapFunc.Kosambi);
  chr4probs = calc_geno_prob_BC(genotype_matrix, pmap_stepped_chr4, rec_frac, 0.002);

  // run scanone and calculate LOD scores
  rss = scanone_hk(chr4probs, pheno_rev, addcovar, intcovar, weights);
  lod = rss_to_lod(rss, rss0, pheno_rev.length);

  Rlod = [[2.6929982418539566424, 2.7693297827432559188],
          [2.9379782999527606080, 3.0172466411684411014],
          [3.1840915638100568685, 3.2660014966187418395],
          [3.4287647426149305829, 3.5129919201415304997],
          [3.6694821272770923315, 3.7556804014506042222],
          [3.9038764270228512032, 3.9916859084516715939],
          [4.1298082894884373673, 4.2188636714475649114],
          [4.3454292000375289717, 4.4353678863719494530],
          [4.5492244603103699774, 4.6396940810516866804],
          [4.7400351587605200621, 4.8307001376323626118],
          [4.9170601145632417683, 5.0076070555286946728],
          [5.0798404252197997266, 5.1699821957850815579],
          [5.2282303214495868815, 5.3177088147325264345],
          [5.3623585007695737659, 5.4509461460661441379],
          [5.4825840465937289991, 5.5700842033494097905],
          [5.5048304458567827169, 5.5920941141387459083],
          [5.9126783332657169012, 5.9913188851926122425],
          [5.7460002773761971184, 5.8037152010232944122],
          [5.4880430810297866628, 5.5360492958554345932],
          [6.2356578237854591862, 6.2811858926794457147],
          [6.5496803619669208274, 6.5900600872237760086],
          [6.5341628908615803084, 6.5739863528609703280],
          [6.5117716954952129527, 6.5509038857206576267],
          [6.6023722689192254620, 6.6432092535522428989],
          [6.7675036917740953868, 6.8123085564427015015],
          [6.8356170798882658346, 6.8837737460322614425],
          [6.8124461683376011933, 6.8629661893191098443],
          [6.7584219420734825690, 6.8050706359456967220],
          [5.8429406617573249605, 5.8509846806332745928],
          [5.8429406617572112737, 5.8509846806331040625],
          [6.1921958515827100200, 6.2028151800516866388],
          [6.2608129039480218125, 6.2735903560371184540],
          [6.2518282898538473091, 6.2647844771433369715],
          [6.3375927749257243704, 6.3536613967542621140],
          [6.3759448064398611677, 6.3953825795654495323],
          [6.3526788728440806153, 6.3752850320694278707],
          [6.3266352043648339532, 6.3504239671348159391],
          [7.4123223942889353566, 7.4178736065271380085],
          [8.0956508280114576337, 8.0833402966446783466],
          [7.6141694380023636768, 7.5959598833436530185],
          [6.3984172087770048165, 6.3758450361173117926],
          [6.3459764643110929683, 6.3205042120330858779],
          [5.1501880315685184542, 5.1249181290803846878],
          [5.1501880315685184542, 5.1249181290802425792],
          [5.1142453795847586662, 5.0903491097162145707],
          [4.8964246052400994813, 4.8767091987206185877],
          [4.8971279759203980575, 4.8792087678744167079],
          [4.7680604426780064387, 4.7586750830427604342],
          [4.7016649024417347391, 4.6908044513285744870],
          [3.7542765878220052400, 3.7312805320328550351],
          [3.7439523874517135482, 3.7192305774863712031],
          [3.7213856217363172618, 3.6949924377738057046],
          [3.6861102832016285902, 3.6581229142744575711],
          [3.6379398651754399907, 3.6084578486241127848],
          [3.5769953691816454011, 3.5461391603212177870],
          [3.5037165197560398155, 3.4716251110952498493],
          [3.4188539936086499438, 3.3856816879698214962],
          [3.3234424848496928462, 3.2893550049558371029],
          [3.2187565821748194139, 3.1839267587572805951],
          [3.1062533482900107629, 3.0708566605360090307],
          [2.9875068298714495540, 2.9517170721538832368],
          [2.8641402787436618382, 2.8281255977832415738],
          [3.0000494815576530527, 2.9620721893376469325],
          [3.1073742464776614725, 3.0677817363130657213],
          [3.1761909305197377762, 3.1354695526172804421],
          [3.1985448843691983711, 3.1572967906851374664],
          [3.1701644932360295570, 3.1290599956371067947],
          [3.0915740937781492903, 3.0512871457530081898],
          [2.9681392270772448683, 2.9292806411657750232],
          [2.8089908744919966921, 2.7720551438737004446],
          [2.6252347685751828976, 2.5905709554400004890],
          [2.4683936316783956499, 2.4356942700829620208],
          [2.4791215894343849868, 2.4467568536793748990],
          [2.5330324717062921991, 2.5024492488153100567],
          [2.5863669700854643452, 2.5577287414414513478],
          [2.6386531303006677263, 2.6121289208022346884],
          [2.6893665848368755178, 2.6651292004894457932],
          [2.7379347116727785760, 2.7161585521745053029],
          [2.7837436210181749630, 2.7646020637306492063],
          [2.8261482696294706329, 2.8098106456425853139],
          [2.8644858396842209913, 2.8511140522349478488],
          [2.8980922921955425409, 2.8878371610481110565],
          [2.9263217218849604251, 2.9193191718842967930],
          [2.9485678219716646709, 2.9449350664296503055],
          [2.9642864467497247460, 2.9641183424266444035],
          [2.9730179809346282127, 2.9763837467815221771],
          [2.9744080359622557808, 2.9813485286441050448],
          [2.9682249393014217276, 2.9787506622884336593],
          [2.9543725937726321717, 2.9684625828925845781],
          [2.9328975671498938027, 2.9504992459306151886],
          [2.9250504073116871950, 2.9436835332420230316]];

  assert(lod.length == Rlod.length);
  foreach(i; 0..lod.length) {
    assert(lod[i].length == Rlod[i].length);
    foreach(j; 0..lod[i].length)
      assert(abs(lod[i][j] - Rlod[i][j]) < 1e-8);
  }

  writeln(" --Scanone for hyper chr 15, with covariates");
  /*******************************
   * Now, with an additive covariate

   data(hyper)
   hyper <- hyper["-X",]
   for(i in seq_along(hyper$geno))
     hyper$geno[[i]]$map <- round(hyper$geno[[i]]$map, 1)
   hyper$pheno[,2] <- sqrt(hyper$pheno[,1])
   acovar <- rep(c(0,1), c(125, 125))
   hyper <- calc.genoprob(hyper, err=0.002, map="haldane", step=2.5)
   outa <- scanone(hyper, phe=1:2, method="hk", chr=15, addcovar=acovar)
   outi <- scanone(hyper, phe=1:2, method="hk", chr=15, addcovar=acovar, intcovar=acovar)

   *
   ******************************/

  // chr 15 map with pseudomarkers
  auto chr15_map = markers_by_chr_sorted[14][1];
  sort(chr15_map); // sort in place
  auto pmap_stepped_chr15 = add_stepped_markers_autosome(chr15_map, 2.5, 0.0);

  // calc_geno_prob
  rec_frac = recombination_fractions(pmap_stepped_chr15, GeneticMapFunc.Haldane);
  auto chr15probs = calc_geno_prob_BC(genotype_matrix, pmap_stepped_chr15, rec_frac, 0.002);

  // covariates
  addcovar = new double[][](genotype_matrix.length, 1);
  foreach(i; 0..addcovar.length) {
    if(i < addcovar.length/2) addcovar[i][0] = 0.0;
    else addcovar[i][0] = 1.0;
  }

  // run scanone and calculate LOD scores
  auto rss_acovar = scanone_hk(chr15probs, pheno_rev, addcovar, intcovar, weights);
  auto rss_null_acovar = scanone_hk_null(pheno_rev, addcovar, weights);
  auto lod_acovar = rss_to_lod(rss_acovar, rss_null_acovar, pheno_rev.length);

  // covariates
  intcovar = addcovar.dup;

  // run scanone and calculate LOD scores
  auto rss_icovar = scanone_hk(chr15probs, pheno_rev, addcovar, intcovar, weights);
  auto lod_icovar = rss_to_lod(rss_icovar, rss_null_acovar, pheno_rev.length);

  Rlod = [[1.05784529080335687468, 1.05912915191220236011, 1.6238828824239135429, 1.60992039473634918068],
          [1.05784529080222000630, 1.05912915191123602199, 1.6238828824223219272, 1.60992039473495651691],
          [0.73338404089736286551, 0.73585085603781408281, 1.0066633433068545855, 0.99732032855365559953],
          [0.76596193133684664645, 0.76828277987533510895, 1.0610478065971165051, 1.05097973029702984604],
          [1.03988094258761520905, 1.04050907891101473979, 1.5516607936859827532, 1.53500943006980605787],
          [1.27292169547899902682, 1.27123922131394806456, 2.0511987066356596188, 2.02781654436742542202],
          [1.28059254114498344279, 1.27881197077113029081, 2.0701381407500321075, 2.04650191029548977895],
          [1.38302113078873389895, 1.37629252662620160663, 2.4749423165510506806, 2.44218478427836771516],
          [1.35290950821831756912, 1.34464207545042313541, 2.5848178036106901345, 2.54927625114538614071],
          [1.35290950822150080057, 1.34464207545340741490, 2.5848178036366107335, 2.54927625117102252261],
          [1.68559880798045469419, 1.67222232081616084542, 2.7152846141002555669, 2.67664853871653463102],
          [1.68559880800000883028, 1.67222232083537392100, 2.7152846141370901023, 2.67664853875291441909],
          [1.68721885224840661976, 1.67292323180512880754, 2.7250207194130098287, 2.68460131058409956495],
          [1.67448792175412108918, 1.65549947849237355513, 2.7102324388212082340, 2.66137535877379605154],
          [1.61996208814252895536, 1.59641785852122097822, 2.5689095488817201840, 2.51362247673370120538],
          [1.51822428962736921676, 1.49075593487015112260, 2.3009530084401603744, 2.24267785155734600266],
          [1.37225325208839876723, 1.34198170171384845162, 1.9467440446368300400, 1.88940068898961044397],
          [1.26825705999817728298, 1.23698440248611518655, 1.7186583981031162693, 1.66349452093930949559],
          [1.27765691151842020190, 1.24621836954770515149, 1.7293194828679361308, 1.67398408583866853405],
          [1.29817100856939759979, 1.26640498476459129051, 1.7477059949850399789, 1.69225742626659325651],
          [1.31286275009802011482, 1.28093164981791574064, 1.7512959864509412000, 1.69627872050565997597],
          [1.31962884035374372615, 1.28774479348487602692, 1.7366713980828762942, 1.68273017824373027906],
          [1.31635720635574671178, 1.28478031721024876788, 1.7013730928699715150, 1.64921352397200848827],
          [1.30120588672798476182, 1.27023541238082771088, 1.6445529525307165386, 1.59488179723308576285],
          [1.27293505536920292798, 1.24289258104494138024, 1.5673995793639505791, 1.52085432863168534823],
          [1.23121702049979830917, 1.20242267478016628957, 1.4731376089270042939, 1.43022089272824359796],
          [1.17682621997278147319, 1.14957136043489072108, 1.3665467581338361924, 1.32758147721355612703],
          [1.11162798403063334263, 1.08614955505359489507, 1.2531448376711296078, 1.21825949571842784280],
          [1.03834562524400553229, 1.01480751145260228441, 1.1383070772566270534, 1.10745087445579315499],
          [1.03222795189708449470, 1.00884971841915671575, 1.1292018932165319711, 1.09866123538293436468],
          [1.28135598200265121704, 1.25629269688553790729, 1.4152035107172196149, 1.38188168574006908784],
          [1.53544029364672951488, 1.50991391355992732315, 1.6933581042912919656, 1.65932534475840043342],
          [1.73100023081690324034, 1.70667009047778606146, 1.8850774365065490201, 1.85303044310006725937],
          [1.75474700048823706311, 1.73075190445214843749, 1.9056565614816918242, 1.87414538352462045623]];


  assert(lod_acovar.length == Rlod.length);
  assert(lod_icovar.length == Rlod.length);

  foreach(i; 0..lod_acovar.length) {
    assert(lod_acovar[i].length == Rlod[i].length/2);
    assert(lod_icovar[i].length == Rlod[i].length/2);
    foreach(j; 0..lod_acovar[i].length) {
      assert(abs(lod_acovar[i][j] - Rlod[i][j]) < 1e-8);
      assert(abs(lod_icovar[i][j] - Rlod[i][j+2]) < 1e-8);
    }
  }
}
