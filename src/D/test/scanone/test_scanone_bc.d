/**
 * Test scanone routines: backcross (hyper data)
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

import qtl.core.map.genetic_map_functions;
import qtl.core.hmm.cross;
import qtl.core.hmm.calcgenoprob;
import qtl.core.scanone.scanone_hk;
import qtl.core.scanone.util;

unittest {
  writeln("Unit test " ~ __FILE__);
}

unittest {
  alias std.path.buildPath buildPath;
  auto dir = to!string(dirName(__FILE__) ~ dirSeparator ~
                       buildPath("..","..","..","..","test","data", "input", "hyper_noX_qtab"));

  // Read founder info
  auto founder_fn = to!string(buildPath(dir, "hyper_noX_founder.qtab"));
  writeln("reading ", founder_fn);
  auto info = read_founder_settings_qtab(founder_fn);
  assert(info["Cross"] == "BC");

  // First read symbol information (the GenotypeSymbolMappers)
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
  auto p_res = read_phenotype_qtab!(Phenotype)(pheno_fn);
  Phenotype[][] pheno = p_res[0];

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

  // form cross
  auto bc = form_cross("BC");

  writeln(" --Scanone for hyper chr 4, no pseudomarkers");
  // calc_geno_prob
  auto rec_frac = recombination_fractions(chr4_map, GeneticMapFunc.Carter_Falconer);
  auto chr4probs = calc_geno_prob(bc, genotype_matrix, chr4_map, rec_frac, 0.001);

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

  auto peak = get_peak_scanone(lod, chr4_map);
  foreach(i; 0..peak.length) {
    writefln("Peak for phenotype %d: max lod = %7.2f at pos = %7.2f", i,
             peak[i][0], peak[i][1].get_position);
  }

  /*******************************
   * in R:

   data(hyper)
   for(i in seq(along=hyper$geno))
     hyper$geno[[i]]$map <- round(hyper$geno[[i]]$map, 1)
   hyper <- hyper["-X",]
   hyper$pheno[,2] <- sqrt(hyper$pheno[,1])
   hyper <- calc.genoprob(hyper, err=0.001, map="c-f")
   out <- scanone(hyper, phe=1:2, method="hk", chr=4)
   paste0("auto Rlod = [", paste0(paste0("[", sprintf("%.20f", out[,3]), ", ", sprintf("%.20f", out[,4]), "]", collapse=",")), "];")

   *
   ******************************/

auto Rlod = [[2.62548658296896064712, 2.70202540976865179800],[5.44687517770546492102, 5.53476802129296174826],[5.48553679133328841999, 5.53358610883989854301],[6.55218021351367951866, 6.59266153452489334086],[6.51458532254321198707, 6.55388312935639305579],[6.83211992829330938548, 6.88385822809547676115],[5.84241390840577423660, 5.85045197132996008804],[5.84241390840554686292, 5.85045197132970429266],[6.27660932143817262840, 6.28877861215488564994],[6.30192834825606951199, 6.32617030368095356607],[8.09368283687513212499, 8.08136362100989913415],[6.38563245406589885533, 6.36301751472171872592],[5.14334162378406745120, 5.11806770484699313783],[5.14334162378395376436, 5.11806770484696471613],[4.89154205861245827691, 4.87182037362637743172],[4.76622761248972892645, 4.75686927167708972775],[3.74256346579477394698, 3.71948558894777647765],[2.79200540246358741570, 2.75701382960559726598],[2.44844415319346353499, 2.41591755884758185857],[2.93480311475536836952, 2.95402485456941121811]];

  assert(lod.length == Rlod.length);
  foreach(i; 0..lod.length) {
    assert(lod[i].length == Rlod[i].length);
    foreach(j; 0..lod[i].length)
      assert(abs(lod[i][j] - Rlod[i][j]) < 1e-10,
             to!string(i) ~ "  " ~ to!string(j) ~ "  " ~
             to!string(lod[i][j]) ~ "  " ~ to!string(Rlod[i][j]) ~
             "  " ~ to!string(log10(abs(lod[i][j] - Rlod[i][j]))));
  }

  writeln(" --Scanone for hyper chr 4, with pseudomarkers");
  // map with pseudomarkers
  auto pmap_stepped_chr4 = add_stepped_markers(chr4_map, 1.0, 0.0);

  // calc_geno_prob
  rec_frac = recombination_fractions(pmap_stepped_chr4, GeneticMapFunc.Kosambi);
  chr4probs = calc_geno_prob(bc, genotype_matrix, pmap_stepped_chr4, rec_frac, 0.002);

  // run scanone and calculate LOD scores
  rss = scanone_hk(chr4probs, pheno_rev, addcovar, intcovar, weights);
  lod = rss_to_lod(rss, rss0, pheno_rev.length);

  peak = get_peak_scanone(lod, pmap_stepped_chr4);
  foreach(i; 0..peak.length) {
    writefln("Peak for phenotype %d: max lod = %7.2f at pos = %7.2f", i,
             peak[i][0], peak[i][1].get_position);
  }

  /*******************************
   * Now, with pseudomarkers; need to round the map here to be like the qtab file

   data(hyper)
   hyper <- hyper["-X",]
   for(i in seq_along(hyper$geno))
     hyper$geno[[i]]$map <- round(hyper$geno[[i]]$map, 1)
   hyper$pheno[,2] <- sqrt(hyper$pheno[,1])
   hyper <- calc.genoprob(hyper, err=0.002, map="kosambi", step=1)
   out <- scanone(hyper, phe=1:2, method="hk", chr=4)
   paste0("Rlod = [", paste0(paste0("[", sprintf("%.20f", out[,3]), ", ", sprintf("%.20f", out[,4]), "]", collapse=",")), "];")

   *
   ******************************/

  Rlod = [[2.69299824185395664244, 2.76932978274325591883],[2.93797829995276060799, 3.01724664116844110140],[3.18409156381005686853, 3.26600149661874183948],[3.42876474261493058293, 3.51299192014153049968],[3.66948212727709233150, 3.75568040145060422219],[3.90387642702285120322, 3.99168590845167159387],[4.12980828948843736725, 4.21886367144756491143],[4.34542920003752897173, 4.43536788637194945295],[4.54922446031036997738, 4.63969408105168668044],[4.74003515876052006206, 4.83070013763236261184],[4.91706011456324176834, 5.00760705552869467283],[5.07984042521979972662, 5.16998219578508155791],[5.22823032144958688150, 5.31770881473252643445],[5.36235850076957376587, 5.45094614606614413788],[5.48258404659372899914, 5.57008420334940979046],[5.50483044585678271687, 5.59209411413874590835],[5.91267833326571690122, 5.99131888519261224246],[5.74600027737619711843, 5.80371520102329441215],[5.48804308102978666284, 5.53604929585543459325],[6.23565782378545918618, 6.28118589267944571475],[6.54968036196692082740, 6.59006008722377600861],[6.53416289086158030841, 6.57398635286097032804],[6.51177169549521295266, 6.55090388572065762673],[6.60237226891922546201, 6.64320925355224289888],[6.76750369177409538679, 6.81230855644270150151],[6.83561707988826583460, 6.88377374603226144245],[6.81244616833760119334, 6.86296618931910984429],[6.75842194207348256896, 6.80507063594569672205],[5.84294066175732496049, 5.85098468063327459276],[5.84294066175721127365, 5.85098468063310406251],[6.19219585158271002001, 6.20281518005168663876],[6.26081290394802181254, 6.27359035603711845397],[6.25182828985384730913, 6.26478447714333697149],[6.33759277492572437041, 6.35366139675426211397],[6.37594480643986116775, 6.39538257956544953231],[6.35267887284408061532, 6.37528503206942787074],[6.32663520436483395315, 6.35042396713481593906],[7.41232239428893535660, 7.41787360652713800846],[8.09565082801145763369, 8.08334029664467834664],[7.61416943800236367679, 7.59595988334365301853],[6.39841720877700481651, 6.37584503611731179262],[6.34597646431109296827, 6.32050421203308587792],[5.15018803156851845415, 5.12491812908038468777],[5.15018803156851845415, 5.12491812908024257922],[5.11424537958475866617, 5.09034910971621457065],[4.89642460524009948131, 4.87670919872061858769],[4.89712797592039805750, 4.87920876787441670785],[4.76806044267800643865, 4.75867508304276043418],[4.70166490244173473911, 4.69080445132857448698],[3.75427658782200524001, 3.73128053203285503514],[3.74395238745171354822, 3.71923057748637120312],[3.72138562173631726182, 3.69499243777380570464],[3.68611028320162859018, 3.65812291427445757108],[3.63793986517543999071, 3.60845784862411278482],[3.57699536918164540111, 3.54613916032121778699],[3.50371651975603981555, 3.47162511109524984931],[3.41885399360864994378, 3.38568168796982149615],[3.32344248484969284618, 3.28935500495583710290],[3.21875658217481941392, 3.18392675875728059509],[3.10625334829001076287, 3.07085666053600903069],[2.98750682987144955405, 2.95171707215388323675],[2.86414027874366183823, 2.82812559778324157378],[3.00004948155765305273, 2.96207218933764693247],[3.10737424647766147245, 3.06778173631306572133],[3.17619093051973777619, 3.13546955261728044206],[3.19854488436919837113, 3.15729679068513746643],[3.17016449323602955701, 3.12905999563710679467],[3.09157409377814929030, 3.05128714575300818979],[2.96813922707724486827, 2.92928064116577502318],[2.80899087449199669209, 2.77205514387370044460],[2.62523476857518289762, 2.59057095544000048903],[2.46839363167839564994, 2.43569427008296202075],[2.47912158943415761314, 2.44675685367909068191],[2.53303247170595113857, 2.50244924881491215274],[2.58636697008500959782, 2.55772874144062711821],[2.63865313030044035258, 2.61212892080183678445],[2.68936658483653445728, 2.66512920048887735902],[2.73793471167277857603, 2.71615855217424950752],[2.78374362101828864979, 2.76460206373084815823],[2.82614826962969800661, 2.80981064564278426587],[2.86448583968433467817, 2.85111405223508995732],[2.89809229219554254087, 2.88783716104842369532],[2.92632172188496042509, 2.91931917188406941932],[2.94856782197177835769, 2.94493506642962188380],[2.96428644674949737237, 2.96411834242672966866],[2.97301798093428715219, 2.97638374678106742977],[2.97440803596202840708, 2.98134852864367871916],[2.96822493930142172758, 2.97875066228849050276],[2.95437259377297323226, 2.96846258289272668662],[2.93289756714966642903, 2.95049924593010359786],[2.92505040731180088187, 2.94368353324213671840]];

  assert(lod.length == Rlod.length);
  foreach(i; 0..lod.length) {
    assert(lod[i].length == Rlod[i].length);
    foreach(j; 0..lod[i].length)
      assert(abs(lod[i][j] - Rlod[i][j]) < 1e-10,
             to!string(i) ~ "  " ~ to!string(j) ~ "  " ~
             to!string(lod[i][j]) ~ "  " ~ to!string(Rlod[i][j]) ~
             "  " ~ to!string(log10(abs(lod[i][j] - Rlod[i][j]))));
  }

  writeln(" --Scanone for hyper chr 15, with covariates");

  // chr 15 map with pseudomarkers
  auto chr15_map = markers_by_chr_sorted[14][1];
  sort(chr15_map); // sort in place
  auto pmap_stepped_chr15 = add_stepped_markers(chr15_map, 2.5, 0.0);

  // calc_geno_prob
  rec_frac = recombination_fractions(pmap_stepped_chr15, GeneticMapFunc.Haldane);
  auto chr15probs = calc_geno_prob(bc, genotype_matrix, pmap_stepped_chr15, rec_frac, 0.002);

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

  peak = get_peak_scanone(lod_acovar, pmap_stepped_chr15);
  foreach(i; 0..peak.length) {
    writefln("Peak for phenotype %d: max lod = %7.2f at pos = %7.2f", i,
             peak[i][0], peak[i][1].get_position);
  }

  writeln(" --Scanone for hyper chr 15, with interactive covariates");

  // covariates
  intcovar = addcovar.dup;

  // run scanone and calculate LOD scores
  auto rss_icovar = scanone_hk(chr15probs, pheno_rev, addcovar, intcovar, weights);
  auto lod_icovar = rss_to_lod(rss_icovar, rss_null_acovar, pheno_rev.length);

  peak = get_peak_scanone(lod_icovar, pmap_stepped_chr15);
  foreach(i; 0..peak.length) {
    writefln("Peak for phenotype %d: max lod = %7.2f at pos = %7.2f", i,
             peak[i][0], peak[i][1].get_position);
  }

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
   paste0("Rlod = [", paste0(paste0("[", sprintf("%.20f", outa[,3]), ", ", sprintf("%.20f", outa[,4]),
          ", ", sprintf("%.20f", outi[,3]), ", ", sprintf("%.20f", outi[,4]),
          "]", collapse=",")), "];")

   *
   ******************************/

  Rlod = [[1.05784529080312950100, 1.05912915191200340814, 1.62388288242379985604, 1.60992039473612180700],[1.05784529080187894579, 1.05912915191075285293, 1.62388288242209455348, 1.60992039473455861298],[0.73338404089736286551, 0.73585085603741617888, 1.00666334330696827237, 0.99732032855328611731],[0.76596193133696033328, 0.76828277987519300041, 1.06104780659711650515, 1.05097973029700142433],[1.03988094258772889589, 1.04050907891135580030, 1.55166079368621012691, 1.53500943007011869668],[1.27292169547888533998, 1.27123922131377753431, 2.05119870663554593193, 2.02781654436731173519],[1.28059254114486975595, 1.27881197077093133885, 2.07013814074980473379, 2.04650191029512029672],[1.38302113078862021212, 1.37629252662594581125, 2.47494231655093699374, 2.44218478427814034148],[1.35290950821820388228, 1.34464207545025260515, 2.58481780361069013452, 2.54927625114527245387],[1.35290950822127342690, 1.34464207545326530635, 2.58481780363661073352, 2.54927625117102252261],[1.68559880798056838103, 1.67222232081596189346, 2.71528461410025556688, 2.67664853871630725735],[1.68559880800012251711, 1.67222232083577182493, 2.71528461413720378914, 2.67664853875328390131],[1.68721885224840661976, 1.67292323180484459044, 2.72502071941300982871, 2.68460131058387219127],[1.67448792175400740234, 1.65549947849197565120, 2.71023243882109454717, 2.66137535877351183444],[1.61996208814241526852, 1.59641785852133466506, 2.56890954888160649716, 2.51362247673381489221],[1.51822428962748290360, 1.49075593486983848379, 2.30095300844027406129, 2.24267785155711862899],[1.37225325208862614090, 1.34198170171387687333, 1.94674404463705741364, 1.88940068898969570910],[1.26825705999806359614, 1.23698440248580254774, 1.71865839810300258250, 1.66349452093899685678],[1.27765691151842020190, 1.24621836954756304294, 1.72931948286793613079, 1.67398408583852642550],[1.29817100856917022611, 1.26640498476402285633, 1.74770599498481260525, 1.69225742626602482233],[1.31286275009813380166, 1.28093164981808627090, 1.75129598645094120002, 1.69627872050583050623],[1.31962884035374372615, 1.28774479348481918350, 1.73667139808298998105, 1.68273017824378712248],[1.31635720635586039862, 1.28478031721022034617, 1.70137309287008520187, 1.64921352397195164485],[1.30120588672798476182, 1.27023541238060033720, 1.64455295253060285177, 1.59488179723294365431],[1.27293505536908924114, 1.24289258104485611511, 1.56739957936383689230, 1.52085432863157166139],[1.23121702049991199601, 1.20242267477988207247, 1.47313760892700429395, 1.43022089272798780257],[1.17682621997278147319, 1.14957136043489072108, 1.36654675813394987927, 1.32758147721364139215],[1.11162798403074702946, 1.08614955505359489507, 1.25314483767135698145, 1.21825949571842784280],[1.03834562524400553229, 1.01480751145288650150, 1.13830707725674074027, 1.10745087445604895038],[1.03222795189708449470, 1.00884971841904302892, 1.12920189321653197112, 1.09866123538279225613],[1.28135598200276490388, 1.25629269688604949806, 1.41520351071733330173, 1.38188168574046699177],[1.53544029364707057539, 1.50991391356012627512, 1.69335810429151933931, 1.65932534475865622881],[1.73100023081690324034, 1.70667009047792817000, 1.88507743650666270696, 1.85303044310032305475],[1.75474700048835074995, 1.73075190445234738945, 1.90565656148191919783, 1.87414538352504678187]];

  assert(lod_acovar.length == Rlod.length);
  assert(lod_icovar.length == Rlod.length);

  foreach(i; 0..lod_acovar.length) {
    assert(lod_acovar[i].length == Rlod[i].length/2);
    assert(lod_icovar[i].length == Rlod[i].length/2);
    foreach(j; 0..lod_acovar[i].length) {
      assert(abs(lod_acovar[i][j] - Rlod[i][j]) < 1e-10,
             to!string(i) ~ "  " ~ to!string(j) ~ "  " ~
             to!string(lod_acovar[i][j]) ~ "  " ~ to!string(Rlod[i][j]) ~
             "  " ~ to!string(log10(abs(lod_acovar[i][j] - Rlod[i][j]))));
      assert(abs(lod_icovar[i][j] - Rlod[i][j+2]) < 1e-10,
             to!string(i) ~ "  " ~ to!string(j+2) ~ "  " ~
             to!string(lod_icovar[i][j]) ~ "  " ~ to!string(Rlod[i][j+2]) ~
             "  " ~ to!string(log10(abs(lod_icovar[i][j] - Rlod[i][j+2]))));
    }
  }
}
