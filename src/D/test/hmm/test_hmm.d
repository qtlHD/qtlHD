/**
 * test reading data and then apply HMM
 */

module test.hmm.test_hmm;

import std.math, std.stdio, std.path;
import std.exception, std.conv;
import std.random;

import qtl.core.primitives;
import qtl.core.phenotype, qtl.core.chromosome, qtl.core.genotype;
import qtl.plugins.qtab.read_qtab;
import qtl.core.map.make_map;


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
  assert(to!string(symbols.decode("HorA")) == "[(0,0), (0,1)]");
  assert(to!string(symbols.decode("HorB")) == "[(0,1), (1,1)]");

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
  Marker[] pmap_stepped, pmap_minimal;
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

    pmap_stepped = add_stepped_markers_autosome(chr[1], 1.0, 0.0);
    pmap_minimal = add_minimal_markers_autosome(chr[1], 1.0, 0.0);

    writefln("\tmarkers: %3d\tpmar (stepped): %3d\tpmar (minimal): %3d", chr[1].length, 
             pmap_stepped.length-chr[1].length, pmap_minimal.length-chr[1].length);
    writeln("----------");
//    foreach(m; pmap_stepped)
//      writefln("%-2s %-20s %5.1f %5s", m.chromosome.name, m.name, m.get_position, isPseudoMarker(m));
//    foreach(m; pmap_stepped)
//      writefln("%-2s %-20s %5.1f %5s", m.chromosome.name, m.name, m.get_position, isPseudoMarker(m));
    writeln("----------");
    //    foreach(m; pmap_minimal) 
    //      writefln("%-2s %-20s %5.1f %5s", m.chromosome.name, m.name, m.get_position, isPseudoMarker(m));
    writeln("----------");
  }

  auto temp = add_stepped_markers_autosome(markers_by_chr_sorted[0][1], 1.0, 0.0);
  foreach(m; temp)
    writefln("%-20s %5.1f %5s", m.name, m.get_position, isPseudoMarker(m));
  


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

  assert(genotype_matrix[19][0] == symbols.decode("H"));
  assert(genotype_matrix[19][2] == symbols.decode("A"));
  assert(genotype_matrix[19][3] == symbols.decode("NA"));
  assert(genotype_matrix[19][169] == symbols.decode("A"));

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
}



void printstuff(BC one)
{
  writeln("BC");
}

void printstuff(F2 one)
{
  writeln("F2");
}



unittest {
  auto bc = new BC;
  printstuff(bc);

  auto f2 = new F2;
  printstuff(f2);

  auto bctg = [bc.A, bc.H];
  writeln(to!string(bctg[0]));
  writeln(to!string(bctg[1]));

  auto f2tg = [f2.A, f2.H, f2.B];
  writeln(to!string(f2tg[0]));
  writeln(to!string(f2tg[1]));
  writeln(to!string(f2tg[2]));

  FounderIndex[] founder = [0,1];
  auto g00 = new TrueGenotype(founder[0],founder[0]);
  auto g01 = new TrueGenotype(founder[0],founder[1]);
  auto g10 = new TrueGenotype(founder[1],founder[0]);
  auto g11 = new TrueGenotype(founder[1],founder[1]);

  auto f2pktg = [[g00], [g01], [g10], [g11]];
  writeln(to!string(f2pktg[0]));
  writeln(to!string(f2pktg[1]));
  writeln(to!string(f2pktg[2]));
  writeln(to!string(f2pktg[3]));
}
