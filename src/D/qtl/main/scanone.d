// Entry point for command line (CLI) scanone tool

import std.getopt;
import std.stdio;
import std.conv;
import std.exception;
import std.file;
import std.string;

import qtl.plugins.qtab.read_qtab;

static string ver = import("VERSION");

string copyright = "; qtlHD project (c) 2012";
string usage = "
  usage: scanone [options] inputfile(s)

  options:

    -v --verbosity    Set verbosity level (default 1)
    -d --debug        Set debug message level (default 0)

  examples:

    Execute scanone with the listeria dataset

      ./scanone -v 1 -d 3 ../../test/data/input/listeria_qtab/listeria_symbol.qtab ../../test/data/input/listeria_qtab/listeria_founder.qtab ../../test/data/input/listeria_qtab/listeria_marker_map.qtab ../../test/data/input/listeria_qtab/listeria_genotype.qtab ../../test/data/input/listeria_qtab/listeria_phenotype.qtab
";

int main(string[] args) {
  writeln("scanone ",strip(ver)," ",copyright);
  if (args.length == 1) {
    writeln(usage);
    return 0;
  }
  writeln(args);
  uint verbosity = 1;
  uint debug_level = 0;
  getopt(args, "v|verbose", (string o, string v) { verbosity = to!int(v); },
               "d|debug", (string o, string d) { debug_level = to!int(d); }
  );

  writeln("Verbosity: ",verbosity);
  writeln("Debug level: ",debug_level);
  auto res = load_qtab(args[1..$]);
  auto s = res[0];
  auto i = res[3];
  auto p = res[4];
  auto o = res[5];
  auto g = res[6];

  if (debug_level > 2) {
    writeln("* Symbol data");
    writeln(s);
    writeln(o);
    writeln("* Individuals");
    writeln(i);
    writeln("* Genotype data");
    writeln(g);
    writeln("* Phenotype data");
    writeln(p);
  }
  return 0;
}
