// Entry point for CLI scanone tool

import std.getopt;
import std.stdio;
import std.conv;
import std.exception;
import std.file;
import std.string;

static string ver = import("VERSION");

string copyright = "; qtlHD project (c) 2012";
string usage = "
  usage: scanone [options] inputfile(s)

  options:

    -v --verbosity    Set verbosity level (default 1)

  examples:

    scanone listeria_symbol.qtab listeria_founder.qtab listeria_marker_map.qtab listeria_genotype.qtab listeria_phenotype.qtab 
";

int main(string[] args) {
  writeln("scanone ",strip(ver)," ",copyright);
  if (args.length == 1) {
    writeln(usage);
    return 0;
  }
  uint verbosity = 1;
  getopt(args, "v|verbosity", (string o, string v) { verbosity = to!int(v); } );

  writeln("Verbosity ",verbosity);
  return 0;
}
