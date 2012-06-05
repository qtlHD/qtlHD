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
  writeln(args);
  uint verbosity = 1;
  uint debug_level = 0;
  getopt(args, "v|verbose", (string o, string v) { verbosity = to!int(v); },
               "d|debug", (string o, string d) { debug_level = to!int(d); }
  );

  writeln("Verbosity: ",verbosity);
  writeln("Debug level: ",debug_level);
  writeln("Input files: ",args[1..$]);
  return 0;
}
