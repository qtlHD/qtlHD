// Entry point for CLI scanone tool

import std.getopt;
import std.stdio;
import std.conv;
import std.exception;

immutable usage = "
scanone qtlHD $version (c) 2012

  usage:

";

int main(string[] args) {
  if (args.length == 1) {
    writeln(usage);
    return 0;
  }
  uint verbosity = 1;
  getopt(args, "v|verbosity", 
    delegate void(string o, string v) { verbosity = to!int(v); } );

  writeln("Verbosity ",verbosity);
  return 0;
}
