/**

  dmd test_rlib.d

*/

module test.rlib.test_rlib;

import std.algorithm;
import std.range;
import std.c.stdio;
import std.stdio;
import std.string;

import std.c.stdlib;

pragma(lib, "R");

extern (C) void Rf_initEmbeddedR(int argc, char **argv);
extern (C) void Rf_endEmbeddedR(int);

void BioLib_R_Init() {
  string args[] = [ "BiolibEmbeddedR", "--gui=none", "--silent", "--no-environ" ];
  char *argv[];
 
  argv.length = args.length;
  foreach (i, s ; args) {
    argv[i] = cast(char *)toStringz(args[i]);
  }
  int argc = args.length;

  writeln("Initialize embedded R (library)");
  setenv("R_HOME","/usr/lib/R",1);
  Rf_initEmbeddedR(argc, argv.ptr);
  writeln("Shutting down R");
  Rf_endEmbeddedR(0);
}

int main()
{

  BioLib_R_Init();

  return 0;
}
