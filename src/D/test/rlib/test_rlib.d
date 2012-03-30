/**

  A minimal example to use the random number generator of R (calling into Rlib.so from D)

  dmd test_rlib.d ; ./test_rlib

*/

module test.rlib.test_rlib;

import std.algorithm;
import std.range;
import std.c.stdio;
import std.stdio;
import std.string;
import std.conv;

import std.c.stdlib;

pragma(lib, "R");

extern (C) void Rf_initEmbeddedR(size_t argc, char **argv);
extern (C) void Rf_endEmbeddedR(size_t);

extern(C){
    double norm_rand();
    double unif_rand();
    double exp_rand();
    void GetRNGstate();
    void PutRNGstate();
}

/**
 * Initialize the R interpreter
 */

void R_Init() {
  version(darwin) { // OSX
    string args[] = [ "BiolibEmbeddedR", "--gui=none", "--silent", "--no-environ", "--no-site-file", "--no-init-file"];
  } else {
    string args[] = [ "BiolibEmbeddedR", "--gui=none", "--silent", "--no-environ"];
  }
  char *argv[];
 
  argv.length = args.length;
  foreach (i, s ; args) {
    argv[i] = cast(char *)toStringz(args[i]);
  }
  auto argc = args.length;

  writeln("Initialize embedded R (library)");
  version(darwin) {
    setenv("R_HOME","/Library/Frameworks/R.framework/Resources/", 1);
  } else version(linux) {
    setenv("R_HOME","/usr/lib/R",1);
  } else {
    throw new Exception("Can not find R libraries on this system");
  }
  Rf_initEmbeddedR(argc, argv.ptr);
}

void R_Close() {
  writeln("Shutting down R");
  PutRNGstate();
  Rf_endEmbeddedR(0);
}

int main()
{

  R_Init();
  
  // call an Rlib function
  GetRNGstate();

  writeln("  - norm_rand: " ~ to!string(norm_rand()));
  writeln("  - norm_rand: " ~ to!string(norm_rand()));
  writeln("  - norm_rand: " ~ to!string(norm_rand()));
  writeln("  - unif_rand: " ~ to!string(unif_rand()));
  writeln("  - unif_rand: " ~ to!string(unif_rand()));
  writeln("  - unif_rand: " ~ to!string(unif_rand())); 

  R_Close();
  return 0;
}
