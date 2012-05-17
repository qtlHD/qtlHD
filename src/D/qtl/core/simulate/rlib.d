/**
  Code to initialize R, for random number generator (calling into Rlib.so from D)
*/

module qtl.core.simulate.rlib;

import std.algorithm;
import std.range;
import std.c.stdio;
import std.stdio;
import std.string;
import std.conv;

import std.c.stdlib;

version(Windows){
  // Windows library binding
  private import std.loader;
  private import arch.windows;
  private import qtl.plugins.renv.libload;
  
  extern (C){
    void function(size_t argc, char **argv) Rf_initEmbeddedR;
    void function(size_t) Rf_endEmbeddedR;
    void function() R_SaveGlobalEnv;
  
    double function() norm_rand;
    double function() unif_rand;
    double function() exp_rand;
    double function(double) rpois;
    void function() GetRNGstate;
    void function() PutRNGstate;
  }

  // load the functions using libload.d
  static this(){
    HXModule lib = load_library("R");
    load_function(norm_rand)(lib,"norm_rand");
    load_function(unif_rand)(lib,"unif_rand");
    load_function(exp_rand)(lib,"exp_rand");
    load_function(Rf_Rf_rpois)(lib,"Rf_Rf_rpois");
    load_function(GetRNGstate)(lib,"GetRNGstate");
    load_function(PutRNGstate)(lib,"PutRNGstate");
    
    load_function(Rf_initEmbeddedR)(lib,"Rf_initEmbeddedR");
    load_function(Rf_endEmbeddedR)(lib,"Rf_endEmbeddedR");
    load_function(R_SaveGlobalEnv)(lib,"R_SaveGlobalEnv");
    writeln("Loaded R functionality");
  }
}else{
  pragma(lib, "R");

  extern (C) void Rf_initEmbeddedR(size_t argc, char **argv);
  extern (C) void Rf_endEmbeddedR(size_t);
  extern (C) void R_SaveGlobalEnv();

  extern(C){
      double norm_rand();
      double unif_rand();
      double exp_rand();
      double Rf_rpois(double);
      void GetRNGstate();
      void PutRNGstate();
  }
}
/**
 * Initialize the R interpreter
 */

void R_Init() {
version(darwin) { // OSX
    string args[] = [ "BiolibEmbeddedR", "--gui=none", "--silent", "--no-environ", "--no-site-file", "--no-init-file"];
  }else version(Windows){ // Win
    string args[] = [ "BiolibEmbeddedR", "--gui=none", "--silent", "--no-environ"];
  }else { // Unix
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
  } else version(Windows) { //Windows version, perhaps a cmdline arg ?
    setEnv("R_HOME","C:/Program Files/R/R-2.14.1/",1);
  }else {
    throw new Exception("Can not find R libraries on this system");
  }
  Rf_initEmbeddedR(argc, argv.ptr);

  // call an Rlib function
  GetRNGstate();
}

void R_Close() {
  writeln("Shutting down R");
  PutRNGstate();
  R_SaveGlobalEnv();
  Rf_endEmbeddedR(0);
}

// Rf_rpois but with output as integer (which it should be anyway)
int rpois(double mu)
{
  return cast(int)Rf_rpois(mu);
}

/* comment out; for now, can only do R_Init() and R_Close() once
unittest {
  R_Init();
  
  writeln("  - norm_rand: " ~ to!string(norm_rand()));
  writeln("  - norm_rand: " ~ to!string(norm_rand()));
  writeln("  - norm_rand: " ~ to!string(norm_rand()));
  writeln("  - unif_rand: " ~ to!string(unif_rand()));
  writeln("  - unif_rand: " ~ to!string(unif_rand()));
  writeln("  - unif_rand: " ~ to!string(unif_rand())); 
  writeln("  ~ rpois(5):     " ~ to!string(rpois(5.0)));
  writeln("  ~ rpois(2.5):   " ~ to!string(rpois(2.5)));
  writeln("  ~ rpois(15.6):  " ~ to!string(rpois(15.6)));
  writeln("  ~ rpois(50.1):  " ~ to!string(rpois(50.1)));

  R_Close();
}
*/