module qtl.plugins.renv.rsession;

import std.stdio;
import std.conv;
import std.process;
import std.regex;
//import std.c.stdlib;
import qtl.plugins.renv.rlib;

/** 
 * Represents an R session object which is used to initialize a new session
 * The reason is that we hope that now seeds are set inside the R environment
 */
class RSession{

  this(){
    initialize();
  }

  ~this(){
    close();
  }

 /**
  * Close down R environment
  */
  void initialize() {
    char* argv[];
    argv ~= "qtlHD".dup.ptr; 
    argv ~= "--gui=none".dup.ptr;
    argv ~= "--silent".dup.ptr;
    argv ~= "--no-environ".dup.ptr;
    int argc = argv.length;

    if (!r_running) {
      writeln("Initialize embedded R (library)");
      //TODO: Execute a single R run giving the R-home
      //setenv("R_HOME","/usr/lib/R",1);
      string output = shell("R CMD config --ldflags");
      writefln("Linker flags returned by R: %s", output);
      auto m = match(output, regex("/bin"));
      rhome = m.pre[2..$];
      writefln("RHome by R: %s", rhome);
      environment["R_HOME"] = rhome;
      Rf_initEmbeddedR(argc, argv.ptr);
      r_running = true;
    }
  }
  
  int callFunction(){
    SEXP e, val;
    int errorOccurred;
    int result = -1;

    Rf_protect(e = Rf_allocVector(LANGSXP, 1));
    SETCAR(e, Rf_install("sum\0".dup.ptr));

    val = R_tryEval(e, R_GlobalEnv, &errorOccurred);

    if(!errorOccurred) {
      Rf_protect(val);
      result = INTEGER(val)[0];
      Rf_unprotect(1);
    }else{
      writeln("An error occurred when calling sum()");
    }
    Rf_unprotect(1); /* e / Assume we have an INTSXP here. */
    return(result);
  }

 /**
  * Close down R environment
  */
  void close() {
    if (r_running) {
      Rf_endEmbeddedR(0);
      r_running = false;
    }
  }
  
private:
  bool r_running;
  string rhome = "/usr/lib/R";
}

unittest{
  RSession r = new RSession();
  r.callFunction();
  //writeln("  - norm_rand: " ~ to!string(norm_rand()));
  //writeln("  - norm_rand: " ~ to!string(norm_rand()));
  //writeln("  - norm_rand: " ~ to!string(norm_rand()));
  r.close();
}
