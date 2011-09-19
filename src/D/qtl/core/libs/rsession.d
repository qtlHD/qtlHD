module qtl.core.libs.rsession;

import std.stdio;
import std.conv;
import std.process;
//import std.c.stdlib;
import qtl.core.libs.rlib;

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
      writeln(output);
      environment["R_HOME"] = rhome;
 //     Rf_initEmbeddedR(argc, argv.ptr);
      r_running = true;
    }
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
  writeln("  - norm_rand: " ~ to!string(norm_rand()));
  writeln("  - norm_rand: " ~ to!string(norm_rand()));
  writeln("  - norm_rand: " ~ to!string(norm_rand()));
  r.close();
}
