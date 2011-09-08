module qtl.core.libs.r;

import std.stdio;
import std.conv;

version (Windows) {
  import qtl.core.libs.libload;
  import std.loader;
  
  extern(C){
    double function(double, double, double, int) dnorm;
    double function(double, double, double, int, int) qf;

    void function() R_SeedsSymbol; // this is the tricky part, we need to combine calling these 3 functions to produce a new seed
    void function(int*) seed_in;
    void function(int*) seed_out;
    
    double function() norm_rand;
    double function() unif_rand;
    double function() exp_rand;
    
  }
  
  //Karl wants to bind:
  //norm_rand; set_seed
  //unif_rand; get_seed
  
  //set_seed ~= seed_in
  //get_seed ~= seed_out
  
  static this(){
    HXModule lib = load_library("R");
    load_function(dnorm)(lib,"Rf_dnorm4");
    load_function(qf)(lib,"Rf_qf");

    load_function(R_SeedsSymbol)(lib,"R_SeedsSymbol");
    load_function(seed_in)(lib,"seed_in");
    load_function(seed_out)(lib,"seed_out");

    load_function(norm_rand)(lib,"norm_rand");
    load_function(unif_rand)(lib,"unif_rand");
    load_function(exp_rand)(lib,"exp_rand");

    writeln("Loaded R functionality");
  }
  
}else{
  pragma(lib, "r");
  
  extern(C){
    double dnorm(double, double, double, int);
    double qf(double, double, double, int, int);
    void R_SeedsSymbol(); // this is the tricky part, we need to combine calling these 3 functions to produce a new seed

    void seed_in(int*);
    void seed_out(int*);

    double norm_rand();
    double unif_rand();
    double exp_rand();

  }
}

unittest{
  writeln("Unit test " ~ __FILE__);
  writeln("  - norm_rand: " ~ to!string(norm_rand()));
  writeln("  - norm_rand: " ~ to!string(norm_rand()));
  writeln("  - norm_rand: " ~ to!string(norm_rand()));
  writeln("  - unif_rand: " ~ to!string(unif_rand()));
  writeln("  - unif_rand: " ~ to!string(unif_rand()));
  writeln("  - unif_rand: " ~ to!string(unif_rand())); 
}
