/**
 * Individual related functions
 */


module qtl.core.individual;

import std.container;
import std.conv;
import std.variant;
import std.random;
import qtl.core.primitives;
import std.algorithm; // to use sort()

import std.stdio;
import std.typecons;

class Individuals {
  mixin ActList!Individual;
  /*
  Individuals opOpAssign(string op)(string name) if (op == "~") {
    list ~= new Individual(name);
    return this;
  }

  Individuals opOpAssign(string op)(Individual ind) if (op == "~") {
    // list ~= ind.dup;
    return this;
  }
  */
}



