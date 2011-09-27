/**
 * Multiple QTL modelling and scanning routine
 **/

module qtl.core.mqm.mqmscan;
 
import std.container;
import qtl.core.primitives;
import std.stdio;
import std.algorithm;
import std.math;
import std.string;
import qtl.core.deprecate.genotype_enum;
import qtl.plugins.input.read_csv;

//Load the R-bindings
import qtl.core.mqm.mqmutils;
