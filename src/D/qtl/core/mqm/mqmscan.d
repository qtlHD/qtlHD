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
import qtl.plugins.csv.read_csv;
import qtl.core.mqm.regression;
import qtl.core.mqm.mqmutils; //Load the R-bindings
