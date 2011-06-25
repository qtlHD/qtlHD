/**
 * Program for converting CSV/CSVr file to XGAP and XGAP binary (xbin)
 */

module qtl.util.csv2xgap;

import qtl.core.primitives;
import qtl.core.chromosome;
import qtl.core.phenotype;
import qtl.core.genotype;
import qtl.core.xgap;
import qtl.plugins.input.read_csvr;
import qtl.plugins.output.write_xgapbin;

import std.stdio;
import std.conv;
import std.string;
import std.path;
import std.file;


void main(string[] args){
  if(args.length != 3){
    writeln("Usage: convert in.csvr out.xbin");
    writeln("  e.g. csv2xgap ../../test/data/input/multitrait.csvr test.xbin");

  }else{
    writeln("reading CSVR (" ~ args[1] ~ ") to XBIN (" ~ args[2] ~ ")");
    auto data = new CsvrReader!RIL(args[1]);
    auto result = new BinaryWriter!(CsvrReader!RIL,RIL)(data,args[2]);
    writefln("Reduced from: %.2f Kb to %.2f Kb", toKb(args[1]), toKb(args[2]));
  }
}

