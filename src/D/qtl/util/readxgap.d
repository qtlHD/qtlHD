/**
 * Program for reading an XGAP binary (xbin)
 */

module qtl.util.readxgap;

import qtl.core.primitives;
import qtl.core.chromosome;
import qtl.core.phenotype;
import qtl.core.deprecate.genotype_enum;
import qtl.plugins.xgap.xgap;
import qtl.plugins.xgap.read_xgapbin;

import std.stdio;
import std.conv;
import std.string;
import std.path;
import std.file;

void print_help(){
  writeln("Usage: readxgap type infile");
  writeln("  infile: XGAP binary (.xbin)");
  writeln("  type:   RIL, F2, BC");
  writeln("  e.g. readxgap RIL multitrait.xbin");
}

/*
 * Function that creates a reader (bin) and provides a template for using different genotype enums
 */

void read_xgap(XType)(string filein){
  auto data = new XbinReader!XType(filein);
}

void main(string[] args){
  if(args.length != 3){
    print_help();
  }else{
    if(args[1] == "RIL") read_xgap!RIL(args[2]);
    if(args[1] == "F2")  read_xgap!F2(args[2]);
    if(args[1] == "BC")  read_xgap!F2(args[2]);
  }
}

