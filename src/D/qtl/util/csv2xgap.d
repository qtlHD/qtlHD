/**
 * Program for converting CSV/CSVr file formats to XGAP and XGAP binary (xbin)
 */

module qtl.util.csv2xgap;

import qtl.core.primitives;
import qtl.core.chromosome;
import qtl.core.phenotype;
import qtl.core.genotype;
import qtl.core.xgap;
import qtl.plugins.input.read_csvr;
import qtl.plugins.input.read_csv;
import qtl.plugins.input.read_interface;
import qtl.plugins.output.write_xgapbin;

import std.stdio;
import std.conv;
import std.string;
import std.path;
import std.file;

void print_help(){
    writeln("Usage: cvs2xgap type infile outfile");
    writeln("  infile: CSV (.csv), CSVr (.csvr)");
    writeln("  type:   RIL, F2, BC");
    writeln("  e.g. csv2xgap RIL ../../test/data/input/multitrait.csvr multitrait.xbin");
    writeln("       csv2xgap F2 ../../test/data/input/hyper.csv hyper.xbin");
}

/*
 * Function that creates a reader (csv or csvr) and provides a template for using different genotype enums
 */

void read_type(XType)(string filein, string fileout){
  if(filein.lastIndexOf(".") > 0){
    string extension = filein[filein.lastIndexOf(".")+1..$];
    GenericReader!XType data;
    if(extension.tolower() == "csv"){
      data = cast(GenericReader!XType)new ReadSimpleCSV!XType(filein);
    }
    if(extension.tolower() == "csvr"){
      data = cast(GenericReader!XType)new CSVrReader!XType(filein);
    }
    writeln("reading CSVR (" ~ filein ~ ") to XBIN (" ~ fileout ~ ")");
    auto result = new BinaryWriter!(GenericReader!XType,XType)(data,fileout);
    writefln("Reduced from: %.2f Kb to %.2f Kb", toKb(filein), toKb(fileout));
  }
}

void main(string[] args){
  if(args.length != 4){
    print_help();
  }else{
    if(args[1] == "RIL") read_type!RIL(args[2],args[3]);
    if(args[1] == "F2") read_type!F2(args[2],args[3]);
    if(args[1] == "BC") read_type!F2(args[2],args[3]);
  }
}

