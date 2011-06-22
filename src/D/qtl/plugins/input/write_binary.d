/**
 * \file write_binary.d - Plugin for writing XBIN files
 *
 * Copyright (c) 2011 Danny Arends
 * Part of the qtlHD package
 *
 *     This program is free software; you can redistribute it and/or
 *     modify it under the terms of the GNU General Public License,
 *     version 3, as published by the Free Software Foundation.
 * 
 *     This program is distributed in the hope that it will be useful,
 *     but without any warranty; without even the implied warranty of
 *     merchantability or fitness for a particular purpose.  See the GNU
 *     General Public License, version 3, for more details.
 * 
 *     A copy of the GNU General Public License, version 3, is available
 *     at http://www.r-project.org/Licenses/GPL-3
 *
 * Written in the D Programming Language (http://www.digitalmars.com/d)
 *
 * - Mac / Linux Unittest:
 * dmd -unittest qtl/plugins/input/write_binary.d qtl/core/*.d  
 *
 * - Win32 (No circular dependancies, or multiple definitions of main allowed):
 *
 * dmd -run cdc.d -lib qtl/core/ -ofCore.lib
 * dmd -run cdc.d -unittest qtl/plugins/input/write_binary.d Core.lib
 *
 *
 **/

module qtl.plugins.input.write_binary;

import qtl.core.primitives;
import qtl.core.chromosome;
import qtl.core.phenotype;
import qtl.core.genotype;

import std.stdio;
import std.conv;
import std.string;
import std.path;

import qtl.plugins.input.read_csvr;

/** 
 * Convert a simple CSVR file containing marker names, chromosome nrs, position, 
 * phenotype and genotype - such as the multitrait.csvr file used in R/qtl.
 *
 * The file is parsed once on class instantiation. Elements can be queried.
 */
class BinaryWriter {
  private File f;
  CsvrReader data;
  
  byte[2] b_footprint = [ 0, 5 ];
  byte[3] b_version = [ 0, 0, 1 ];
  
  void myWrite(T)(T x,File f,bool bin=true){
    if(bin){
      f.rawWrite(x);
    }else{
      f.writeln(x);
    }
  }
  
  void write_matrix(T)(T[] towrite, File to){
    uint[2] sizes =[towrite.length,towrite[0].values.length];
    myWrite(sizes,to);
    foreach(T e;towrite){
      myWrite(e.values,to);
    }
  }
  
  void write_binary(in string filename){
    f = File(filename,"wb");
    myWrite(b_footprint,f);
    myWrite(b_version,f);
    myWrite(b_footprint,f);
    uint[1] nmatrix = [ 2 ];
    myWrite(nmatrix,f);
    write_matrix!(NPhenotype!double)(data.phenotypes,f);
    myWrite(b_footprint,f);
    write_matrix!(GeneticMarker!string)(data.genotypes,f);
    myWrite(b_footprint,f);
    f.close();
  }

  this(CsvrReader indata, in string outfilename){
    data = indata
    write_binary(outfilename);
  }
}

unittest {
  writeln("Unit test " ~ __FILE__);
  alias std.path.join join;
  auto infn = dirname(__FILE__) ~ sep ~ join("..","..","..","..","..","test","data","input","multitrait.csvr");
  auto outfn = dirname(__FILE__) ~ sep ~ join("..","..","..","..","..","test","data","input","multitrait.xbin");
  writeln("  - reading CSVR " ~ infn ~" to " ~ outfn);
  auto data = new CsvrReader(infn);
  BinaryWriter(data,outfn);
}

void main(string[] args){ 
  if(args.length != 3){
    writeln("Usage: convert in.csvr out.xbin");
  }else{
    writeln("reading CSVR (" ~ args[1] ~ ") to XBIN (" ~ args[2] ~ ")");
    auto data = new CsvrReader(args[1]);
    BinaryWriter(data,args[2]);
  }
}

