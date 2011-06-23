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

import qtl.plugins.input.binary_types;
import qtl.plugins.input.read_csvr;

/** 
 * Convert a simple CSVR file containing marker names, chromosome nrs, position, 
 * phenotype and genotype - such as the multitrait.csvr file used in R/qtl.
 *
 * The file is parsed once on class instantiation. Elements can be queried.
 */
class BinaryWriter(XType) {
  private File f;
  CsvrReader!XType data;
 
  void myWrite(T)(T x,File f,bool bin = true, MatrixType t = MatrixType.EMPTY){
    if(bin){
      f.rawWrite(x);
    }else{
      f.writeln(x);
    }
  }
  
  void write_matrix(T)(T[][] towrite, File outfile, MatrixType t = MatrixType.EMPTY){
    uint[1] type = [10];
    uint[] sizes =[towrite.length,towrite[0].length];
    switch(t){
      case MatrixType.EMPTY:
        sizes ~= cast(uint)0;
        break;
      case MatrixType.INTMATRIX:
        sizes ~= int.sizeof;
        break;        
      case MatrixType.DOUBLEMATRIX:
        sizes ~= double.sizeof;
        break;        
      case MatrixType.FIXEDCHARMATRIX:
        string s = to!string(towrite[0][0]);
        sizes ~= char.sizeof * s.length;
        break;        
      case MatrixType.VARCHARMATRIX:
        for(int r=0;r < towrite.length;r++){
          for(int c=0;c < towrite[r].length;c++){
            sizes ~= char.sizeof * to!string(towrite[r][c]).length;
          }
        }
        break;        
      default:
        break;
    }
    myWrite(type,outfile);
    myWrite(sizes,outfile);
    foreach(T[] e;towrite){
      myWrite(e,outfile);
    }
  }
  
  string[][] getMarkerInfoMatrix(Marker[] markers){
    string[][] return_matrix;
    foreach(Marker m;markers){
      string[] row;
      row ~= m.name;
      row ~= m.chromosome.name;
      row ~= to!string(m.position);
      return_matrix ~= row;
    }
    return return_matrix;
  }
  
  void write_binary(in string filename){
    f = File(filename,"wb");
    myWrite(b_footprint,f);
    myWrite(b_version,f);
    uint[1] nmatrix = [ 3 ];
    myWrite(nmatrix,f);
    write_matrix!(Phenotype!double)(data.phenotypes, f, MatrixType.DOUBLEMATRIX);
    myWrite(b_footprint,f);
    write_matrix!(Genotype!XType)(data.genotypes,f, MatrixType.FIXEDCHARMATRIX);
    myWrite(b_footprint,f);
    auto markermatrix = getMarkerInfoMatrix(data.markers);
    write_matrix!(string)(markermatrix,f, MatrixType.VARCHARMATRIX);
    myWrite(b_footprint,f);
    f.close();
  }

  this(CsvrReader!XType indata, in string outfilename){
    data = indata;
    write_binary(outfilename);
  }
}

unittest {
  writeln("Unit test " ~ __FILE__);
  alias std.path.join join;
  auto infn = dirname(__FILE__) ~ sep ~ join("..","..","..","..","..","test","data","input","multitrait.csvr");
  auto outfn = dirname(__FILE__) ~ sep ~ join("..","..","..","..","..","test","data","input","multitrait.xbin");
  writeln("  - reading CSVR " ~ infn ~" to " ~ outfn);
  auto data = new CsvrReader!RIL(infn);
  auto result = new BinaryWriter!RIL(data,outfn);
}

void main(string[] args){ 
  if(args.length != 3){
    writeln("Usage: convert in.csvr out.xbin");
  }else{
    writeln("reading CSVR (" ~ args[1] ~ ") to XBIN (" ~ args[2] ~ ")");
    auto data = new CsvrReader!RIL(args[1]);
    auto result = new BinaryWriter!RIL(data,args[2]);
  }
}

