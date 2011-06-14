/**
 * \file csvrtobinary.d - Plugin for converting CSVR files to XBIN
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
 * Test: dmd -unittest qtl/plugins/input/csvrtobinary.d qtl/core/*.d
 *
 **/

module qtl.plugins.input.csvrtobinary;

import qtl.core.primitives;
import qtl.core.chromosome;
import qtl.core.phenotype;
import qtl.core.genotype;

import std.stdio;
import std.conv;
import std.string;
import std.path;

struct GeneticMarker(T) {
  string      name;
  string      chr;
  bool        autosomal;
  double      dist;
  T[]         values;
}

struct NPhenotype(T) {
  string      name;
  T[]         values;
}

/** 
 * Convert a simple CSVR file containing marker names, chromosome nrs, position, 
 * phenotype and genotype - such as the multitrait.csvr file used in R/qtl.
 *
 * The file is parsed once on class instantiation. Elements can be queried.
 */
class CsvrToBinary {
  private File f;
  int individuals = 0;
  int phenotypes = 0;
  int markers = 0;
  
  NPhenotype!double phenome[];
  GeneticMarker!char genome[];
  
  bool check_individuals(string[] items){
    if(individuals==0){
      individuals = (items.length-3);
      return true;
    }else{
      return(individuals==(items.length-3));
    }
  }
  
  bool is_phenotype(string[] items){
   return(items[1] == "" && items[2] == "");
  }
  
  void read_csvr(in string filename){
    f = File(filename,"r");
    string line;
    int linecount;
    while((line = f.readln()) != ""){
      auto items = split(line,",");
      if(!check_individuals(items))throw new Exception("Not enough items on line " ~ to!string(linecount));
      if(is_phenotype(items)){
        writeln("Phenotype: " ~ items[0]);
        NPhenotype!double phenotype;
        phenotype.name = items[0];
        for(uint i=0;i<3;i++) {
          if(items[i+3] != "-"){
            phenotype.values ~= to!double(items[i+3]);
          }else{
            phenotype.values ~= PHENOTYPE_NA;
          }
        }
        phenome ~= phenotype;
        phenotypes++;
      }else{
        writeln("Marker: " ~ items[0]);
        GeneticMarker!char marker;
        marker.name = items[0];
        marker.chr = items[1];
        if(items[1] == "X"){
          marker.autosomal=true;
        }else{
          marker.autosomal=false;
        }
        marker.dist = to!double(items[2]);
        for(int i=0;i<individuals;i++) {
          marker.values ~= to!string(items[i+3]);
        }
        genome ~= marker;
        markers++;
      }
      linecount++;
    }
    writefln("Read %s with %d phenotypes and %d markers measured at %d individuals",filename,phenotypes,markers,individuals);
    f.close();
  }
  
  void write_binary(in string filename){
    f = File(filename,"wb");

    f.close();
  }

  this(in string infilename, in string outfilename){
    read_csvr(infilename);
    write_binary(outfilename);
  }
}

unittest {
  writeln("Unit test " ~ __FILE__);
  alias std.path.join join;
  auto infn = dirname(__FILE__) ~ sep ~ join("..","..","..","..","..","test","data","input","multitrait.csvr");
  auto outfn = dirname(__FILE__) ~ sep ~ join("..","..","..","..","..","test","data","output","multitrait.xbin");
  writeln("  - converting CSVR " ~ infn ~ " to " ~ outfn);
  auto data = new CsvrToBinary(infn, outfn);

}

void main(string[] args){ 
  if(args.length != 3){
    writeln("Usage: convert in.csvr out.xbin");
  }else{
    writeln("Converting CSVR (" ~ args[1] ~ ") to XBIN (" ~ args[2] ~ ")");
    auto data = new CsvrToBinary(args[1],args[2]);
  }
}
