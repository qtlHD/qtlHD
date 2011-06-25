/**
 * \file read_csvr.d - Plugin for reading CSVR files
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
 * dmd -unittest qtl/plugins/input/read_csvr.d qtl/core/*.d  
 *
 * - Win32 (No circular dependancies, or multiple definitions of main allowed):
 *
 * dmd -run cdc.d -lib qtl/core/ -ofCore.lib
 * dmd -run cdc.d -unittest qtl/plugins/input/read_csvr.d Core.lib
 *
 **/

module qtl.plugins.input.read_csvr;

import qtl.core.primitives;
import qtl.core.chromosome;
import qtl.core.phenotype;
import qtl.core.genotype;

import std.stdio;
import std.conv;
import std.string;
import std.path;

/** 
 * Loads a simple CSVR file containing marker names, chromosome nrs, position, 
 * phenotype and genotype - such as the multitrait.csvr file used in R/qtl.
 *
 * The file is parsed once on class instantiation. Elements can be queried.
 */
class CsvrReader(XType){
  private File f;
  int nindividuals = 0;
  int nphenotypes = 0;
  int nmarkers = 0;
  
  string[] phenotypenames;
  Marker[] markers;
  Chromosome[string] chromosomes;
  Phenotype!double[][] phenotypes;
  Genotype!XType[][] genotypes;

  
  bool check_individuals(string[] items){
    if(nindividuals==0){
      nindividuals = (items.length-3);
      return true;
    }else{
      return(nindividuals==(items.length-3));
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
      line = strip(line);
      auto items = split(line,",");
      if(!check_individuals(items))throw new Exception("Not enough items on line " ~ to!(string)(linecount));
      if(is_phenotype(items)){
        //Phenotype
        debug writeln("Phenotype: " ~ items[0]);
        Phenotype!double[] ps;
        phenotypenames ~= items[0];
        for(uint i=0;i<nindividuals;i++){
          ps ~= set_phenotype!double(items[i+3]);
        }
        phenotypes ~= ps;
        nphenotypes++;
      }else{
        debug writeln("Marker: " ~ items[0]);
        //CHR
        if (!(items[1] in chromosomes)){
          chromosomes[items[1]] = get_chromosome_with_id(items[1]);
        }

        //Marker
        Marker m = new Marker(to!double(items[2]),items[0]); //
        m.chromosome = chromosomes[items[1]];
        markers ~= m;

        //Genotype
        Genotype!XType[] gs;
        for(int i=0;i<nindividuals;i++) {
          gs ~= set_genotype!(XType)(items[i+3]);
        }
        genotypes ~= gs;
        nmarkers++;
      }
      linecount++;
    }
    writefln("Read %s with %d phenotypes and %d markers measured at %d individuals",filename,nphenotypes,nmarkers,nindividuals);
    f.close();
  }

  this(in string infilename){
    read_csvr(infilename);
  }
}

unittest{
  writeln("Unit test " ~ __FILE__);
  alias std.path.join join;
  auto infn = dirname(__FILE__) ~ sep ~ join("..","..","..","..","..","test","data","input","multitrait.csvr");
  writeln("  - reading CSVR " ~ infn);
  auto data = new CsvrReader!RIL(infn);
}

