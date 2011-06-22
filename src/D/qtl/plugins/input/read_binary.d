/**
 * \file read_binary.d - Plugin for writing XBIN files
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
 * dmd -unittest qtl/plugins/input/read_binary.d qtl/core/*.d  
 *
 * - Win32 (No circular dependancies, or multiple definitions of main allowed):
 *
 * dmd -run cdc.d -lib qtl/core/ -ofCore.lib
 * dmd -run cdc.d -unittest qtl/plugins/input/read_binary.d Core.lib
 *
 *
 **/
 
 module qtl.plugins.input.read_binary;

import qtl.core.primitives;
import qtl.core.chromosome;
import qtl.core.phenotype;
import qtl.core.genotype;

import std.stdio;
import std.conv;
import std.string;
import std.path;

class XbinReader(XType){
  private File f;
  int nindividuals = 0;
  int nphenotypes = 0;
  int nmarkers = 0;
  
  this(in string filename){
  
  }
}

unittest{
  writeln("Unit test " ~ __FILE__);
  alias std.path.join join;
  auto infn = dirname(__FILE__) ~ sep ~ join("..","..","..","..","..","test","data","input","multitrait.xbin");
  writeln("  - reading XBIN " ~ infn);
  auto data = new XbinReader(infn);
}
