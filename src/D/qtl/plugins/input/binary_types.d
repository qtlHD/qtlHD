/**
 * \file binary_types.d - Plugin for writing XBIN files
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
 **/
module qtl.plugins.input.binary_types; 
 
 enum MatrixType : uint { 
  EMPTY = 0, 
  INTMATRIX = 1, 
  DOUBLEMATRIX = 2, 
  FIXEDCHARMATRIX = 3, 
  VARCHARMATRIX = 4
};

byte[2] b_footprint = [ 0, 5 ];
byte[3] b_version = [ 0, 0, 1 ];

