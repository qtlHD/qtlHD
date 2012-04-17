/**
 * Data manipulation
 */

module qtl.core.util.data_manip;

import std.conv, std.stdio, std.string, std.exception;
import qtl.core.primitives;
import qtl.core.phenotype;
import qtl.core.genotype;

// return boolean vector indicating whether any
bool[] is_any_phenotype_missing(T)(Phenotype!T[][] pheno)
{
  auto ret = new bool[](pheno.length);
  foreach(i; 0..pheno.length) {
    ret[i] = false;
    foreach(j; 0..pheno[i].length) {
      if(isNA(pheno[i][j])) ret[i] = true;
      break;
    }
  }
  return ret;
}

// omit individuals from phenotype data
Phenotype!T[][] omit_ind_from_phenotypes(T)(Phenotype!T[][] pheno, bool[] to_omit)
{
  if(pheno.length != to_omit.length)
    throw new Exception("no. individuals in pheno (" ~ to!string(pheno.length) ~ 
                        ") doesn't match length of to_omit (" ~ to!string(to_omit.length) ~ ")");

  Phenotype!T[][] ret;

  foreach(i; 0..to_omit.length) {
    if(!to_omit[i])
      ret ~= pheno[i];
  }

  return ret;
}

// omit individuals from genotype matrix
GenotypeCombinator[][] omit_ind_from_genotypes(GenotypeCombinator[][] geno, bool[] to_omit)
{
  if(geno.length != to_omit.length)
    throw new Exception("no. individuals in geno (" ~ to!string(geno.length) ~ 
                        ") doesn't match length of to_omit (" ~ to!string(to_omit.length) ~ ")");

  GenotypeCombinator[][] ret;

  foreach(i; 0..to_omit.length) {
    if(!to_omit[i])
      ret ~= geno[i];
  }

  return ret;
}

