/**
 * Interface class for the csv and csvr readers
 **/

module qtl.plugins.input.read_interface;

import qtl.core.primitives;
import qtl.core.chromosome;
import qtl.core.phenotype;
import qtl.core.genotype;

abstract class GenericReader(XType){
  string[] phenotypenames;
  Marker[] markers;
  Chromosome[string] chromosomes;
  Phenotype!double[][] phenotypes;
  Genotype!XType[][] genotypes;
  
  this(){
  
  }
  
  this(in string infilename){
  
  }
}