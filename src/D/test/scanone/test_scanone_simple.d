/**
 * Test scanone routines: simplest case
 */

module test.scanone.test_scanone_simple;

import std.stdio;
import std.conv;
import std.string;
import std.path;
import std.algorithm;
import std.math;
import std.algorithm;
import std.container;
import std.conv;
import std.typecons;

import qtl.plugins.qtab.read_qtab;


import qtl.core.primitives;
import qtl.core.chromosome;
import qtl.core.phenotype;
import qtl.core.marker;
import qtl.core.genotype;
import qtl.core.phenotype;
import qtl.core.map.map;
import qtl.core.map.make_map;

import qtl.core.map.genetic_map_functions;
import qtl.core.hmm.cross;
import qtl.core.hmm.calcgenoprob;
import qtl.core.scanone.scanone_hk;
import qtl.core.scanone.util;
import qtl.core.util.data_manip;


unittest {
  writeln("Unit test " ~ __FILE__);

  writeln(" --Loading data");
  // form file names
  alias std.path.buildPath buildPath;
  auto dir = to!string(dirName(__FILE__) ~ dirSeparator ~
                       buildPath("..","..","..","..","test","data", "input", "listeria_qtab"));
  string[] files = ["listeria_founder.qtab", "listeria_symbol.qtab",
                    "listeria_genotype.qtab", "listeria_phenotype.qtab",
                    "listeria_marker_map.qtab"];
  foreach(ref file; files)
    file = dir ~ dirSeparator ~ file;

  // load founder info
  auto founder_info = read_founder_settings_qtab(files[0]);
  // cross type -> Cross class
  auto cross_class = form_cross(founder_info["Cross"]);

  // load symbols and genotype data
  auto symbols = read_genotype_symbol_qtab(files[1]);
  auto g_res = read_genotype_qtab(files[2], symbols);
  auto genotype_matrix = g_res[1];

  // load phenotype data
  auto p_res = read_phenotype_qtab!(Phenotype!double)(files[3]);
  Phenotype!double[][] pheno = p_res[0];

  // load marker map
  auto markers = read_marker_map_qtab!(Marker)(files[4]);

  // omit individuals with missing phenotype
  auto ind_to_omit = is_any_phenotype_missing(pheno);
  auto n_to_omit = count(ind_to_omit, true);
  writeln(" --Omitting ", n_to_omit, " individuals with missing phenotype");
  genotype_matrix = omit_ind_from_genotypes(genotype_matrix, ind_to_omit);
  pheno = omit_ind_from_phenotypes(pheno, ind_to_omit);

  // split markers into chromsomes; sort chromosomes
  auto markers_by_chr = sort_chromosomes_by_marker_id(get_markers_by_chromosome(markers));

  // add pseudomarkers at 1.0 spacing
  auto pmar_stepped_by_chr = add_stepped_markers(markers_by_chr, 1.0);
  auto pmar_minimal_by_chr = add_minimal_markers(markers_by_chr, 1.0);
}
