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

  // load all qtab files
  auto all_data = load_qtab(files);
  auto founder_info = all_data[1];
  auto markers = all_data[2];
  auto pheno = all_data[4];
  auto genotype_matrix = all_data[6];

  auto cross_class = form_cross(founder_info["Cross"]);

  // omit individuals with missing phenotype
  auto ind_to_omit = is_any_phenotype_missing(pheno);
  auto n_to_omit = count(ind_to_omit, true);
  writeln(" --Omitting ", n_to_omit, " individuals with missing phenotype");
  genotype_matrix = omit_ind_from_genotypes(genotype_matrix, ind_to_omit);
  pheno = omit_ind_from_phenotypes(pheno, ind_to_omit);

  // split markers into chromsomes; sort chromosomes
  auto markers_by_chr = sort_chromosomes_by_marker_id(get_markers_by_chromosome(markers));

  // drop X chromosome
  markers_by_chr = markers_by_chr[0..($-1)];

  // add pseudomarkers at 1.0 cM spacing
  auto pmar_by_chr = add_stepped_markers(markers_by_chr, 1.0);

  // inter-marker recombination fractions
  auto rec_frac = recombination_fractions(pmar_by_chr, GeneticMapFunc.Haldane);

  // empty covariate matrices
  auto addcovar = new double[][](0, 0);
  auto intcovar = new double[][](0, 0);
  auto weights = new double[](0);

  // null model
  auto rss0 = scanone_hk_null(pheno, addcovar, weights);

  // to store LOD curves and peaks for all chromosomes
  double[][] lod;
  Tuple!(double, Marker)[][] peaks;

 // calc genoprob for each chromosome, then scanone
  foreach(i, chr; pmar_by_chr) {
    auto genoprobs = calc_geno_prob(cross_class, genotype_matrix, chr[1], rec_frac[i][0], 0.002);
    auto rss = scanone_hk(genoprobs, pheno, addcovar, intcovar, weights);
    auto lod_this_chr = rss_to_lod(rss, rss0, pheno.length);
    lod ~= lod_this_chr;

    auto peak_this_chr = get_peak_scanone(lod_this_chr, chr[1]);
    peaks ~= peak_this_chr;
  }

  // print peaks
  double threshold = 2;
  writeln(" --Peaks with LOD > ", threshold, ":");
  foreach(peak; peaks) {
    foreach(j; 0..peak.length) {
      if(peak[j][0] > threshold)
        writefln(" ----Chr %-2s : peak for phenotype %d: max lod = %7.2f at pos = %7.2f", peak[j][1].chromosome.name, j,
                 peak[j][0], peak[j][1].get_position);
    }
  }
}
