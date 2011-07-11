/**
 * Test scanone routines, using listeria (CSV) set
 */

module test.scanone.test_scanone;

import qtl.core.primitives;
import qtl.core.chromosome;
import qtl.core.phenotype;
import qtl.core.genotype;
import qtl.core.marker;
import qtl.core.map;
import qtl.plugins.input.read_csv;
import std.stdio;
import std.conv;
import std.string;
import std.path;
import std.algorithm;

// The comments are based on the Ruby/Biolib-R/qtl integration

static bool VERBOSE = false;

unittest {

  writeln("Unit test " ~ __FILE__);
  alias std.path.join join;
  auto fn = dirname(__FILE__) ~ sep ~ join("..","..","..","..","test","data","input","listeria.csv");
  writeln("  - reading CSV " ~ fn);
  Marker m2 = new Marker(4.8);
  auto markers = [ m2 ];
  auto data = new ReadSimpleCSV!F2(fn);
  assert(data.markers.length == 133, to!string(data.markers.length));
  assert(data.phenotypenames[0] == "T264");
  assert(data.markers[0].name == "D10M44");
  assert(data.markers[0].id == 0);
  assert(data.markers[1].id == 1);
  // Check chromosomes
  assert(data.chromosomes.length == 20, to!string(data.chromosomes.length));
  assert(data.chromosomes["X"].id == ID_UNKNOWN);
  assert(data.chromosomes["7"].id == 7);
  assert(data.markers[2].position == 24.84773, "Marker position not matching");
  // Check phenotype
  assert(data.phenotypes[29][0].value == PHENOTYPE_NA, to!string(data.phenotypes[29][0].value));
  assert(data.phenotypes[30][0].value == 74.417);
  // Check genotype
  assert(data.genotypes[1][0].value == F2.NA);
  assert(data.genotypes[1][1].value == F2.B);
  assert(data.individuals.length == 120);

  // Markers per chromosome
  /*
    => [["1", 13], ["10", 5], ["11", 6], ["12", 6], ["13", 12], ["14", 4], ["15", 8], ["16", 4], ["17", 4], ["18", 4], ["19", 4], ["2", 6], ["3", 6], ["4", 4], ["5", 13], ["6", 13], ["7", 6], ["8", 6], ["9", 7], ["X", 2]]
  */
  auto c_mslist = get_markers_by_chromosome(data.markers);
  foreach(c_ms ; c_mslist) {
    auto c = c_ms[0];  // first part of Tuple
    auto ms = c_ms[1]; // second part of Tuple
    // look for X chromosome with 2 markers
    if (c.name == "X") {
      assert(ms.length == 2, to!string(ms.length));
      assert(ms[0].position == 0, to!string(ms[0].position));
      assert(ms[0].toString == "DXM186~0", ms[0].toString);
      // writeln(ms[1]);
      assert(ms[1].name == "DXM64", ms[1].name);
      assert(ms[1].toString == "DXM64~42.3459", ms[1].toString);
    }
    if (VERBOSE) {
      writeln(c.name);
      foreach (m; ms) {
        writeln(m.name,'-',m.position);
      }
    }
  }
  // Find marker by name
  Marker m = find!("a.name == \"D10M44\"")(data.markers)[0];
  assert(m.name == "D10M44");
  assert(m.id == 0);
  // Find marker by id
  Marker m1 = find!("a.id == 1")(data.markers)[0];
  assert(m1.id == 1);
  assert(m1.position == 0.99675);
  Marker m14 = find!("a.id == 14")(data.markers)[0];
  assert(m14.chromosome.name == "2");

  // Calculate map size X Chromosome
  auto cx = find!("a[0].name == \"X\"")(c_mslist)[0];
  auto msx = cx[1];
  assert(msx.length == 2);
  assert(to!string(map_size(msx))=="42.3459",to!string(map_size(msx)));
  // Calculate map size all Chromosomes
  double totalsize = 0.0;
  foreach(c_ms ; c_mslist) {
    totalsize += map_size(c_ms[1]);
  }
  assert(to!string(totalsize) == "1104.29");

  // Chromosome 1
  auto c1 = find!("a[0].name == \"1\"")(c_mslist)[0];
  auto ms1 = c1[1].sort;
  writeln(ms1);
  auto rfs = recombination_fractions(ms1);
  assert(rfs.length==12);
  assert(to!string(rfs[0])=="0.00986881");
  assert(to!string(rfs[1])=="0.189685");
  // Chromosome X
  writeln(recombination_fractions(msx));
  assert(to!string(recombination_fractions(msx)[0])=="0.285633");
  // expand map

  /*

  >> expanded_map = map.expand(2.5)
  >> expanded_map.size
  => 1104.28582
  >> expanded_map.positions('1')[0..2]
  => [0, 0.99675, 2.5]
  >> expanded_map.recombination_fractions('1').size
  => 49 
  >> expanded_map.recombination_fractions('1')[0..2]
  => [0.0098688058543826, 0.0148087716805754, 0.024385287749643]

=end


$: << '..'
require 'qtl/rqtl'
require 'test/unit'

TESTDATADIR = '../../../test/data/qtl'
LISTERIAFN  = TESTDATADIR+'/listeria.csv'

# require 'test/unit/testcase'

class TestBiolibRQtl < Test::Unit::TestCase

  def setup
    @qtl = QTL.new(LISTERIAFN)
  end

  def test_info
    d = @qtl.data
    assert_equal(:f2,d.type)
    assert_equal(120,d.individuals.size)
    assert_equal(120,d.nind)
    assert_equal(1,d.nphe)
    assert_equal(133,d.totmar)
    assert_equal(20,d.nchr)
    assert_equal([["1", 13], ["10", 5], ["11", 6], ["12", 6], ["13", 12], ["14", 4], ["15", 8], ["16", 4], ["17", 4], ["18", 4], ["19", 4], ["2", 6], ["3", 6], ["4", 4], ["5", 13], ["6", 13], ["7", 6], ["8", 6], ["9", 7], ["X", 2]],d.nmar.sort)
  end

  def test_markers
    d = @qtl.data
    assert_equal('D10M44',d.marker(0).name)
    assert_equal('2',d.marker(14).chromosome)
    assert_equal(133,d.markers.size)
  end

=begin

Chromosome info

  >> d.chromosomes.size
  => 20
  >> d.chromosomes.autosomes.size
  => 19

=end

  def test_chromosomes
    d = @qtl.data
    assert_equal(20,d.chromosomes.size)
    assert_equal(19,d.chromosomes.autosomes.size)
    # assert_equal('X',d.chromosomes.x.name)
    # assert_equal(13,d.chromosomes[1].markers.size)
  end

=begin

Phenotype column info

  >> d.phenotypenames[0].name
  => 'T264'
  >> d.phenotypenames.size
  => 1

=end

  def test_phenotypecolumns
    d = @qtl.data
    assert_equal('T264',d.phenotypenames[0].name)
    assert_equal(1,d.phenotypenames.size)
  end

=begin

Get statistics

  >> d.perc_phenotyped
  => [96.7]
  >> d.perc_genotyped
  => 88.5

=end

  def test_phenotyped
    d = @qtl.data
    assert_equal([96.7],d.perc_phenotyped)
  end

  def test_genotyped
    d = @qtl.data
    assert_equal(88.5,d.perc_genotyped)
  end

=begin

Get phenotype information for the listeria set (only one 
defined here.

  >> d.phenotype(0,0)
  => 118.317
 
  >> d.phenotype(1,0)
  => 264

Get all phenotypes for the first individual - note this set has only
one phenotype

  >> d.phenotypes[0]
  => [118.317]

Get all phenotypes for all individuals

  >> d.phenotypes[0..4]
  => [[118.317], [264], [194.917], [264], [145.417]]

=end

=begin
R/qtl's genotype matrix for the listeria set contains

  T264  D10M44     D1M3     D1M7
             1        1         1
             0  0.99675  24.84773
  118.317    B        B         B
  264        -        B         B
  194.917    -        H         H
  264        B        B         H

also we test on:

        D1M291    D1M209    D1M155
      84.93474  92.68394  93.64344
             H         H         H
             H         H         H
             B         B         B
             B         B         B
             H         H         H
             B         B         B
             H         H         H
             H         H         H
             H         -         H

which translates internally to:

        D1M291    D1M209    D1M155
             2         2         2
             2         2         2
             3         3         3
             3         3         3
             2         2         2
             3         3         3
             2         2         2
             2         2         2
             2        NA         2

Raw names read from input file

  >> d.genotypes.namesread.sort
  => [["-", 1840], ["A", 3701], ["B", 3387], ["C", 128], ["H", 6904]]

Validated names

  >> d.genotypes.names.sort
  => ["A", "B", "C", "H"]

  >> d.genotypes.na.sort
  => ["-", "NA"]

  >> d.genotypes.alleles.sort
  => ["A", "B"]

First we get the original values from the Listeria .csv file
(individual/marker):

  >> d.individuals[1].genotypes[2].value
  => "B"

  >> d.markers['D1M291'].mid
  => 10

  >> d.individuals[1].genotypes[d.markers['D1M291'].mid].value
  => 'H'

The same, but nicer, query by individual/marker 

  >> d.genotype(0,0)
  => 'B'

  >> d.genotype(1,0)
  => '-'

  >> d.genotype(2,0)
  => '-'

  >> d.genotype(3,0)
  => 'B'

  >> d.genotype(1,2)
  => 'B'

  >> d.genotype(1,'D1M291')
  => 'H'

  >> d.genotype(2,'D2M493')
  => '-'

  >> d.genotype(0,'D13M106')
  => 'A'

  >> d.genotype(2,'D13M106')
  => 'H'

Here we create an adapter for translating genotype information into an
input object suitable for use by R/qtl. Normally a user won't do this, 
as a qtl object is handled transparently.

  >> r = RQtlInputAdaptor.new(d)

  >> r.genotype(1,'D1M291')
  => 2

  >> r.genotype(8,'D1M209')
  => 'NA'

  >> r.genotype(3,'D1M155')
  => 3

  >> r.genotype(8,'D10M44')
  => 1

Scanone expects a flat Array where the first five elements represent the 
first marker genotypes for 5 individuals:

  >> r1 = RQtlScanoneAdaptor.new(d)
  >> r1.use_individuals.size
  => 116
  >> r1.scanone_ingenotypematrix[0..4]
  => [3,0,0,3,2]

Note this Ruby implementation of R/qtl really retains the original dataset for
querying and uses adaptors for querying modified or 'derived' datasets. The
RQtlInputAdaptor above converts genotype names from:

  >> d.genotypes.names
  => ["A", "B", "C", "H"]

to the indexed values used by R/qtl internally:

  >> r.genotypes.names
  => [1,3,4,2]

  >> r.genotypes.alleles
  => [1,3]

  >> r.genotypes.na
  => ["NA","NA"]

  >> r.phenotypes[0..4]
  => [[118.317], [264], [194.917], [264], [145.417]]

Similarly the RQtlScanoneAdaptor modifies the input dataset to make it suitable
for single QTL mapping with R/qtl. For example it drops the 4 'NA' phenotypes.
Note, again, the adaptor is normally not seen by the end user.

=end

  def test_rqtl_input_adaptor
    d = @qtl.data
    assert_equal('B',d.genotype(0,0))
    r = RQtlInputAdaptor.new(d)
    assert_equal(3,r.genotype(0,0))
    assert_equal(2,r.genotype(1,'D1M291'))
    assert_equal([1,3,4,2],r.genotypes.names)
  end


  def test_rqtl_scanone_adaptor
    d = @qtl.data
    assert_equal(["NA"],d.phenotypes[29])
    assert_equal(120,d.phenotypes.size)
    r1 = RQtlScanoneAdaptor.new(d)
    # p r1.use_individuals
    assert_equal(116,r1.use_individuals.size)
    assert_equal([3,0,0,3,2],r1.scanone_ingenotypematrix[0..4])
    assert_equal([1, 2, 2, 2, 2, 1, 1, 1, 2, 2, 1],r1.scanone_ingenotypematrix[15410..15420])
  end

=begin

We are going to scan for QTL's. The first R equivalent here is:

    > mr = scanone(listeria,method='mr')
      Warning message:
      Dropping 4 individuals with missing phenotypes.
        in: checkcovar(cross, pheno.col, addcovar, intcovar, perm.strata,  

    > mr
              chr      pos         lod
      D10M44    1  0.00000 0.457256443
      D1M3      1  0.99675 0.687377692
      D1M75     1 24.84773 0.244500303
      D1M215    1 40.41361 0.070120158
      (...)
      D19M117  19 16.36398 0.436675089
      D19M65   19 32.82935 0.007363598
      D19M10   19 44.49432 0.000000000
      DXM186    X  0.00000 0.678435461
      DXM64     X 42.34593 0.002756074

Execute single QTL mapping using R/qtl's marker regression passing in
the @qtl dataset. We need the RQTL convenience class to map against biolib:

    >> require 'qtl/rqtl'
    >> rqtl = RQTL.new(qtl)

Now execute single QTL mapping with scanone using multiple regression analysis

    >> mr = rqtl.scanone_mr()

Which returns and array of markers containing the marker name, chromosome,
position and lod score:

    >> mr[0].lod
    => 0.457256442704442
    >> mr[0].name
    => 'D10M44'
    >> mr[0].chromosome
    => '1'
    >> mr[0].position
    => 0
    >> mr[0].to_a
    => ["D10M44", "1", 0, 0.457256442704442]

    >> mr[3].lod
    => 0.0701201576892285
    >> mr[3].name
    => 'D1M215'
    >> mr[3].chromosome
    => '1'
    >> mr[3].position
    => 40.41361 

Now for the multiple imputation method by Sen and Churchill (2001). In R:

    > data(listeria)
    > gp =  sim.geno(listeria,step=2.5,n.draws=8)
    > out.imp <- scanone(gp, method="imp")
    > out.imp[74,]
                 chr  pos       lod
      c2.loc52.5   2 52.5 0.9199136
    > out.imp[555,]
    > cX.loc17.5   X 17.5 0.3018425

The Ruby equivalent is (FIXME: this is not working for some reason)

    >> rqtl.expand_markers!(2.5)
    !>> mr = rqtl.scanone_imp(rqtl.sim_geno(8))
    !>> mr[74].to_a
    !=> ["c2.loc52.5", "2", 52.2, 0.9199136]
    !>> mr[555].to_a
    !=> []

=end

  def test_scanone
    rqtl = RQTL.new(@qtl)
    mr = rqtl.scanone_mr()
    # roundoff, otherwise complains:
    assert_equal(0.457,(mr[0].lod*1000).to_i/1000.0)
  end
end
*/

}
