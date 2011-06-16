/**
 * HMM module
 */

module qtl.core.hmm;

// things I think I really need
import qtl.core.primitives;
import qtl.core.genotype;
import qtl.core.cross;
import std.stdio;
import std.math;

// things for the unit tests 
import qtl.plugins.input.read_csv;
import std.path;


// calculate QTL genotype probabilities
double[size_t][size_t][T] calcGenoprob(T)(Cross cross, double[] rec_frac, double error_prob)
{
  
}

unittest {
  writeln("Unit test " ~ __FILE__);
  alias std.path.join join;
  auto fn = dirname(__FILE__) ~ sep ~ join("..","..","..","..","test","data","input","listeria.csv");
  writeln("  - read CSV " ~ fn);
  auto data = new ReadSimpleCSV!F2(fn);
  auto cross = new F2Cross(data.genotypes);
}

double[size_t][Genotype!T] forwardEquations(T)(Cross cross, double[] rec_frac, double error_prob)
{
  double[size_t][Genotype!T] alpha;

}

double[size_t][Genotype!T] backwardEquations(T)(Cross cross, double[] rec_frac, double error_prob)
{
  double[size_t][Genotype!T] beta;

}

