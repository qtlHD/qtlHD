/**
 * Utility functions for permutation tests
 **/

module qtl.core.util.permutation;

import std.stdio;
import std.math;
import std.random;

// returns shuffled version of vector (0, 2, ..., n-1)
size_t[] random_permutation(in size_t n, ref Random gen)
{
  auto output = new size_t[](n);
  foreach(i, ref o; output) o = i;

  randomShuffle(output, gen);
  return output;
}

unittest {
  writeln("Unit test " ~ __FILE__);

  Random gen;
  gen.seed(unpredictableSeed);

  write("Randomize (0, 1, ..., 9): ");
  auto x = random_permutation(10, gen);
  foreach(o; x)
    writef("%d ", o);
  writeln;
}

// returns matrix, with each column being a shuffled version of vector (0, 2, ..., n-1)
size_t[][] random_permutation(in size_t n_perm, in size_t n_ind, ref Random gen)
{
  auto output = new size_t[][](n_perm, n_ind);

  foreach(i; 0..n_perm)
    output[i] = random_permutation(n_ind, gen);

  return output;
}

unittest {
  Random gen;
  gen.seed(unpredictableSeed);

  writeln("\nFour random permutations of (0, 1, ..., 9):");
  auto x = random_permutation(4, 10, gen);
  foreach(i; 0..4) {
    write("    ");
    foreach(o; x[i])
      writef("%d ", o);
    writeln;
  }
}

// stratified permutations
// shuffle 0..strata.length within strata defined by strata
size_t[] random_stratified_permutation(in char strata[], ref Random gen)
{
  size_t[size_t][char] stratav;
  size_t[char] stratav_length;
  auto output = new size_t[](strata.length);

  // create strata indices
  foreach(i, s; strata) {
    stratav_length[s]++;
    stratav[s][stratav_length[s]-1] = i;
  }

  // randomize within each stratum
  foreach(s; stratav.keys) {
    auto x = random_permutation(stratav[s].length, gen);
    foreach(i, si; stratav[s])
      output[si] = stratav[s][x[i]];
  }

  return output;
}

unittest {
  Random gen;
  gen.seed(unpredictableSeed);

  auto strata = ['A', 'A', 'B', 'B', 'A', 'A', 'B', 'B'];
  write("\nstrata: ");
  foreach(s; strata)
    writef("%s ", s);
  writeln;

  auto x = random_stratified_permutation(strata, gen);
  write("perm:   ");
  foreach(xx; x)
    writef("%d ", xx);
  writeln;

  // check that the perms were within strata
  auto y = strata.dup;
  foreach(i; 0..x.length)
    y[i] = strata[x[i]];

  assert(y == strata);


  strata = ['A', 'A', 'A', 'A', 'B', 'B', 'B', 'B', 'C', 'C', 'C', 'C'];
  write("\nstrata: ");
  foreach(s; strata)
    writef(" %s ", s);
  writeln;

  x = random_stratified_permutation(strata, gen);
  write("perm:   ");
  foreach(xx; x)
    writef("%2d ", xx);
  writeln;

  // check that the perms were within strata
  y = strata.dup;
  foreach(i; 0..x.length)
    y[i] = strata[x[i]];

  assert(y == strata);
}

// stratified permutations, multiply
// shuffle 0..strata.length within strata defined by strata
size_t[][] random_stratified_permutation(in size_t n_perm, in char strata[], ref Random gen)
{
  size_t[size_t][char] stratav;
  size_t[char] stratav_length;
  auto output = new size_t[][](n_perm, strata.length);

  // create strata indices
  foreach(i, s; strata) {
    stratav_length[s]++;
    stratav[s][stratav_length[s]-1] = i;
  }

  // randomize within each stratum
  foreach(s; stratav.keys) {
    auto x = random_permutation(n_perm, stratav[s].length, gen);
    foreach(perm; 0..n_perm) {
      foreach(i, si; stratav[s])
        output[perm][si] = stratav[s][x[perm][i]];
    }
  }

  return output;
}

unittest {
  Random gen;
  gen.seed(unpredictableSeed);

  auto strata = ['A', 'A', 'A', 'A', 'B', 'B', 'B', 'B', 'C', 'C', 'C', 'C'];
  write("\nstrata: ");
  foreach(s; strata)
    writef(" %s ", s);
  writeln;

  auto x = random_stratified_permutation(4, strata, gen);
  foreach(i; 0..4) {
    writef("perm %d: ", i+1);
    foreach(xx; x[i])
      writef("%2d ", xx);
    writeln;
  }

  // check that the perms were within strata
  auto y = strata.dup;
  foreach(i; 0..4) {
    foreach(j; 0..x[i].length)
      y[j] = strata[x[i][j]];
    assert(y == strata);
  }
}