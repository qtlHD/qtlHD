# Mapping inputs

In general, file standards for QTL mapping are a muddle. The basic premises,
however, are pretty straightforward. Here, for qtlHD we propose a 'final'
tab-delimited format that will also have a binary (and perhaps XML) counterpart
representation. A choice for tab-delimited (or TSV, tab separated value) files
is made to secure human readability (and editability!) of data. Since TSV files
contain no 'grammar' we describe content and set options through section header
information. In principle we attempt to DRY (do not repeat yourself) and sparse
data (no inclusion of zero values). The idea is to have a simple
straightforward standard that can be applied to 95% of cases, is easily
readable by humans, and easily parsable by software.

In addition, the new tab delimited format should be versioned, extensible,
handle names with spaces (except tabs), and allow some validation of section
sizes. Section data is either TSV, or JSON. JSON is used for tree like data
structures.  Finally, we support both multi-file and single-file
representations of the same data. The single file represenstation simply being
a concatenation of other files.

Note: The qtlHD project is the single authority for the file standard. qtlHD
acts as a reference implementation.

Each file is described individually in

* Symbol file        - symbols used to describe (observed) genotypes, 
                       markers, founders, individuals and phenotypes
* Founder file       - Founder and cross layout
* Genotype file      - marker x individual genotyping
* Marker map file    - marker (sex) chromosome positions
* Phenotype file     - individual x phenotypes

Before describing the formats we describe parser logic, and the general
file/section header layout.

# Parser logic

Any line that starts with a Hash # symbol is considered a comment. Meta-data
for the parser starts with # --- Command - so hash and three dashes.

A comment may start also later in a line. This is signified by multiple
spaces, a hash and a space. Not tabs after are allowed. Not all sections
will allow these types of comments.

# Header layout

A header should be easily recognizable, and give a version. Every file 
or file section should start with a commented line (using spaces)

      # --- qtlHD-in-x.x Sectionname Description

where x.x is the file format version and name represents the section name. The 
description is a free string. Valid section headers are, for example

      # --- qtlHD-in-0.9 Genotype Mouse CC JAX
      # --- qtlHD-in-0.9 Phenotype Mouse CC JAX - Our experiment/date

This line is used to recognise and start parsing sections. After the section
header line, a header is recognized by 

      # --- Header begin

and ends with

      # --- Header end

Valid data sections are defined by each format. Typically a data section starts
and ends with

      # --- Data Name begin
      # --- Data Name end

This is an extensible data format. To avoid getting all different formats,
again, try to fit your data model in an existing format, before adding a new
one. And note again that qtlHD is the reference implementation.

# The symbol file

The symbol file substitutes symbols with values. A block should start/end with

      # --- Data Type begin/end

A possible founder symbol table

      # --- Data Founder begin
      A/J 0
      C57BL/6J 1
      129S1/SvImJ 2
      ...
      # --- Data Founder end

When leaving out the numbers, they will be created automatically. So above
is equal to:

      # --- Data Founder begin
      A/J
      C57BL/6J
      129S1/SvImJ
      ...
      # --- Data Founder end

Likewise, individuals, phenotypes, markers and chromosomes can have symbol to
value substitutions.

Observed genotypes have a somewhat more complex symbol to value substitution.
To translate a genotype symbol to multiple possible (observed) true genotypes,
a possible genotype symbol table:

      # --- Data Genotype begin
      NA - as None        
      A as 0,0
      B BB as 1,1
      C CC as 2,2         # aliases
      AB as 1,0           # directional
      BA as 0,1
      AC CA as 0,2 2,0    # not directional
      D as 3              # expands to 3,3
      E as A/J            # Founder symbol expands to 0,0
      ...
      # --- Data Genotype end

A symbol is always a string, and can contain spaces (but no tabs). Note that
remarks are also allowed, starting with # and a space.  The 'as' keyword splits
aliased symbols and multiple genotypes (make sure there is no symbol named
'as'!). The genotype numbers refer to the founders. If founders have symbols,
these symbols can be used too. Symbols are simply expanded to their numeric
values.

# The founder file

A founder file contains information about the founders useful for algorithms.
Possibilities are cross description and allele frequencies.

# The genotype file

The genotype file contains markers (columns) x individuals (rows), giving the
(observed) genotypes. The values can either be numbers, or symbols (but not
both). Valid genotypes could be (tab delimited):

        # --- Data Observed begin
        #     M1   M2    M3   M4    M5      M6   
        Ind1   A   AorB   1   1,1   0,1   0,1|1,1
        Ind2   H   AorH   NA  1,1   1,1     1,1
        ...
        # --- Data Observed end

where the last marker of Ind1 represents, for example, an AorH. Note the use of
NA - which is the standard name in qtlHD for a missing value. Also note the
treatment of marker names. Marker names are normally listed in the symbol
table, or are deduced from the marker map file. The names in the columns are
normally ignored (as they start with a Hash #). 

# The marker map file

The marker map contains a list of markers and their chromosome + locations.
E.g.

        # --- Data Location begin
        #   Chr   Pos
        M1   I   100.0
        M2   I   2000.0
        M3   X   100.0
        M4   X   1500.0
        ...
        # --- Data Location end

where, again, markernames and chromosomes can be symbols.

# The phenotype file

The phenotype file contains a list of individuals and their phenotypic values.
E.g.

        # --- Data Phenotypes begin
        #    Sex   P1      P2
        Ind1   X   100.0   45
        Ind2   I   2000.0  46
        ...
        # --- Data Phenotypes end

Where the phenotype names are listed in the symbols file. The optional 
name row is ignored by qtlHD as it starts with a Hash symbol.

Phenotypic values can have types. E.g. floating point (which is the default
type), integer, binary, discrete, and perhaps ranges. These are defined in the
header, using the phenotype names:

        # --- Type Phenotypes begin
        Sex Discrete                        # types are automatically found from data
        P1  Float                           # double precision value with floating point (default)
        P2  Integer
        P3  Ranges 0..1 1..4 4..112 112..$  # valid ranges
        P4  Discrete  T,N,S                 # predefined discrete types
        P5  Integer 1..$                    # only values larger than 1
        ...
        # --- Type Phenotypes end

Note that types are not symbols. A symbol, defined in the symbol file,
represents a simple value expansion. A type says someting about possible
values. When a type is illegal (say a float for an int), or a type falls
outside a predefined set, which can be checked with P2, P3, P4, and P5, the
software should throw an error.

