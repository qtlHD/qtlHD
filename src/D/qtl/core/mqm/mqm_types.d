/**
 * Genotype module (deprecated)
 */

module qtl.core.mqm.mqm_types;

import std.conv;
import std.stdio;
import qtl.core.primitives;

/**
 * Enum (C) mqm genotypes
 *
 */
 
enum MQMCrossType { CUNKNOWN = 'U', CF2 = 'F', CBC = 'B', CRIL = 'R' };
enum MQMCofactorType { MNOCOF = '0', MCOF ='1', MSEX = '2', MQTL = '3' };
enum MQMMarker { MAA = '0', MH = '1', MBB = '2', MNOTAA    = '3', MNOTBB = '4',  MMISSING = '9', MUNUSED = '-'};
