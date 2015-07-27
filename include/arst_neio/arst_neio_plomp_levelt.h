/** @file */
#ifndef HEADER_PLOMP_LEVELT_H_INCLUDED
#define HEADER_PLOMP_LEVELT_H_INCLUDED

#include "arst_neio/arst_neio_gene.h"

/// Returns the Plomp-Levelt dissonance of a pair of Harmo
double disso(const Harmo& h1, const Harmo& h2, const Generic& cg);
/// Returns the Plomp-Levelt dissonance of a pair of Pitch (sum on their pairs of Harmo)
double disso(const Pitch& p1, const Pitch& p2, const Generic& cg);
/// Returns the Plomp-Levelt dissonance of a pair of Chord (sum on their pairs of Harmo)
double disso(const Chord& c1, const Chord& c2, const Generic& cg);
/// Returns the Plomp-Levelt dissonance of a pair of Chords (sum on their pairs of Harmo)
double disso(const ChordS& t1, const ChordS& t2, const Generic& cg); 
/// Returns the Plomp-Levelt dissonance of a pair of Scale (sum on their pairs of Harmo)
double disso(const Scale& s1, const Scale& s2, const Generic& cg); 

/// Returns the Plomp-Levelt consonance of a pair of Harmo (-disso up to a const)
double conso(Harmo h1, Harmo h2, Generic cg); //double d1, double w1, double d2, double w2, double w);

/// Optimizes the Plompt-Levelt dissonance
void plomp_levelt(const std::string& o_f_n);

#endif
