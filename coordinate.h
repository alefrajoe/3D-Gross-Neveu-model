#ifndef COORDINATE_H
#define COORDINATE_H

#include "macro.h"

// coordinate stuff
void PrintCoordinate(const int * coor_array);
void ConvertLessicographicIntoCoordinate(const int number, int * coor_array);
void ConvertCoordinateIntoLessicographic(const int const * array, int * number);
void InitializePPArray(int site, int * destination);
void InitializeMMArray(int site, int * destination);
void StartEtapp(int site, double * destination);
void StartEtamm(int site, double * destination);
void StartEtapbc(int site, double * destination);
void StartEpsilonPi(int site, int fermion_site, double * destination);

// coordinate fermions
inline int ParityCoordinate(const int const * array){int parity = 0; for(int d=0; d<DIM; d++) parity += array[d]; return (parity % 2);};
void ConvertFermionLessicographicIntoLessicographic(const int fermion, int * stereo);
void ConvertFermionLessicographicIntoCoordinate(const int number, int * coor_array);
void ConvertLessicographicIntoFermionLessicographic(const int stereo, int * fermion);
void ConvertCoordinateIntoFermionLessicographic(const int const * array, int * number);
void InitializeFermionPPArray(int site, int * destination);
void InitializeFermionMMArray(int site, int * destination);
void InitializeFastFermionPPArray(int sitepp, int * destination, const int const * fermionsite);
#endif