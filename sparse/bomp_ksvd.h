#ifndef BOMP_KSVD_H
#define BOMP_KSVD_H

// Input:
// dict       = dictionary with nbAtoms x atomSize atoms (rows x columns), in row-major form
// signals    = array of nbSignals x atomSize signals, in row-major form
// maxEntries = maximum number of atoms for reconstruction
// error      = error threshold for reconstruction
//
// Output:
// nbEntries  = number of atom indices and weights per signal
// indices    = array of maxEntries indices per signal
// values     = array of maxEntries weights per signal

void BOMP(
    const double* dict, int atomSize, int nbAtoms,
    const double* signals, int nbSignals,
    int maxEntries, double error,
    int* nbEntries, int* indices, double* values);

void KSVD(
    double* dict, int atomSize, int nbAtoms,
    const double* signals, int nbSignals,
    int maxEntries, double error,
    int nbIters);

#endif