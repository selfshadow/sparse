#include "bomp_ksvd.h"

#include <memory.h>
#include <math.h>


inline double Dot(const double* a, const double* b, int len)
{
    double r = 0.0;
    for (int i = 0; i < len; i++)
        r += a[i]*b[i];
    return r;
}

// Forward substitution:
// Solve M.x = y for x, with lower-triangular matrix M
inline void ForeSub(const double* M, const double* y, double* x, int len, int dim)
{
    for (int i = 0; i < len; i++)
        x[i] = (y[i] - Dot(&M[i*dim], x, i))/M[i + i*dim];
}

// Back substitution, with transpose:
// Solve (M^T).x = y for x, with lower-triangular matrix M
inline void BackSubT(const double* M, const double* y, double* x, int len, int dim)
{
    for (int i = len - 1; i >= 0; i--)
    {
        double dot = 0.0;
        for (int j = i + 1; j <= len - 1; j++)
            dot += M[i + j*dim]*x[j];
        x[i] = (y[i] - dot)/M[i + i*dim];
    }
}

// Returns the index to the entry in v with the highest absolute value
inline int FindAbsMax(const double* v, int len)
{
    int    max_i = 0;
    double max_x = 0.0;
    for (int i = 0; i < len; i++)
    {
        double x = fabs(v[i]);
        if (max_x < x)
        {
            max_x = x;
            max_i = i;
        }
    }
    return max_i;
}

inline void GramMatrix(const double* M, double* G, int nbColumns, int nbRows)
{
    // TODO: Just compute the upper triangle and replicate?
    for (int i = 0; i < nbRows; i++)
    for (int j = 0; j < nbRows; j++)
    {
        float dot = 0.0;
        for (int k = 0; k < nbColumns; k++)
            dot += M[i*nbColumns + k]*M[j*nbColumns + k];
        G[i*nbRows + j] = dot;
    }
}

inline void Normalize(double* v, int len)
{
    double dot = 0.0;
    for (int i = 0; i < len; i++)
        dot += v[i]*v[i];

    // TODO: Needed?
    if (dot == 0.0)
        return;

    double m = 1.0/sqrt(dot);
    for (int i = 0; i < len; i++)
        v[i] *= m;
}


struct OMPContext
{
    OMPContext(int atomSize_, int nbAtoms_, int maxEntries_) :
        atomSize(atomSize_),
        nbAtoms(nbAtoms_),
        maxEntries(maxEntries_)
    {
        G      = new double[nbAtoms*nbAtoms];
        alpha0 = new double[nbAtoms];
        alpha  = new double[nbAtoms];
        used   = new bool  [nbAtoms];
        L      = new double[maxEntries*maxEntries];

        temp0  = new double[maxEntries];
        temp1  = new double[maxEntries];
    };

    ~OMPContext()
    {
        delete[] G;
        delete[] alpha0;
        delete[] alpha;
        delete[] used;
        delete[] L;

        delete[] temp0;
        delete[] temp1;
    }

    int atomSize;
    int nbAtoms;
    int maxEntries;

    double* G;
    double* alpha0;
    double* alpha;
    bool*   used;
    double* L;

    double* temp0;
    double* temp1;
};


void BOMPInternal(
    OMPContext& c,
    const double* dict, const double* signals, int nbSignals,
    int* entries, int* indices, double* values)
{
    GramMatrix(dict, c.G, c.atomSize, c.nbAtoms);

    for (int s = 0; s < nbSignals; s++)
    {
        for (int i = 0; i < c.nbAtoms; i++)
        {
            // Compute alpha0 = D.x
            c.alpha0[i] = Dot(&dict[i*c.atomSize], &signals[s*c.atomSize], c.atomSize);
            c.alpha[i]  = c.alpha0[i];
            c.used[i]   = false;
        }

        c.L[0] = 1.0;

        int nbEntries = 0;

        for (int i = 0; i < c.maxEntries; i++)
        {
            int aIndex = FindAbsMax(c.alpha, c.nbAtoms);

            // Exit if atom has already been used
            if (c.used[aIndex])
                break;

            c.used[aIndex] = true;

            if (nbEntries > 0)
            {
                for (int j = 0; j < nbEntries; j++)
                    c.temp0[j] = c.G[aIndex*c.nbAtoms + indices[j]];

                ForeSub(c.L, c.temp0, c.temp1, nbEntries, c.maxEntries);

                double* Lrow = &c.L[nbEntries*c.maxEntries];
                for (int j = 0; j < nbEntries; j++)
                    *Lrow++ = c.temp1[j];
                *Lrow = sqrt(1 - Dot(c.temp1, c.temp1, nbEntries));
            }

            indices[nbEntries++] = aIndex;

            for (int j = 0; j < nbEntries; j++)
                c.temp0[j] = c.alpha0[indices[j]];

            ForeSub (c.L, c.temp0, c.temp1, nbEntries, c.maxEntries); // TODO: This step can be done incrementally
            BackSubT(c.L, c.temp1, values,  nbEntries, c.maxEntries);

            for (int j = 0; j < c.nbAtoms; j++)
            {
                double dot = 0.0;
                for (int m = 0; m < nbEntries; m++)
                    dot += c.G[j*c.nbAtoms + indices[m]]*values[m];
                c.alpha[j] = c.alpha0[j] - dot;
            }
        }

        *entries++ = nbEntries;
        indices   += c.maxEntries;
        values    += c.maxEntries;
    }
}


struct KSVDContext
{
    KSVDContext(int atomSize_, int nbAtoms_, int maxEntries_, int nbSignals_) :
        atomSize(atomSize_),
        nbAtoms(nbAtoms_),
        maxEntries(maxEntries_),
        nbSignals(nbSignals_),
        nbRefs(0)
    {
        GammaE = new int[nbSignals];
        GammaI = new int[nbSignals*maxEntries];
        GammaV = new double[nbSignals*maxEntries];
        GammaJ = new double[nbSignals];

        refs   = new int[nbSignals];

        atom   = new double[atomSize];
    }

    ~KSVDContext()
    {
        delete[] GammaE;
        delete[] GammaI;
        delete[] GammaV;
        delete[] GammaJ;

        delete[] refs;

        delete[] atom;
    }

    int atomSize;
    int nbAtoms;
    int nbSignals;
    int maxEntries;

    int*    GammaE; // Entry counts (one per signal)
    int*    GammaI; // Atom indices (up to maxEntries per signal)
    double* GammaV; // Atom weights (up to maxEntries per signal)
    double* GammaJ; // Atom weights for the current (J-th) atom

    int     nbRefs;
    int*    refs;   // Indices for the signals referencing the current atom

    double* atom;   // Revised atom
};


void FindReferences(KSVDContext& k, int j)
{
    k.nbRefs = 0;

    for (int i = 0; i < k.nbSignals; i++)
    {
        // TODO: Avoid this search
        for (int m = 0; m < k.GammaE[i]; m++)
        {
            int atom = k.GammaI[i*k.maxEntries + m];
            if (atom == j)
            {
                k.refs[k.nbRefs] = i;
                k.GammaJ[k.nbRefs] = k.GammaV[i*k.maxEntries + m];
                k.nbRefs++;
            }
        }
    }
}

void UpdateAtom(KSVDContext& k, const double* dict, const double* signals)
{
    for (int r = 0; r < k.nbRefs; r++)
    {
        int sIndex = k.refs[r];

        for (int i = 0; i < k.atomSize; i++)
        {
            double x = 0.0;
            for (int m = 0; m < k.GammaE[sIndex]; m++)
            {
                int aIndex = k.GammaI[sIndex*k.maxEntries + m];
                x += dict[aIndex*k.atomSize + i]*k.GammaV[sIndex*k.maxEntries + m];
            }

            k.atom[i] += k.GammaJ[r]*(signals[sIndex*k.atomSize + i] - x);
        }
    }
}

void UpdateWeights(KSVDContext& k, const double* dict, const double* signals)
{
    for (int r = 0; r < k.nbRefs; r++)
    {
        int sIndex = k.refs[r];

        for (int i = 0; i < k.atomSize; i++)
        {
            double x = 0.0;
            for (int m = 0; m < k.GammaE[sIndex]; m++)
            {
                int aIndex = k.GammaI[sIndex*k.maxEntries + m];
                x += dict[aIndex*k.atomSize + i]*k.GammaV[sIndex*k.maxEntries + m];
            }

            k.GammaJ[r] += k.atom[i]*(signals[sIndex*k.atomSize + i] - x);
        }
    }
}


void KSVDStep(
    KSVDContext& k, OMPContext& c,
    double* dict, const double* signals)
{
    BOMPInternal(c, dict, signals, k.nbSignals, k.GammaE, k.GammaI, k.GammaV);

    // TODO: Use random atom permutation
    for (int j = 0; j < k.nbAtoms; j++)
    {
        // Find signals referencing the J-th atom
        FindReferences(k, j);

        // Clear the atom
        for (int i = 0; i < k.atomSize; i++)
            dict[j*k.atomSize + i] = 0.0;

        // TODO: replace atom when unreferenced
        if (k.nbRefs == 0)
            continue;

        // Compute revised atom
        for (int i = 0; i < k.atomSize; i++)
            k.atom[i] = 0.0;
        UpdateAtom(k, dict, signals);  // TODO: Cache computation for UpdateWeights later
        Normalize(k.atom, k.atomSize);

        // Compute new weights
        for (int i = 0; i < k.nbRefs; i++)
            k.GammaJ[i] = 0.0;
        UpdateWeights(k, dict, signals);

        // Update dictionary with atom
        for (int i = 0; i < k.atomSize; i++)
            dict[j*k.atomSize + i] = k.atom[i];

        // Update signal weights
        for (int i = 0; i < k.nbRefs; i++)
        {
            int sIndex = k.refs[i];

            // TODO: Avoid this search
            for (int m = 0; m < k.GammaE[sIndex]; m++)
            {
                int aIndex = k.GammaI[sIndex*k.maxEntries + m];
                if (aIndex == j)
                    k.GammaV[sIndex*k.maxEntries + m] = k.GammaJ[i];
            }
        }
    }
}


void BOMP(
    const double* dict, int atomSize, int nbAtoms,
    const double* signals, int nbSignals, int maxEntries,
    int* nbEntries, int* indices, double* values)
{
    OMPContext c(atomSize, nbAtoms, maxEntries);
    BOMPInternal(c, dict, signals, nbSignals, nbEntries, indices, values);
}


void KSVD(
    double* dict, int atomSize, int nbAtoms,
    const double* signals, int nbSignals,
    int maxEntries,
    int nbIters)
{
    int dictSize = atomSize*nbAtoms;

    KSVDContext k(atomSize, nbAtoms, maxEntries, nbSignals);
    OMPContext  c(atomSize, nbAtoms, maxEntries);

    for (int i = 0; i < nbIters; i++)
        KSVDStep(k, c, dict, signals);
}