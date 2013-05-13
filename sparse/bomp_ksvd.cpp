#include "bomp_ksvd.h"

#include <math.h>


inline void Set(double* v, double x, int len)
{
    for (int i = 0; i < len; i++)
        v[i] = x;
}

inline void Copy(double* r, const double* v, int len)
{
    for (int i = 0; i < len; i++)
        r[i] = v[i];
}

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
    int    maxI = 0;
    double maxX = 0.0;

    for (int i = 0; i < len; i++)
    {
        double x = fabs(v[i]);
        if (maxX < x)
        {
            maxX = x;
            maxI = i;
        }
    }

    return maxI;
}

inline void GramMatrix(const double* M, double* G, int nbColumns, int nbRows)
{
    // TODO: Just compute the upper triangle and replicate?
    for (int i = 0; i < nbRows; i++)
    for (int j = 0; j < nbRows; j++)
    {
        double dot = 0.0;
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
    OMPContext(int atomSize_, int nbAtoms_, int maxEntries_, double epsilon_) :
        atomSize(atomSize_),
        nbAtoms(nbAtoms_),
        maxEntries(maxEntries_),
        epsilon(epsilon_*epsilon_)
    {
        G      = new double[nbAtoms*nbAtoms];
        alpha0 = new double[nbAtoms];
        alpha  = new double[nbAtoms];
        beta   = new double[nbAtoms];
        used   = new bool  [nbAtoms];
        L      = new double[maxEntries*maxEntries];
        temp0  = new double[maxEntries];
        temp1  = new double[maxEntries];

        // Set a sensible lower bound to avoid numerical instability
        if (epsilon < 1e-10)
            epsilon = 1e-10;
    };

    ~OMPContext()
    {
        delete[] G;
        delete[] alpha0;
        delete[] alpha;
        delete[] beta;
        delete[] used;
        delete[] L;
        delete[] temp0;
        delete[] temp1;
    }

    int atomSize;
    int nbAtoms;

    int maxEntries;
    double epsilon;

    double* G;
    double* alpha0;
    double* alpha;
    double* beta;
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

    const double* signal = signals;

    for (int s = 0; s < nbSignals; s++)
    {
        for (int i = 0; i < c.nbAtoms; i++)
        {
            // Compute alpha0 = D.x
            c.alpha0[i] = Dot(&dict[i*c.atomSize], signal, c.atomSize);
            c.alpha[i]  = c.alpha0[i];
            c.used[i]   = false;
        }

        c.L[0] = 1.0;

        double eps = Dot(signal, signal, c.atomSize);
        double delta = 0.0;

        int nbEntries = 0;

        for (int i = 0; i < c.maxEntries; i++)
        {
            // Exit if we've hit the threshold
            if (eps <= c.epsilon)
                break;

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
                c.beta[j] = 0.0;
                for (int m = 0; m < nbEntries; m++)
                    c.beta[j] += c.G[j*c.nbAtoms + indices[m]]*values[m];
                c.alpha[j] = c.alpha0[j] - c.beta[j];
            }

            eps += delta;
            delta = 0.0;
            for (int m = 0; m < nbEntries; m++)
                delta += c.beta[indices[m]]*values[m];
            eps -= delta;
        }

        *entries++ = nbEntries;
        indices   += c.maxEntries;
        values    += c.maxEntries;
        signal    += c.atomSize;
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
        used   = new bool[nbSignals];
        atom   = new double[atomSize];

        for (int i = 0; i < nbSignals; i++)
            used[i] = false;
    }

    ~KSVDContext()
    {
        delete[] GammaE;
        delete[] GammaI;
        delete[] GammaV;
        delete[] GammaJ;
        delete[] refs;
        delete[] used;
        delete[] atom;
    }

    int atomSize;
    int nbAtoms;
    int nbSignals;
    int maxEntries;

    int nbRefs;

    int*    GammaE; // Entry counts (one per signal)
    int*    GammaI; // Atom indices (up to maxEntries per signal)
    double* GammaV; // Atom weights (up to maxEntries per signal)
    double* GammaJ; // Atom weights for the current (J-th) atom
    int*    refs;   // Indices for the signals referencing the current atom
    bool*   used;   // Track signals that have replaced unused atoms in the dictionary
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

void ReplaceAtom(KSVDContext& k, const double* dict, const double* signals, double* atom)
{
    int    maxI = 0;
    double maxE = 0.0;

    for (int s = 0; s < k.nbSignals; s++)
    {
        double error = 0.0;

        for (int i = 0; i < k.atomSize; i++)
        {
            double x = 0.0;
            for (int m = 0; m < k.GammaE[s]; m++)
            {
                int aIndex = k.GammaI[s*k.maxEntries + m];
                x += dict[aIndex*k.atomSize + i]*k.GammaV[s*k.maxEntries + m];
            }

            double v = signals[s*k.atomSize + i] - x;
            error += v*v;
        }

        if (maxE < error && !k.used[s])
        {
            maxE = error;
            maxI = s;
        }
    }

    // Mark signal so it's not used to replace multiple atoms
    k.used[maxI] = true;

    for (int i = 0; i < k.atomSize; i++)
        atom[i] = signals[maxI*k.atomSize + i];
    Normalize(atom, k.atomSize);
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

        // If the atom isn't referenced, replace it with the worst-represented signal
        if (k.nbRefs == 0)
        {
            ReplaceAtom(k, dict, signals, &dict[j*k.atomSize]);
            continue;
        }

        // Clear the current atom in the dictionary
        Set(&dict[j*k.atomSize], 0.0, k.atomSize);

        // Compute revised atom
        Set(k.atom, 0.0, k.atomSize);
        UpdateAtom(k, dict, signals);  // TODO: Cache computation for UpdateWeights later
        Normalize(k.atom, k.atomSize);

        // Compute new weights
        Set(k.GammaJ, 0.0, k.nbRefs);
        UpdateWeights(k, dict, signals);

        // Update dictionary with atom
        Copy(&dict[j*k.atomSize], k.atom, k.atomSize);

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
    const double* signals, int nbSignals,
    int maxEntries, double epsilon,
    int* nbEntries, int* indices, double* values)
{
    OMPContext c(atomSize, nbAtoms, maxEntries, epsilon);
    BOMPInternal(c, dict, signals, nbSignals, nbEntries, indices, values);
}


void KSVD(
    double* dict, int atomSize, int nbAtoms,
    const double* signals, int nbSignals,
    int maxEntries, double epsilon,
    int nbIters)
{
    int dictSize = atomSize*nbAtoms;

    KSVDContext k(atomSize, nbAtoms, maxEntries, nbSignals);
    OMPContext  c(atomSize, nbAtoms, maxEntries, epsilon);

    for (int i = 0; i < nbIters; i++)
        KSVDStep(k, c, dict, signals);
}