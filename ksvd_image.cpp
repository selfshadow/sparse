#include "sparse/bomp_ksvd.h"

#include "stb/stb_image.h"
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb/stb_image_write.h"

#include <math.h>

typedef unsigned char byte;


double clamp(double v, double min, double max)
{
    if (v > max) return max;
    if (v < min) return min;
    return v;
}


void CreateDCTDictionary(double* atoms, int blockW, int nbDCT)
{
    const int blockSize = blockW*blockW;
    const int nbAtoms   = nbDCT*nbDCT;
    const double pi     = 3.14159265;

    double* dct1D = new double[nbDCT*blockW];
    int offset = 0;
    for (int k = 0; k < nbDCT; k++)
    {
        double mean = 0.0;
        for (int i = 0; i < blockW; i++)
        {
            double value = cos(i*k*pi/nbDCT);
            dct1D[offset + i] = value;
            mean += value;
        }
        mean /= blockW;

        if (k > 0)
        {
            for (int i = 0; i < blockW; i++)
                dct1D[offset + i] -= mean;
        }
        offset += blockW;
    }

    double* atom = atoms;
    for (int j = 0; j < nbDCT; j++)
    for (int i = 0; i < nbDCT; i++)
    {
        double m = 0.0;

        int k = 0;
        for (int h = 0; h < blockW; h++)
        for (int w = 0; w < blockW; w++)
        {
            double a = dct1D[i*blockW + h];
            double b = dct1D[j*blockW + w];
            double v = a*b;

            atom[k++] = v;
            m += v*v;
        }

        // Normalise
        m = 1.0/sqrt(m);
        for (int k = 0; k < blockSize; k++)
            atom[k] *= m;

        atom += blockSize;
    }

    delete[] dct1D;
}


void ImageToSignals(const byte* image, double* signals, int w, int h, int blockW)
{
    int offset = 0;
    for (int y = 0; y < h; y += blockW)
    for (int x = 0; x < w; x += blockW)
    {
        for (int sy = 0; sy < blockW; sy++)
        for (int sx = 0; sx < blockW; sx++)
            signals[offset++] = image[(x + sx) + (y + sy)*w]/255.0;
    }
}

void SignalsToImage(double* signals, byte* image, int w, int h, int blockW)
{
    for (int y = 0; y < h; y += blockW)
    for (int x = 0; x < w; x += blockW)
    {
        for (int sy = 0; sy < blockW; sy++)
        for (int sx = 0; sx < blockW; sx++)
        {
            double signal = *signals++;
            signal = clamp(signal*255 + 0.5, 0, 255);
            image[(x + sx) + (y + sy)*w] = (byte)signal;
        }
    }
}

void UpdateSignals(
    const double* atoms, int blockSize, double* signals, int nbSignals,
    int maxEntries, int* nbEntries, int* indices, double* values)
{
    for (int j = 0; j < nbSignals; j++)
    {
        for (int i = 0; i < blockSize; i++)
            signals[i] = 0.0;

        int num = nbEntries[j];
        for (int k = 0; k < num; k++)
        {
            for (int i = 0; i < blockSize; i++)
                signals[i] += atoms[indices[k]*blockSize + i]*values[k];
        }

        signals += blockSize;
        indices += maxEntries;
        values  += maxEntries;
    }
}


int main(int argc, char* argv[])
{
    int w, h, comp;

    if (argc != 6)
    {
        printf("Syntax: %s <input image> <output image> atoms epsilon steps\n", argv[0]);
        return -1;
    }

    const int nbDCT = 11;
    const int nbAtoms = nbDCT*nbDCT;

    const int maxEntries = atoi(argv[3]);
    if (!(maxEntries >= 1 && maxEntries <= nbAtoms))
    {
        printf("Atom count must be in the range [1, %d]\n", nbAtoms);
        return -1;
    }

    const double epsilon = atof(argv[4]);
    if (epsilon < 0)
    {
        printf("Epsilon must be >= 0\n");
        return -1;
    }

    const int nbIters = atoi(argv[5]);
    if (nbIters < 0)
    {
        printf("Number of K-SVD steps must be >= 0\n");
        return -1;
    }

    // Load the image
    byte* image = stbi_load(argv[1], &w, &h, &comp, 0);

    if (!image)
    {
        printf("Failed to load input image %s \n", argv[1]);
        return -1;
    }

    if (w & 7 || h & 7)
    {
        printf("Wrong image dimensions (not a multiple of 8)\n");
        stbi_image_free(image);
        return -1;
    }

    if (comp != 1)
    {
        printf("Wrong image format (needs to be greyscale)\n");
        stbi_image_free(image);
        return -1;
    }

    const int blockW    = 8;
    const int blockSize = blockW*blockW;
    int nbSignals = w*h/blockSize;

    // Convert the image into 8x8 blocks ('signals')
    double* signals = new double[w*h];
    ImageToSignals(image, signals, w, h, blockW);

    // Create an initial over-complete DCT dictionary
    double* atoms = new double[nbAtoms*blockSize];
    CreateDCTDictionary(atoms, blockW, nbDCT);

    // Evolve the DCT dictionary via K-SVD
    KSVD(atoms, blockSize, nbAtoms, signals, nbSignals, maxEntries, epsilon, nbIters);

    int*    nbEntries = new int[nbSignals];
    int*    indices   = new int[nbSignals*maxEntries];
    double* values    = new double[nbSignals*maxEntries];

    // Approximate the signals using the resulting dictionary
    BOMP(atoms, blockSize, nbAtoms, signals, nbSignals,
        maxEntries, epsilon,
        nbEntries, indices, values);

    // Update the signals with compressed version
    UpdateSignals(atoms, blockSize, signals, nbSignals, maxEntries,
        nbEntries, indices, values);

    // Write out the reconstructed image
    SignalsToImage(signals, image, w, h, blockW);
    stbi_write_png(argv[2], w, h, comp, image, 0);

    delete[] signals;
    delete[] atoms;
    delete[] nbEntries;
    delete[] indices;
    delete[] values;

    stbi_image_free(image);

    return 0;
}