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


void CreateDCTDictionary(double* atoms, int blockW, int nbDCT, int nbChannels)
{
    const int blockSize = blockW*blockW*nbChannels;
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

            for (int c = 0; c < nbChannels; c++)
            {
                atom[k++] = v;
                m += v*v;
            }
        }

        // Normalise
        m = 1.0/sqrt(m);
        for (int k = 0; k < blockSize; k++)
            atom[k] *= m;

        atom += blockSize;
    }

    delete[] dct1D;
}


void ImageToSignals(const byte* image, double* signals, int w, int h, int blockW, int blockH)
{
    int offset = 0;
    for (int y = 0; y < h; y += blockH)
    for (int x = 0; x < w; x += blockW)
    {
        for (int sy = 0; sy < blockH; sy++)
        for (int sx = 0; sx < blockW; sx++)
            signals[offset++] = image[(x + sx) + (y + sy)*w]/255.0;
    }
}

void SignalsToImage(double* signals, byte* image, int w, int h, int blockW, int blockH)
{
    for (int y = 0; y < h; y += blockH)
    for (int x = 0; x < w; x += blockW)
    {
        for (int sy = 0; sy < blockH; sy++)
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
    if (argc != 6)
    {
        printf("Syntax: %s <input image> <output image> atoms error steps\n", argv[0]);
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

    const double error = atof(argv[4]);
    if (error < 0)
    {
        printf("Error must be >= 0\n");
        return -1;
    }

    const int nbIters = atoi(argv[5]);
    if (nbIters < 0)
    {
        printf("Number of K-SVD steps must be >= 0\n");
        return -1;
    }

    // Load the image
    int w, h, nbChannels;
    byte* image = stbi_load(argv[1], &w, &h, &nbChannels, 0);

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

    const int stride    = w*nbChannels;
    const int blockH    = 8;
    const int blockW    = blockH*nbChannels;
    const int blockSize = blockW*blockH;
    const int nbSignals = stride*h/blockSize;

    // Convert the image into 8x8 blocks ('signals')
    double* signals = new double[stride*h];
    ImageToSignals(image, signals, stride, h, blockW, blockH);

    // Create an initial over-complete DCT dictionary
    double* atoms = new double[nbAtoms*blockSize];
    CreateDCTDictionary(atoms, blockH, nbDCT, nbChannels);

    // Evolve the DCT dictionary via K-SVD
    KSVD(atoms, blockSize, nbAtoms, signals, nbSignals, maxEntries, error, nbIters);

    int*    nbEntries = new int[nbSignals];
    int*    indices   = new int[nbSignals*maxEntries];
    double* values    = new double[nbSignals*maxEntries];

    // Approximate the signals using the resulting dictionary
    BOMP(atoms, blockSize, nbAtoms, signals, nbSignals,
        maxEntries, error,
        nbEntries, indices, values);

    // Update the signals with the compressed version
    UpdateSignals(atoms, blockSize, signals, nbSignals, maxEntries,
        nbEntries, indices, values);

    // Write out the reconstructed image
    SignalsToImage(signals, image, stride, h, blockW, blockH);
    stbi_write_png(argv[2], w, h, nbChannels, image, 0);

    delete[] signals;
    delete[] atoms;
    delete[] nbEntries;
    delete[] indices;
    delete[] values;

    stbi_image_free(image);

    return 0;
}