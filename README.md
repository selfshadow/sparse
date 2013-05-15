sparse
======

A sparse encoding testbed.

### Current Features
* BOMP: _Batch Orthogonal Matching Pursuit_ [1]
* KSVD: _Approximate K-SVD_ dictionary training algorithm [1]

### Examples
**ksvd_image**: a simple test app demonstrating how to use the BOMP and KSVD functions, with the syntax:

    ksvd_image <input image> <output image> atoms error steps

For instance:
    
    ksvd_image data/barbara.png babs_out.png 4 0.01 10

### Building
Use [premake](http://industriousone.com/premake) to generate the appropriate files for your system:

    premake4 gmake    # for GNU makefiles using GCC
    premake4 vs2010   # for a Visual Studio 2010 solution

### Acknowledgements
Robin Green and Manny Ko for their GDC presentation [2] and helpful correspondence.  
Ron Rubinstein for [OMP-Box](http://www.cs.technion.ac.il/~ronrubin/software.html) and [KSVD-Box](http://www.cs.technion.ac.il/~ronrubin/software.html).  
Sean Barrett for [stblib](http://code.google.com/p/stblib/).
    
### References
[1] Rubinstein, R., Zibulevsky, M. and Elad, M. ["Efficient Implementation of the K-SVD Algorithm using Batch Orthogonal Matching Pursuit"](http://www.cs.technion.ac.il/~ronrubin/research.html). Technical Report - CS Technion, April 2008.  
[2] Green, R., Ko, M. ["Orthogonal Matching Pursuit and K-SVD for Sparse Encoding"](https://docs.google.com/file/d/0B0AdsFXmX1ORRVFvajgyNWN4OE0/edit). Game Developer's Conference, March 2013.
