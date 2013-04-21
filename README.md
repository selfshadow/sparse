sparse
======

A sparse encoding testbed.

### Current Features
* BOMP: _Batch Orthogonal Matching Pursuit_ [1]
* KSVD: _Approximate K-SVD_ dictionary training algorithm [1]

### Examples
**ksvd_image**: a simple test app demonstrating how to use the BOMP and KSVD functions, with the syntax:

    ksvd_image <input image> <output image> atoms epsilon steps
    
### Acknowledgements
Sean Barrett for [stblib](http://code.google.com/p/stblib/).
    
### References

[1] Rubinstein, R., Zibulevsky, M. and Elad, M. "Efficient Implementation of the K-SVD Algorithm using Batch Orthogonal Matching Pursuit". Technical Report - CS Technion, April 2008.
