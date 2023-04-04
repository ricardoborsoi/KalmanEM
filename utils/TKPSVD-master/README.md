Tensor-based Kronecker Product Singular Value Decomposition (TKPSVD) for Matlab&copy;/Octave&copy;
--------------------------------------------------------------------------------------------------

The Tensor-based Kronecker Product Singular Value Decomposition (TKPSVD) decomposes an arbitrary tensor A into a unique linear combination of d Kronecker products, characterized by an easy quantification of the relative approximation error. Furthermore, if A is a general symmetric tensor (symmetric, persymmetric, centrosymmetric,...) or has a shifted-index structure (Toeplitz, Hankel, constant fibres), then the Kronecker factors are guaranteed to inherit this structure.

1. Functions
------------

* [B,sigmas] = tkpsvd(A,n) or [B,sigmas] = tkpsvd(A,n,R)

Use this function to compute the TKPSVD for tensor A with Kronecker product factors dimensions as specified in n.

* demo.m

Small demo that illustrates the use of most functions in this pacakge.


2. Reference
------------

"A constructive arbitrary-degree Kronecker product decomposition of tensors"

Authors: Kim Batselier, Ngai Wong
