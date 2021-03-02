# FIF2
Multidimensional Fast Iterative Filtering for the decompostion of 2D non-stationary signals [1,2].

Please refer to "Example_v3.m" for an example of how to use the code.

It is based on FFT, which makes FIF2 to be really fast [2,3]. This implies that it is required a periodical extension at the boundaries.

To overcome this limitation we can preextend the signal under investigation [4].

Please cite our works:

[1] A. Cicone, H. Zhou. "Multidimensional Iterative Filtering method for the decomposition of high-dimensional non-stationary signals".
    Cambridge Core in Numerical Mathematics: Theory, Methods and Applications, Volume 10, Issue 2, Pages 278-298, 2017.
    doi:10.4208/nmtma.2017.s05

[2] S. Sfarra, A. Cicone, B. Yousefi, S. Perilli, L. Robol, X. P.V. Maldague. 
    "Maximizing the detection of thermal imprints in civil engineering composites after a thermal stimulus - The contribution of an
    innovative mathematical pre-processing tool: the 2D Fast Iterative Filtering algorithm. Philosophy, comparisons, numerical, qualitative
    and quantitative results". 2021. Submitted

[3] A. Cicone, H. Zhou. "Numerical Analysis for Iterative Filtering with New Efficient Implementations Based on FFT". 
    Numerische Mathematik, 147 (1), pages 1-28, 2021. doi: 10.1007/s00211-020-01165-5

[4] A. Stallone, A. Cicone, M. Materassi. "New insights and best practices for the successful use of Empirical Mode Decomposition, 
    Iterative Filtering and derived algorithms". Scientific Reports, Volume 10, article number 15161, 2020. doi: 10.1038/s41598-020-72193-2
