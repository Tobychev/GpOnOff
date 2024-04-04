Gammapy analysis scrips 
======================

Some scripts for analysing HESS data with gammapy, primarily used for loose ImPACT data, which have the quirk that they have to be manually exported and lack a background IRF. For this reason some extra work is required.

The two scrips are StandardReflectedAnalysis and StandardRingAnalysis (will eventually) reduce runs into gammapy datasets, while SpectralAnalysis will do a forward folding spectral fit on the reflected dataset.
