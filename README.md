# ZHadronAnalysis
Code for studying Z-tagged hadron yields.

Currently works for 2018 Pb+Pb, 2017 pp, and associated MC.
TODO implement 2016 p+Pb.

To run the analysis, first run the code in AnalysisCode. This will process output from the CERN grid (which has minimal event selection and analysis performed) and create intermediate trees.
Then you can run the Macros which will generate histograms and enable plots to be made.
