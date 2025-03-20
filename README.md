# ErrorEstimationForH1
- Target: HybridH1vsMixed
- The adaptive algorithm is selected within the EstimateError function (line 287 of the Methods/Solver.cpp file) To apply algorithm 3 you must set the index of the node that corresponds to the singularity
- The domain type: square centered or not centered at the origin or L-shaped domain, is selected within the Configure function (line 29 of the Methods/InputTreatment.cpp file)
- The estimated errors and effectiveness indices of the hp-adaptive process are published in the file Projects/HybridH1/.../adaptivityresults.txt
