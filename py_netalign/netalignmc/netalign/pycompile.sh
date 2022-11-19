c++ -O3 -shared -std=c++11 -fopenmp -fPIC $(python3 -m pybind11 --includes) matchingHalfApproxDominatingNew.c parseInputFiles.c sparseImpl.cpp netAlignExactMatch.cpp netAlignImplMP.cpp netAlignImplMR.cpp netAlignImpl.cpp  netAlignPY.cpp -o netAlignPY.so

mv netAlignPY.so ../../src
