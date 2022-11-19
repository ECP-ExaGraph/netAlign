#include "netAlignKernel.h"

int main(int argc, char** argv)
{

    omp_set_num_threads(1);
    netAlign(argc, argv);
    return 0;
    
}

