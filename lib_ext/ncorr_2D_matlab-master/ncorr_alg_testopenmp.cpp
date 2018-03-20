// This function tests if OpenMP is actually working correctly.

#include <mex.h>
#include <vector>
#include "standard_datatypes.h"
#include "ncorr_datatypes.h"
#include "ncorr_lib.h"
            
#ifdef NCORR_OPENMP
    #include <omp.h>              // openmp header
#endif

// ----------------------------------------------------//
// Main Class -----------------------------------------//
// ----------------------------------------------------//

class class_testopenmp {
public:
    // Constructor
    class_testopenmp(mxArray *plhs [ ],const mxArray *prhs [ ]);
    
    // Methods
    void analysis();

private:
    // Properties
    // Inputs: None

    // Outputs:
    bool *enabled_openmp;
    
    // Other variables:
    int total_threads;
    std::vector<char> vec_enabled_thread; // Do NOT use vector<bool> - it is not safe to concurrently write to vector<bool>
};

class_testopenmp::class_testopenmp(mxArray *plhs[ ],const mxArray *prhs[ ]) {
    // Get inputs ---------------------------------------------------// 
    // None
    
    // Set total threads --------------------------------------------//
    total_threads = 4; // Any number greater than 1 should suffice
    
    // Thread enabled vector ----------------------------------------//
    vec_enabled_thread.resize(total_threads,0); // One for each thread, initialized to false
                
    // OpenMP Setup -------------------------------------------------//
    #ifdef NCORR_OPENMP
        // Set number of threads
        omp_set_num_threads(total_threads);
    #endif
                
    // Form/set outputs ---------------------------------------------//
    // output 1: enabled_openmp
    plhs[0] = mxCreateLogicalMatrix(1,1);
        
    // Get outputs --------------------------------------------------//
    // output 1: enabled_openmp
    enabled_openmp = mxGetLogicals(plhs[0]);
}

// ----------------------------------------------------//
// Main Class Methods ---------------------------------//
// ----------------------------------------------------//

void class_testopenmp::analysis() {     
    // Initialize enabled_openmp to true
    *enabled_openmp = true; 
            
    // Enter parallel region - anything inside here needs to be threadsafe
    #ifdef NCORR_OPENMP
    #pragma omp parallel 
    {       
    #endif
            
        #ifdef NCORR_OPENMP
            // Get thread number 
            int num_thread = omp_get_thread_num();
        #else
            // Set to zero if openmp is not enabled
            int num_thread = 0;
        #endif
                      
        // Each thread needs to set vec_enabled_thread to true
        vec_enabled_thread[num_thread] = true;

    #ifdef NCORR_OPENMP
    }    
    #endif
       
    // Check vec_enabled_thread to make sure all seeds processed correctly
    for (int i=0; i<total_threads; i++) {
        if (!vec_enabled_thread[i]) {
            *enabled_openmp = false; 
        }
    }
}

void mexFunction(int nlhs,mxArray *plhs[ ],int nrhs,const mxArray *prhs[ ]) {
    if (nrhs == 0 && nlhs == 1) {
        // Create testopenmp 
        class_testopenmp testopenmp(plhs,prhs);
        
        // Run analysis
        testopenmp.analysis();
    } else {
        // Thread safe because it is single threaded up to this point
        mexErrMsgTxt("Incorrect number of inputs or outputs.\n");
    }
}
