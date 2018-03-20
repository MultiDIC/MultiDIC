// This function forms a boundary given a mask input

#include <mex.h>
#include <math.h>
#include <vector>
#include <algorithm>
#include "standard_datatypes.h"
#include "ncorr_datatypes.h"
#include "ncorr_lib.h"

// ----------------------------------------------------//
// Main Class -----------------------------------------//
// ----------------------------------------------------//

class class_formboundary {
public:
    // Constructor
    class_formboundary(mxArray *plhs[ ],const mxArray *prhs [ ]);
    
    // Methods
    void analysis();

private:
    // Properties
    // Inputs:
    class_integer_array point_init;       // standard datatype
    int direc;                            // standard datatype
    class_logical_array mask;             // standard datatype
    
    // Outputs: 
    mxArray *mat_boundary;
    class_double_array boundary;
};

class_formboundary::class_formboundary(mxArray *plhs[ ],const mxArray *prhs[ ]) {    
    // Get inputs -------------------------------------------------------//
    // input 1: point_init
    get_integer_array(point_init,prhs[0]);
    // input 2: direc
    get_integer_scalar(direc,prhs[1]);   
    // input 3: mask
    get_logical_array(mask,prhs[2]);   

    // Check inputs
    if (point_init.width == 2 && point_init.height == 1) {       
        // Form/set outputs ---------------------------------------------//
        // output 1: boundary
        // Size isn't known beforehand, so just set to 0 for now. 
        // This will be reset after boundary has been calculated
        mat_boundary = mxCreateDoubleMatrix(0, 0, mxREAL); 
        plhs[0] = mat_boundary;
        
        // Get outputs --------------------------------------------------//
        // output 1: boundary
        get_double_array(boundary,mat_boundary);
    } else {
        mexErrMsgTxt("Initiation point must be size 1x2.\n");
    }
}

// ----------------------------------------------------//
// Main Class Methods ---------------------------------//
// ----------------------------------------------------//

void class_formboundary::analysis() {    
    // Initialize vec_boundary and vec_point_init
    std::vector<std::vector<int> > vec_boundary;
    std::vector<int> vec_point_init(2,0);
    vec_point_init[0] = point_init.value[0];
    vec_point_init[1] = point_init.value[1];
            
    // Call ncorr_lib form_boundary function
    form_boundary(vec_boundary,vec_point_init,mask,direc);
        
    // Resize boundary then transfer data
    boundary.alloc((int)vec_boundary.size(),2);
        
    // Copy elements and convert to double
    for (int i=0; i<boundary.height; i++) {
        boundary.value[i] = (double)vec_boundary[i][0];
        boundary.value[i+boundary.height] = (double)vec_boundary[i][1];
    }
    
    // Set values
    mxSetPr(mat_boundary,boundary.value);
    mxSetM(mat_boundary,boundary.height);
    mxSetN(mat_boundary,boundary.width);
}

void mexFunction(int nlhs,mxArray *plhs[ ],int nrhs,const mxArray *prhs[ ]) {
    // Make sure number of input and output arguments are correct
    if (nrhs == 3 && nlhs == 1) {
        // Create formboundary
        class_formboundary formboundary(plhs,prhs);
        
        // Run analysis and assign outputs
        formboundary.analysis();
    } else {
        mexErrMsgTxt("Only three inputs and one output arguments allowed.\n");
    }
}
