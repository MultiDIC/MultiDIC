// This function forms regions given a mask (logical array) input. 
// It calls the form_regions function from ncorr_lib

#include <mex.h>
#include <vector>
#include "standard_datatypes.h" 
#include "ncorr_datatypes.h" 
#include "ncorr_lib.h"

// ----------------------------------------------------//
// Main Class -----------------------------------------//
// ----------------------------------------------------//

class class_formregions {
public:
    // Constructor
    class_formregions(mxArray *plhs [ ],const mxArray *prhs [ ]);
    
    // Methods
    void analysis();

private:
    // Properties
    // Inputs:
    class_logical_array mask;   // standard datatype
    int cutoff;                 // standard datatype
    bool preservelength;        // standard datatype

    // Outputs:
    mxArray *mat_regions;
    bool *removed;

    // Other variables:
    int numfields;
};

class_formregions::class_formregions(mxArray *plhs [ ],const mxArray *prhs[ ]) {
    // Get inputs -------------------------------------//
    // input 1: mask
    get_logical_array(mask,prhs[0]);     
    // input 2: cutoff
    get_integer_scalar(cutoff,prhs[1]);  
    // input 3: preservelength - this parameters tells whether to preserve the width of the mask when forming the noderange/nodelist
    get_logical_scalar(preservelength,prhs[2]); 
        
    // Form/set Outputs -------------------------------//
    // output 1: regions
    mwSize dims[2] = {0,0};
    numfields = 7;
    const char *fieldnames[] = {"nodelist","noderange","leftbound","rightbound","upperbound","lowerbound","totalpoints"};
    mat_regions = mxCreateStructArray(2,dims,numfields,fieldnames);
    plhs[0] = mat_regions;
    
    // output 2: removed - tells user if small regions were removed due to the cutoff      
    mxArray *mat_removed = mxCreateLogicalMatrix(1,1);
    plhs[1] = mat_removed;    
    
    // Get outputs ------------------------------------//
    // output 1: regions
    // This output is a structure whose elements' size is not known beforehand, so get them in the analysis function
    
    // output 2: removed
    removed = mxGetLogicals(mat_removed);
}

// ----------------------------------------------------//
// Main Class Methods ---------------------------------//
// ----------------------------------------------------//

void class_formregions::analysis() {            
    // Initialize removed to false
    *removed = false; // Set to true if small regions were removed
    
    // Initialize regions to send to form_regions - used vec_struct_region
    std::vector<vec_struct_region> regions;
    
    // Call ncorr_lib form_regions function
    form_regions(regions,*removed,mask,cutoff,preservelength);
    
    // Resize structure and store info
    mwSize dims[2] = {1,(mwSize)regions.size()};
    mxSetDimensions(mat_regions,dims,2);
    mxSetData(mat_regions,mxCalloc(dims[0]*dims[1]*numfields,sizeof(mxArray *))); // Allocate memory for resized structure
    for (int i=0; i<(int)regions.size(); i++) {
        // Form fields
        mxArray *mat_nodelist = mxCreateDoubleMatrix(regions[i].height_nodelist, regions[i].width_nodelist, mxREAL);
        mxArray *mat_noderange = mxCreateDoubleMatrix(regions[i].height_nodelist, 1, mxREAL);
        mxArray *mat_leftbound = mxCreateDoubleMatrix(1, 1, mxREAL);
        mxArray *mat_rightbound = mxCreateDoubleMatrix(1, 1, mxREAL);
        mxArray *mat_upperbound = mxCreateDoubleMatrix(1, 1, mxREAL);
        mxArray *mat_lowerbound = mxCreateDoubleMatrix(1, 1, mxREAL);
        mxArray *mat_totalpoints = mxCreateDoubleMatrix(1, 1, mxREAL);

        // Get Data
        class_double_array nodelist;
        class_double_array noderange; 
        get_double_array(nodelist,mat_nodelist);
        get_double_array(noderange,mat_noderange);
        double *leftbound = mxGetPr(mat_leftbound);
        double *rightbound = mxGetPr(mat_rightbound);
        double *upperbound = mxGetPr(mat_upperbound);
        double *lowerbound = mxGetPr(mat_lowerbound);
        double *totalpoints = mxGetPr(mat_totalpoints);

        // Copy Data - make sure to cast to double since Matlab uses double 
        for (int j=0; j<regions[i].height_nodelist; j++) {   
            noderange.value[j] = (double)regions[i].noderange[j];
            for (int k=0; k<regions[i].noderange[j]; k++) {   
                nodelist.value[j+k*nodelist.height] = (double)regions[i].nodelist[j+k*regions[i].height_nodelist];
            }
        }
        *leftbound = (double)regions[i].leftbound;
        *rightbound = (double)regions[i].rightbound;
        *upperbound = (double)regions[i].upperbound;
        *lowerbound = (double)regions[i].lowerbound;
        *totalpoints = (double)regions[i].totalpoints;

        // Add fields to structure
        mxSetFieldByNumber(mat_regions,i,0,mat_nodelist);
        mxSetFieldByNumber(mat_regions,i,1,mat_noderange);
        mxSetFieldByNumber(mat_regions,i,2,mat_leftbound);
        mxSetFieldByNumber(mat_regions,i,3,mat_rightbound);
        mxSetFieldByNumber(mat_regions,i,4,mat_upperbound);
        mxSetFieldByNumber(mat_regions,i,5,mat_lowerbound);
        mxSetFieldByNumber(mat_regions,i,6,mat_totalpoints);    
    }
}

void mexFunction(int nlhs,mxArray *plhs[ ],int nrhs,const mxArray *prhs[ ]) {
    if (nrhs == 3 && nlhs == 2) {
        // Create formregions
        class_formregions formregions(plhs,prhs);
        
        // Run analysis and assign outputs
        formregions.analysis();
    } else {
        mexErrMsgTxt("Incorrect number of input and output arguments.\n");
    }
}
