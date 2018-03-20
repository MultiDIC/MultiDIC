// This function forms unioned region given a ROI and mask (logical array) input. 
// It calls the form_union function from ncorr_lib

#include <mex.h>
#include <vector>
#include <algorithm> 
#include "standard_datatypes.h" 
#include "ncorr_datatypes.h" 
#include "ncorr_lib.h"

// ----------------------------------------------------//
// Main Class -----------------------------------------//
// ----------------------------------------------------//

class class_formunion {
public:
    // Constructor
    class_formunion(mxArray *plhs [ ],const mxArray *prhs [ ]);
    
    // Methods
    void analysis();

private:
    // Properties
    // Inputs:
    std::vector<ncorr_class_roi> roi;   // ncorr datatype
    class_logical_array mask_union;     // standard datatype

    // Outputs:
    mxArray *mat_region_union;

    // Other variables:
    int numfields;
};

class_formunion::class_formunion(mxArray *plhs [ ],const mxArray *prhs[ ]) {
    // Get inputs -------------------------------------//
    // input 1: ROI
    get_rois(roi,prhs[0]); 
    // input 2: mask_union
    get_logical_array(mask_union,prhs[1]); 
    
    // Form/set Outputs -------------------------------//
    // output 1: region
    mwSize dims[2] = {0,0};
    numfields = 7;
    const char *fieldnames[] = {"nodelist","noderange","leftbound","rightbound","upperbound","lowerbound","totalpoints"};
    mat_region_union = mxCreateStructArray(2,dims,numfields,fieldnames);
    plhs[0] = mat_region_union;    
    
    // Get Outputs ------------------------------------//
    // output 1: region
    // This output is a structure whose elements' size is not known beforehand, so get them in the analysis function
}

// ----------------------------------------------------//
// Main Class Methods ---------------------------------//
// ----------------------------------------------------//

void class_formunion::analysis() {
    // Initialize region_union to send to form_union - use vec_struct_region
    std::vector<vec_struct_region> region_union;
    region_union.resize(roi[0].region.size()); // Note that size is preserved, even if the unioned region is empty
    
    // Call ncorr_lib form_union function - Do not do this in-place since it's possible for unioned region to be larger
    form_union(region_union,roi[0].region,mask_union,false);
        
    // Resize structure and store info
    mwSize dims[2] = {1,(mwSize)roi[0].region.size()};
    mxSetDimensions(mat_region_union,dims,2); 
    mxSetData(mat_region_union,mxCalloc(dims[0]*dims[1]*numfields,sizeof(mxArray *))); // Allocate memory for resized structure
    for (int i=0; i<(int)roi[0].region.size(); i++) {    
        // Form fields
        mxArray *mat_nodelist = mxCreateDoubleMatrix(region_union[i].height_nodelist,region_union[i].width_nodelist,mxREAL);
        mxArray *mat_noderange = mxCreateDoubleMatrix(region_union[i].height_nodelist,1,mxREAL);
        mxArray *mat_leftbound = mxCreateDoubleMatrix(1,1,mxREAL);
        mxArray *mat_rightbound = mxCreateDoubleMatrix(1,1,mxREAL);
        mxArray *mat_upperbound = mxCreateDoubleMatrix(1,1,mxREAL);
        mxArray *mat_lowerbound = mxCreateDoubleMatrix(1,1,mxREAL);
        mxArray *mat_totalpoints = mxCreateDoubleMatrix(1,1,mxREAL);

        // Get data
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
        for (int j=0; j<region_union[i].height_nodelist; j++) {
            noderange.value[j] = (double)region_union[i].noderange[j];
            for (int k=0; k<region_union[i].noderange[j]; k++) {
                nodelist.value[j+k*nodelist.height] = (double)region_union[i].nodelist[j+k*region_union[i].height_nodelist];
            }
        }
        *leftbound = (double)region_union[i].leftbound;
        *rightbound = (double)region_union[i].rightbound;
        *upperbound = (double)region_union[i].upperbound;
        *lowerbound = (double)region_union[i].lowerbound;
        *totalpoints = (double)region_union[i].totalpoints;

        // Add fields to structure
        mxSetFieldByNumber(mat_region_union,i,0,mat_nodelist);
        mxSetFieldByNumber(mat_region_union,i,1,mat_noderange);
        mxSetFieldByNumber(mat_region_union,i,2,mat_leftbound);
        mxSetFieldByNumber(mat_region_union,i,3,mat_rightbound);
        mxSetFieldByNumber(mat_region_union,i,4,mat_upperbound);
        mxSetFieldByNumber(mat_region_union,i,5,mat_lowerbound);
        mxSetFieldByNumber(mat_region_union,i,6,mat_totalpoints);    
    }
}

void mexFunction(int nlhs,mxArray *plhs[ ],int nrhs,const mxArray *prhs[ ]) {
    if (nrhs == 2 && nlhs == 1) {
        // Create formunion
        class_formunion formunion(plhs,prhs);
        
        // Run analysis and assign outputs
        formunion.analysis();
    } else {
        mexErrMsgTxt("Only two input and one output arguments.\n");
    }
}
