// This function extrapolates data by expanding and filtering it

#include <mex.h>
#include <vector>
#include "standard_datatypes.h"
#include "ncorr_datatypes.h"
#include "ncorr_lib.h"

// ----------------------------------------------------//
// Main Class -----------------------------------------//
// ----------------------------------------------------//

class class_extrapdata {
public:
    // Constructor
    class_extrapdata(mxArray *plhs[ ],const mxArray *prhs [ ]);
    
    // Methods:
    void analysis();
    
private:
    // Properties
    // Inputs:
    class_double_array plot_data;                    // standard datatype
    std::vector<ncorr_class_roi> roi;                // ncorr datatype
    int border_extrap;                               // standard datatype
    
    // Outputs:
    std::vector<class_double_array> plot_extrap;
};

class_extrapdata::class_extrapdata(mxArray *plhs[ ],const mxArray *prhs[ ]) {
    // Get inputs ---------------------------------------------------//
    // input 1: plot_data
    get_double_array(plot_data,prhs[0]);
    // input 2: ROI
    get_rois(roi,prhs[1]);
    // input 3: border_extrap
    get_integer_scalar(border_extrap,prhs[2]);
        
    // Form/set outputs ---------------------------------------------//
    // output 1: plot_extrap
    int numdims = 2;
    mwSize dims[2] = {1,(mwSize)roi[0].region.size()};
    mxArray *mat_plot_extrap = mxCreateCellArray(numdims, dims);
    plhs[0] = mat_plot_extrap;    
    // Form one plot_data per region
    for (int i=0; i<(int)roi[0].region.size(); i++) {
        // Get dimensions of plot_extrap
        int height_plot_extrap = 0;
        int width_plot_extrap = 0;        
        if (roi[0].region[i].totalpoints > 0) {
            height_plot_extrap = (roi[0].region[i].lowerbound-roi[0].region[i].upperbound+1)+2*border_extrap;
            width_plot_extrap = (roi[0].region[i].rightbound-roi[0].region[i].leftbound+1)+2*border_extrap;
        }           
        
        // Create array
        mxArray *mat_plot_extrap_buf = mxCreateDoubleMatrix(height_plot_extrap, width_plot_extrap, mxREAL);
        // Add array to cell
        mxSetCell(mat_plot_extrap,i,mat_plot_extrap_buf);   

        // Get outputs -------------------------------------------------//
        // output 1: plot_extrap
        class_double_array plot_extrap_template;
        get_double_array(plot_extrap_template,mat_plot_extrap_buf);
        plot_extrap.push_back(plot_extrap_template);
    }
}

// ----------------------------------------------------//
// Main Class Methods ---------------------------------//
// ----------------------------------------------------//

void class_extrapdata::analysis() {    
    // Cycle over regions 
    for (int i=0; i<(int)roi[0].region.size(); i++) {   
        // Fill plot_extrap with values from plot_data. These points will not be altered. 
        // Only the region outside this region will be altered/extrapolated.        
        if (roi[0].region[i].totalpoints > 0) {
            for (int j=0; j<roi[0].region[i].noderange.height; j++) {
                int x = j+border_extrap;
                for (int k=0; k<roi[0].region[i].noderange.value[j]; k+=2) {
                    for (int l=roi[0].region[i].nodelist.value[j+k*roi[0].region[i].nodelist.height]; l<=roi[0].region[i].nodelist.value[j+(k+1)*roi[0].region[i].nodelist.height]; l++) {
                        int y = l-roi[0].region[i].upperbound+border_extrap;
                        plot_extrap[i].value[y+x*plot_extrap[i].height] = plot_data.value[l+(j+roi[0].region[i].leftbound)*plot_data.height];
                    }
                }
            }

            // Form inverseregion
            ncorr_class_inverseregion inverseregion(roi[0].region[i],border_extrap);

            // Expand and filter data
            expand_filt(plot_extrap[i],inverseregion);
        }
    }
}

void mexFunction(int nlhs,mxArray *plhs[ ],int nrhs,const mxArray *prhs[ ]) {
    if (nrhs == 3 && nlhs == 1) {
        // From extrapdata 
        class_extrapdata extrapdata(plhs,prhs);
        
        // Run analysis and assign outputs
        extrapdata.analysis();
    } else {
        mexErrMsgTxt("Incorrect number of inputs or outputs.\n");
    }
}
