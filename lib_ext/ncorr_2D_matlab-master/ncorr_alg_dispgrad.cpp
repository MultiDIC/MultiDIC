// This function calculates displacement gradients given displacement field inputs

#include <mex.h>
#include <math.h>
#include <vector>
#include "standard_datatypes.h"
#include "ncorr_datatypes.h"
#include "ncorr_lib.h"

// ----------------------------------------------------//
// Main Class -----------------------------------------//
// ----------------------------------------------------//

class class_dispgrad {
public:
    // Constructor
    class_dispgrad(mxArray *plhs [ ],const mxArray *prhs [ ]);
    
    // Methods
    void analysis();

private:
    // Properties
    // Inputs:
    class_double_array plot_u;                // standard datatype
    class_double_array plot_v;                // standard datatype
    std::vector<ncorr_class_roi> roi;         // ncorr datatype
    int radius_strain;                        // standard datatype
    double pixtounits;                        // standard datatype
    int spacing;                              // standard datatype
    bool subsettrunc;                         // standard datatype
    int num_img;                              // standard datatype
    int total_imgs;                           // standard datatype
        
    // Outputs:
    class_double_array plot_dudx;
    class_double_array plot_dudy;
    class_double_array plot_dvdx;
    class_double_array plot_dvdy;
    class_logical_array plot_validpoints;
    double *outstate;
    
    // Other variables:
    class_waitbar waitbar;
};

class_dispgrad::class_dispgrad(mxArray *plhs[ ],const mxArray *prhs[ ]){
    // Get inputs ------------------------------------------------------//
    // input 1: u plot
    get_double_array(plot_u,prhs[0]);
    // input 2: v plot
    get_double_array(plot_v,prhs[1]);
    // input 3: ROI
    get_rois(roi,prhs[2]);
    // input 4: strain radius
    get_integer_scalar(radius_strain,prhs[3]);
    // input 6: pixtounits - this is units/pixel used as a conversion
    get_double_scalar(pixtounits,prhs[4]);
    // input 5: spacing
    get_integer_scalar(spacing,prhs[5]);
    // input 7: subsettrunc
    get_logical_scalar(subsettrunc,prhs[6]);
    // input 8: image id
    get_integer_scalar(num_img,prhs[7]);
    // input 9: total images
    get_integer_scalar(total_imgs,prhs[8]);
    
    // Check inputs - rudimentary check
    if (plot_u.width == roi[0].mask.width && plot_u.height == roi[0].mask.height &&
        plot_v.width == roi[0].mask.width && plot_v.height == roi[0].mask.height) {
        // Set cirroi ---------------------------------------------------//
        // Only allocate one cirroi
        roi[0].set_cirroi(radius_strain,1);
        
        // Form/set outputs ---------------------------------------------//
        // output 1: plot_dispgrad
        // Form displacement gradient structure
        mwSize dims[2] = {1,1};
        int numfields = 5;
        const char *fieldnames[] = {"plot_dudx","plot_dudy","plot_dvdx","plot_dvdy","plot_validpoints"};
        plhs[0] = mxCreateStructArray(2,dims,numfields,fieldnames);

        // Form fields
        mxArray *mat_plot_dudx = mxCreateDoubleMatrix(roi[0].mask.height,roi[0].mask.width,mxREAL);
        mxArray *mat_plot_dudy = mxCreateDoubleMatrix(roi[0].mask.height,roi[0].mask.width,mxREAL);
        mxArray *mat_plot_dvdx = mxCreateDoubleMatrix(roi[0].mask.height,roi[0].mask.width,mxREAL);
        mxArray *mat_plot_dvdy = mxCreateDoubleMatrix(roi[0].mask.height,roi[0].mask.width,mxREAL);
        mxArray *mat_plot_validpoints = mxCreateLogicalMatrix(roi[0].mask.height,roi[0].mask.width); 

        // Add fields to structure
        // add dudx:
        mxSetFieldByNumber(plhs[0],0,0,mat_plot_dudx);
        // add dudy:
        mxSetFieldByNumber(plhs[0],0,1,mat_plot_dudy);
        // add dvdx:
        mxSetFieldByNumber(plhs[0],0,2,mat_plot_dvdx);
        // add dvdy:
        mxSetFieldByNumber(plhs[0],0,3,mat_plot_dvdy);
        // add validpoints:
        mxSetFieldByNumber(plhs[0],0,4,mat_plot_validpoints);        

        // output 2: outstate
        plhs[1] = mxCreateDoubleMatrix(1,1,mxREAL);
        
        // Get outputs --------------------------------------------------//
        // output 1: plot_dispgrad
        // dudx:
        get_double_array(plot_dudx,mat_plot_dudx);
        // dudy:
        get_double_array(plot_dudy,mat_plot_dudy);
        // dvdx:
        get_double_array(plot_dvdx,mat_plot_dvdx);
        // dvdy:
        get_double_array(plot_dvdy,mat_plot_dvdy);
        // validpoints:
        get_logical_array(plot_validpoints,mat_plot_validpoints);
        // output 2: outstate
        outstate = mxGetPr(plhs[1]);
    } else {
        mexErrMsgTxt("Mask and displacements are not the same size.\n");
    }
}

// ----------------------------------------------------//
// Main Class Methods ---------------------------------//
// ----------------------------------------------------//

void class_dispgrad::analysis() {
    // Initialize outstate to cancelled
    *outstate = (double)CANCELLED; 

    // Set up waitbar ----------------------------------------------//
    int computepoints = 0;
    for (int i=0; i<(int)roi[0].region.size(); i++) {
        computepoints += roi[0].region[i].totalpoints;
    }
    waitbar.start(num_img,total_imgs,computepoints);    
    
    // Begin Analysis ----------------------------------------------//
    std::vector<double> mat_LS(9,0);
    std::vector<double> u_vec_LS(3,0);
    std::vector<double> v_vec_LS(3,0);
    
    // Cycle over each region and calculate displacement gradients
    for (int i=0; i<(int)roi[0].region.size(); i++) {
        // Update cirroi
        roi[0].update_cirroi(i,0);
        for (int j=0; j<roi[0].region[i].noderange.height; j++) {
            int x = j+roi[0].region[i].leftbound;
            for (int k=0; k<roi[0].region[i].noderange.value[j]; k+=2) {
                for (int l=roi[0].region[i].nodelist.value[j+k*roi[0].region[i].nodelist.height]; l<=roi[0].region[i].nodelist.value[j+(k+1)*roi[0].region[i].nodelist.height]; l++) {
                    int y = l;
                                        
                    // Get cirroi corresponding to thread 0 since this is not multithreaded
                    roi[0].get_cirroi(x,y,i,subsettrunc,0);    
                    
                    // Iterate over subset nodes
                    // Initialize to zero first
                    std::fill(mat_LS.begin(),mat_LS.end(),0.0);
                    std::fill(u_vec_LS.begin(),u_vec_LS.end(),0.0);
                    std::fill(v_vec_LS.begin(),v_vec_LS.end(),0.0);
                    for (int m=0; m<roi[0].cirroi[0].region.noderange.height; m++) {
                        int idx_x_LS = m+x-roi[0].cirroi[0].radius;
                        for (int n=0; n<roi[0].cirroi[0].region.noderange.value[m]; n+=2) {
                            for (int p=roi[0].cirroi[0].region.nodelist.value[m+n*roi[0].cirroi[0].region.nodelist.height]; p<=roi[0].cirroi[0].region.nodelist.value[m+(n+1)*roi[0].cirroi[0].region.nodelist.height]; p++) {
                                int idx_y_LS = p;

                                double x_LS = ((double)m-(double)roi[0].cirroi[0].radius);
                                double y_LS = ((double)p-(double)y);
                                double u_LS = plot_u.value[idx_y_LS+idx_x_LS*plot_u.height];
                                double v_LS = plot_v.value[idx_y_LS+idx_x_LS*plot_v.height];

                                // Do matrix first
                                mat_LS[0] += pow(x_LS,2);
                                mat_LS[3] += x_LS*y_LS;
                                mat_LS[4] += pow(y_LS,2);
                                mat_LS[6] += x_LS;
                                mat_LS[7] += y_LS;

                                // Do vec next
                                u_vec_LS[0] += x_LS*u_LS;
                                u_vec_LS[1] += y_LS*u_LS;
                                u_vec_LS[2] += u_LS;

                                v_vec_LS[0] += x_LS*v_LS;
                                v_vec_LS[1] += y_LS*v_LS;
                                v_vec_LS[2] += v_LS;
                            }
                        }
                    }
                    
                    // Fill symmetric parts of matrix
                    mat_LS[1] = mat_LS[3];
                    mat_LS[2] = mat_LS[6];
                    mat_LS[5] = mat_LS[7];
                    mat_LS[8] = (double)roi[0].cirroi[0].region.totalpoints;

                    // Find new plane parameters
                    // Solve for new parameters via cholesky decomp (from Golub and Van Loan)
                    // Lower triangle of mat_LS overwritten with parameters used in cholesky decomp
                    // v_vec_LS and u_vec_LS are overwritten with displacement gradient info
                    bool positivedef = true;
                    cholesky(mat_LS,positivedef,3);
                    if (positivedef) {
                        // Ax = b
                        // GG'x = b, where G is lower triangular
                        // Gy = b -> G'x = y
                        // Step 1: solve for y with forward substitution; y is stored in gradient_buffer
                        forwardsub(u_vec_LS,mat_LS,3);
                        forwardsub(v_vec_LS,mat_LS,3);
                        
                        // Step 2: solve for x with back substitution
                        backwardsub(u_vec_LS,mat_LS,3);
                        backwardsub(v_vec_LS,mat_LS,3);
                        
                        // Normalize values to account for spacing and displacement unit conversion
                        for (int i=0; i<3; i++) {
                            u_vec_LS[i] /= (double)(spacing+1)*pixtounits;
                            v_vec_LS[i] /= (double)(spacing+1)*pixtounits;
                        }
                    
                        // Now store displacement gradients
                        plot_dudx.value[y+x*plot_dudx.height] = u_vec_LS[0];
                        plot_dudy.value[y+x*plot_dudy.height] = u_vec_LS[1];
                        plot_dvdx.value[y+x*plot_dvdx.height] = v_vec_LS[0];
                        plot_dvdy.value[y+x*plot_dvdy.height] = v_vec_LS[1];
                        plot_validpoints.value[y+x*plot_validpoints.height] = true;
                    }
                    
                    // Update, check, and increment waitbar
                    if (!waitbar.updateandcheck()) {
                        // Waitbar was cancelled
                        return;
                    }
                    waitbar.increment();
                }
            }
        }
    }
    
    // At this point analysis has been completed successfully
    *outstate = (double)SUCCESS;
}

void mexFunction(int nlhs,mxArray *plhs[ ],int nrhs,const mxArray *prhs[ ]) {
    if (nrhs == 9 && nlhs == 2) {
        // Create dispgrad 
        class_dispgrad dispgrad(plhs,prhs);
        
        // Run analysis
        dispgrad.analysis();
    } else {
        mexErrMsgTxt("Incorrect number of inputs or outputs.\n");
    }
}
