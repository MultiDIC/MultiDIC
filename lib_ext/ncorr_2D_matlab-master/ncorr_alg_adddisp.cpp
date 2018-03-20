// This function adds displacement plots together

#include <mex.h>
#include <math.h>
#include <vector>
#include <list>
#include "standard_datatypes.h"
#include "ncorr_datatypes.h"
#include "ncorr_lib.h"

// ----------------------------------------------------//
// Interp plot input ----------------------------------//
// ----------------------------------------------------//

void get_plot_interp(std::vector<std::vector<class_double_array> > &plot_interp,const mxArray *prhs) {    
    // Check input - it's a cells containing cells. The outer cell corresponds to the image, while the inner cells correspond to regions
    if (mxIsClass(prhs,"cell") && mxGetM(prhs) == 1) {
        for (int i=0; i<(int)mxGetN(prhs); i++) {
            // Get content of cell corresponding to image number
            mxArray *mat_plot_interp_img = mxGetCell(prhs,i);
            std::vector<class_double_array> plot_interp_img_template;
            if (mat_plot_interp_img != 0 && mxIsClass(mat_plot_interp_img,"cell") && mxGetM(mat_plot_interp_img) == 1) {
                for (int j=0; j<(int)mxGetN(mat_plot_interp_img); j++) {
                    // Get content of cell corresponding to region
                    mxArray *mat_plot_interp_region = mxGetCell(mat_plot_interp_img,j);
                    if (mat_plot_interp_region != 0) {
                        // Form interp plot
                        class_double_array plot_interp_region_template;
                        get_double_array(plot_interp_region_template,mat_plot_interp_region);
                        
                        // Store interp plot per region
                        plot_interp_img_template.push_back(plot_interp_region_template);
                    } else {
                        mexErrMsgTxt("Some cell contents are empty.\n");
                    }
                }
                // Store interp plot per img
                plot_interp.push_back(plot_interp_img_template);
            } else {
                mexErrMsgTxt("Some cell contents are empty.\n");
            }
        }
    } else {
        mexErrMsgTxt("Interp data must be a row vector of class 'cell'.\n");
    }
}

// ----------------------------------------------------//
// Main Class -----------------------------------------//
// ----------------------------------------------------//

class class_adddisp {
public:
    // Constructor
    class_adddisp(mxArray *plhs [ ],const mxArray *prhs [ ]);
    
    // Methods:
    void analysis();

private:
    // Properties
    // Inputs:
    std::vector<std::vector<class_double_array> > plots_u_interp;            // local datatype
    std::vector<std::vector<class_double_array> > plots_v_interp;            // local datatype
    std::vector<ncorr_class_roi> rois_interp;                                // ncorr datatype        
    int border_interp;                                                       // standard datatype
    int spacing;                                                             // standard datatype
    int num_img;                                                             // standard datatype
    int total_imgs;                                                          // standard datatype
        
    // Outputs:
    class_double_array plot_u_added;
    class_double_array plot_v_added;
    class_logical_array plot_validpoints;
    double *outstate;
        
    // Other variables:
    class_waitbar waitbar;
};

class_adddisp::class_adddisp(mxArray *plhs[ ],const mxArray *prhs[ ]) {
    // Get inputs -------------------------------------------------//
    // input 1: plots_u_interp
    get_plot_interp(plots_u_interp,prhs[0]);
    // input 2: plots_v_interp
    get_plot_interp(plots_v_interp,prhs[1]);
    // input 3: rois_interp - used for determining bounds when forward propagating
    get_rois(rois_interp,prhs[2]);
    // input 4: border_interp
    get_integer_scalar(border_interp,prhs[3]);
    // input 5: spacing
    get_integer_scalar(spacing,prhs[4]);
    // input 6: num_img
    get_integer_scalar(num_img,prhs[5]);
    // input 7: total_imgs
    get_integer_scalar(total_imgs,prhs[6]);

    // Check inputs - check sizes - UPDATE THIS LATER!!!!!   
        
    // Form/set outputs -------------------------------------------//
    // output 1: plot_adddisp
    // Form deformation structure
    mwSize def_dims[2] = {1,1};
    int def_numfields = 3;
    const char *def_fieldnames[] = {"plot_u_added","plot_v_added","plot_validpoints"};
    plhs[0] = mxCreateStructArray(2, def_dims, def_numfields, def_fieldnames);

    // Form fields
    mxArray *mat_plot_u_added = mxCreateDoubleMatrix(rois_interp[0].mask.height, rois_interp[0].mask.width, mxREAL);
    mxArray *mat_plot_v_added = mxCreateDoubleMatrix(rois_interp[0].mask.height, rois_interp[0].mask.width, mxREAL);
    mxArray *mat_plot_validpoints = mxCreateLogicalMatrix(rois_interp[0].mask.height, rois_interp[0].mask.width);

    // Add fields to structure
    // Add u:
    mxSetFieldByNumber(plhs[0],0,0,mat_plot_u_added);
    // Add v:
    mxSetFieldByNumber(plhs[0],0,1,mat_plot_v_added);
    // Add valid points:
    mxSetFieldByNumber(plhs[0],0,2,mat_plot_validpoints);

    // output 2: outstate
    plhs[1] = mxCreateDoubleMatrix(1,1,mxREAL);
    
    // Get outputs -----------------------------------------------//
    // output 1: plot_adddisp
    // u:
    get_double_array(plot_u_added,mat_plot_u_added);
    // v:
    get_double_array(plot_v_added,mat_plot_v_added);
    // plot_validpoints
    get_logical_array(plot_validpoints,mat_plot_validpoints);
    // output 2: outstate
    outstate = mxGetPr(plhs[1]);
}

// ----------------------------------------------------//
// Main Class Methods ---------------------------------//
// ----------------------------------------------------//

void class_adddisp::analysis() {
    // Initialize outstate to cancelled
    *outstate = (double)CANCELLED; 
    
    // Set up waitbar ----------------------------------------------//
    int computepoints = 0; 
    for (int i=0; i<(int)rois_interp[0].region.size(); i++) {
        computepoints += rois_interp[0].region[i].totalpoints;
    }
    waitbar.start(num_img,total_imgs,computepoints);  
    
    // Cycle over each region in the first ROI
    for (int i=0; i<(int)rois_interp[0].region.size(); i++) {
        // Cycle over each point in each region in the first ROI
        for (int j=0; j<rois_interp[0].region[i].noderange.height; j++) {
            int x_ref_reduced = j+rois_interp[0].region[i].leftbound;
            for (int k=0; k<rois_interp[0].region[i].noderange.value[j]; k+=2) {
                for (int m=rois_interp[0].region[i].nodelist.value[j+k*rois_interp[0].region[i].nodelist.height]; m<=rois_interp[0].region[i].nodelist.value[j+(k+1)*rois_interp[0].region[i].nodelist.height]; m++) {
                    int y_ref_reduced = m;     
                    
                    // For each point (x,y), forward propagate this point 
                    // by using the bspline coefficients corresponding to i
                    bool propsuccess = true; // Set this to false if the point cant be propagated
                    
                    // Cycle over each plot and propagate until finished
                    double x_cur_reduced = (double)x_ref_reduced; 
                    double y_cur_reduced = (double)y_ref_reduced;
                    double u_interp = 0;
                    double v_interp = 0;
                    double u_interp_buf = 0; 
                    double v_interp_buf = 0;
                    for (int n=0; n<(int)rois_interp.size(); n++) {
                        // Get interpolated u and v displacement values
                        if (plots_u_interp[n][i].height > 0 && plots_u_interp[n][i].width > 0 && plots_v_interp[n][i].height > 0 && plots_v_interp[n][i].width > 0 &&
                            interp_qbs(u_interp_buf,x_cur_reduced,y_cur_reduced,plots_u_interp[n][i],rois_interp[n].mask,rois_interp[n].region[i].leftbound,rois_interp[n].region[i].upperbound,border_interp) == SUCCESS && 
                            interp_qbs(v_interp_buf,x_cur_reduced,y_cur_reduced,plots_v_interp[n][i],rois_interp[n].mask,rois_interp[n].region[i].leftbound,rois_interp[n].region[i].upperbound,border_interp) == SUCCESS) {
                            // Update displacements
                            u_interp += u_interp_buf;
                            v_interp += v_interp_buf;
                            
                            // Update x_cur and y_cur
                            x_cur_reduced += u_interp_buf/(spacing+1);
                            y_cur_reduced += v_interp_buf/(spacing+1);
                        } else {
                            propsuccess = false;
                            break;
                        }
                    }
                    
                    // If point was successfully propagated then store it
                    if (propsuccess) {
                        // Store final u_interp and v_interp value
                        plot_u_added.value[y_ref_reduced + x_ref_reduced*plot_u_added.height] = u_interp;
                        plot_v_added.value[y_ref_reduced + x_ref_reduced*plot_v_added.height] = v_interp;
                        
                        // Set this is a valid point
                        plot_validpoints.value[y_ref_reduced + x_ref_reduced*plot_validpoints.height] = true;
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
    if (nrhs == 7 && nlhs == 2) {
        // Creat adddisp
        class_adddisp adddisp(plhs,prhs);
        
        // Run analysis
        adddisp.analysis();
    } else {
        mexErrMsgTxt("Incorrect number of inputs.\n");
    }
}
