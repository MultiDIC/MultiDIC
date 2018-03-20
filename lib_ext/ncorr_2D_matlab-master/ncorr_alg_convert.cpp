// This function converts from the old configuration to the new one (i.e. Lagrangian to Eulerian or vice-versa)
// Uses "old" and "new" because it's not certain whether it's lagrangian to eulerian or eulerian to lagrangian

#include <mex.h>
#include <math.h>
#include <queue>
#include <vector>
#include <list>
#include "standard_datatypes.h"
#include "ncorr_datatypes.h"
#include "ncorr_lib.h"

#define DISTANCECUTOFF 0.005
#define GRADNORMCUTOFF 0.00001
#define COUNTERCUTOFF 10

// ----------------------------------------------------//
// Interp plot input ----------------------------------//
// ----------------------------------------------------//

void get_plot_interp(std::vector<class_double_array> &plot_interp,const mxArray *prhs) {    
    // Check input    
    if (mxIsClass(prhs,"cell") && mxGetM(prhs) == 1) {
        for (int i=0; i<(int)mxGetN(prhs); i++) {
            mxArray *mat_plot_interp = mxGetCell(prhs,i);
            if (mat_plot_interp != 0) {
                // Form interpdata
                class_double_array plot_interp_template;
                get_double_array(plot_interp_template,mat_plot_interp);
                
                // Store interp plot
                plot_interp.push_back(plot_interp_template);
            } else {
                mexErrMsgTxt("Some cell contents are empty.\n");
            }
        }
    } else {
        mexErrMsgTxt("Interpdata must be a row vector of class 'cell'.\n");
    }
}

// ----------------------------------------------------//
// convertseedinfo input ------------------------------//
// ----------------------------------------------------//

struct local_struct_convertseedinfo {
    // Constructor
    local_struct_convertseedinfo(){    
        num_region_new = 0;
        num_region_old = 0;
    }
    
    // Properties
    class_double_array paramvector;
    int num_region_new;
    int num_region_old;
};

void get_convertseedinfo(std::vector<local_struct_convertseedinfo> &convertseedinfo,const mxArray *prhs){    
    // Check input    
    if (mxIsClass(prhs,"struct") && mxGetM(prhs) == 1) {
        for (int i=0; i<(int)mxGetN(prhs); i++) {
            mxArray *mat_paramvector = mxGetField(prhs,i,"paramvector");        
            mxArray *mat_num_region_new = mxGetField(prhs,i,"num_region_new");    
            mxArray *mat_num_region_old = mxGetField(prhs,i,"num_region_old");
            if (mat_paramvector != 0 &&
                mat_num_region_new != 0 &&
                mat_num_region_old != 0) {
                // Form convertseedinfo
                local_struct_convertseedinfo convertseedinfo_template;
                get_double_array(convertseedinfo_template.paramvector,mat_paramvector);
                get_integer_scalar(convertseedinfo_template.num_region_new,mat_num_region_new);    
                get_integer_scalar(convertseedinfo_template.num_region_old,mat_num_region_old);    

                // Test sizes
                if (convertseedinfo_template.paramvector.width != 7 || convertseedinfo_template.paramvector.height != 1) {
                    mexErrMsgTxt("'paramvector' is supposed to be a 1x7 vector.\n");
                }
                
                // Store seedinfo
                convertseedinfo.push_back(convertseedinfo_template);
            } else {
                mexErrMsgTxt("Some fields do not exist for convertseedinfo.\n");
            }
        }
    } else {
        mexErrMsgTxt("convertseedinfo must be a row vector of class 'cell'.\n");
    }
}

// ----------------------------------------------------------------------//
// Queue ----------------------------------------------------------------//
// ----------------------------------------------------------------------//

struct comp_queue {
    bool operator ()(std::vector<double> const& a, std::vector<double> const& b) const {
        // sixth element is the distance
        return a[6] > b[6];
    }
};

typedef std::priority_queue<std::vector<double>, std::vector<std::vector<double> >,  comp_queue> heap;

// ----------------------------------------------------//
// Main Class -----------------------------------------//
// ----------------------------------------------------//

class class_convert {
public:
    // Constructor
    class_convert(mxArray *plhs [ ],const mxArray *prhs [ ]);
    
    // Methods:
    void analysis();

private:
    // Properties
    // Inputs:
    std::vector<class_double_array> plot_u_interp_old;                // standard datatype
    std::vector<class_double_array> plot_v_interp_old;                // standard datatype                                       
    std::vector<ncorr_class_roi> roi_old;                             // ncorr datatype
    std::vector<ncorr_class_roi> roi_new;                             // ncorr datatype 
    std::vector<local_struct_convertseedinfo> convertseedinfo;        // local datatype   
    int spacing;                                                      // standard datatype          
    int border_interp;                                                // standard datatype
    int num_img;                                                      // standard datatype
    int total_imgs;                                                   // standard datatype

    // Outputs:
    class_double_array plot_u_new;
    class_double_array plot_v_new;
    class_logical_array plot_validpoints;
    double *outstate;        
        
    // Other variables:  
    class_waitbar waitbar;
    class_double_array plot_calcpoints;

    // Methods -----------------------------------------//
    void analyzepoint(heap &queue,const int &x_new,const int &y_new,const std::vector<double> &paramvector_init,const int &num_region_new,const int &num_region_old);    
    OUT calcpoint(std::vector<double> &paramvector,const int &x_new,const int &y_new,const std::vector<double> &paramvector_init,const int &num_region_old);
    OUT iterativesearch(std::vector<double> &defvector,double &u_old,double &v_old,double &distance,const int &x_new,const int &y_new,const std::vector<double> &defvector_init,const int &num_region_old);
    OUT newton(std::vector<double> &defvector,double &u_old,double &v_old,double &distance,double &gradnorm,const int &x_new,const int &y_new,const std::vector<double> &defvector_init,const int &num_region_old);
    OUT interpqbs_convert(std::vector<double> &interpvector,const std::vector<double> &defvector,const int &num_region_old);
};

class_convert::class_convert(mxArray *plhs[ ],const mxArray *prhs[ ]) {
    // Get inputs -------------------------------------------------//
    // input 1: plot_u_interp_old
    get_plot_interp(plot_u_interp_old,prhs[0]);
    // input 2: plot_v_interp_old
    get_plot_interp(plot_v_interp_old,prhs[1]);
    // input 3: roi_old
    get_rois(roi_old,prhs[2]);
    // input 4: roi_new
    get_rois(roi_new,prhs[3]);
    // input 5: convertseedinfo
    get_convertseedinfo(convertseedinfo,prhs[4]);
    // input 6: spacing    
    get_integer_scalar(spacing,prhs[5]);
    // input 7: border_interp
    get_integer_scalar(border_interp,prhs[6]);
    // input 8: num_img    
    get_integer_scalar(num_img,prhs[7]);
    // input 9: total_imgs
    get_integer_scalar(total_imgs,prhs[8]);

    // Check inputs - check sizes - UPDATE THIS LATER!!!!!   
    
    // Calc points -----------------------------------------------//
    plot_calcpoints.alloc(roi_new[0].mask.height,roi_new[0].mask.width);
        
    // Form/get outputs ------------------------------------------//
    // output 1: plot_new
    // Form deformation structure
    mwSize def_dims[2] = {1,1};
    int def_numfields = 3;
    const char *def_fieldnames[] = {"plot_u_new","plot_v_new","plot_validpoints"};
    plhs[0] = mxCreateStructArray(2, def_dims, def_numfields, def_fieldnames);

    // Form fields
    mxArray *mat_plot_u_new = mxCreateDoubleMatrix(roi_new[0].mask.height, roi_new[0].mask.width, mxREAL);
    mxArray *mat_plot_v_new = mxCreateDoubleMatrix(roi_new[0].mask.height, roi_new[0].mask.width, mxREAL);
    mxArray *mat_plot_validpoints = mxCreateLogicalMatrix(roi_new[0].mask.height,roi_new[0].mask.width);

    // Add fields to structure
    // add u:
    mxSetFieldByNumber(plhs[0],0,0,mat_plot_u_new);
    // add v:
    mxSetFieldByNumber(plhs[0],0,1,mat_plot_v_new);
    // add valid points:
    mxSetFieldByNumber(plhs[0],0,2,mat_plot_validpoints);

    // output 2: outstate
    plhs[1] = mxCreateDoubleMatrix(1,1,mxREAL);
    
    // Get outputs -----------------------------------------------//
    // output 1: plot_new
    // u_new:
    get_double_array(plot_u_new,mat_plot_u_new);
    // v_new:
    get_double_array(plot_v_new,mat_plot_v_new);
    // plot_validpoints:
    get_logical_array(plot_validpoints,mat_plot_validpoints);
    // output 2: outstate
    outstate = mxGetPr(plhs[1]);
}

// ----------------------------------------------------//
// Main Class Methods ---------------------------------//
// ----------------------------------------------------//

void class_convert::analysis() {
    // Initialize outstate to cancelled
    *outstate = (double)CANCELLED;

    // Set up waitbar ----------------------------------------------//
    int computepoints = 0;
    for (int i=0; i<(int)roi_new[0].region.size(); i++) {
        computepoints += roi_new[0].region[i].totalpoints;
    }
    waitbar.start(num_img,total_imgs,computepoints);   
    
    // Begin Analysis ----------------------------------------------//
    // Iterate over the seeds
    for (int i=0; i<(int)convertseedinfo.size(); i++) {
        // Initialize queue
        heap queue;

        // Add seed to queue 
        std::vector<double> paramvector_seed(7,0); // [x_new y_new x_old y_old u_old v_old distance]
        paramvector_seed[0] = convertseedinfo[i].paramvector.value[0];
        paramvector_seed[1] = convertseedinfo[i].paramvector.value[1];
        paramvector_seed[2] = convertseedinfo[i].paramvector.value[2];
        paramvector_seed[3] = convertseedinfo[i].paramvector.value[3];
        paramvector_seed[4] = convertseedinfo[i].paramvector.value[4];
        paramvector_seed[5] = convertseedinfo[i].paramvector.value[5];
        paramvector_seed[6] = convertseedinfo[i].paramvector.value[6];
        queue.push(paramvector_seed);

        // Inactivate seed point
        plot_calcpoints.value[(int)paramvector_seed[1]+(int)paramvector_seed[0]*plot_calcpoints.height] = true; 
        plot_validpoints.value[(int)paramvector_seed[1]+(int)paramvector_seed[0]*plot_validpoints.height] = true; 

        // Enter While Loop - Exit when queue is empty
        while (!queue.empty()) {                        
            // 1) Load point from queue
            // 2) Delete point from queue
            // 3) Add data to plots 
            // 4) Analyze four surrounding points
            
            // Step 1: load 
            std::vector<double> paramvector_init = queue.top();

            // Step 2: delete
            queue.pop();

            // Step 3: add data to plots
            plot_u_new.value[(int)paramvector_init[1]+(int)paramvector_init[0]*plot_u_new.height] = -paramvector_init[4];
            plot_v_new.value[(int)paramvector_init[1]+(int)paramvector_init[0]*plot_v_new.height] = -paramvector_init[5];
            
            // Step 4: analyze four surrounding points
            analyzepoint(queue,(int)paramvector_init[0],(int)paramvector_init[1]-1,paramvector_init,convertseedinfo[i].num_region_new,convertseedinfo[i].num_region_old);  
            analyzepoint(queue,(int)paramvector_init[0]+1,(int)paramvector_init[1],paramvector_init,convertseedinfo[i].num_region_new,convertseedinfo[i].num_region_old); 
            analyzepoint(queue,(int)paramvector_init[0],(int)paramvector_init[1]+1,paramvector_init,convertseedinfo[i].num_region_new,convertseedinfo[i].num_region_old); 
            analyzepoint(queue,(int)paramvector_init[0]-1,(int)paramvector_init[1],paramvector_init,convertseedinfo[i].num_region_new,convertseedinfo[i].num_region_old);               
        
            // Update and check waitbar
            if (!waitbar.updateandcheck()) {
                // Waitbar was cancelled
                return;
            }
        }
    }    
    
    // At this point analysis has been completed successfully
    *outstate = SUCCESS; 
}

void class_convert::analyzepoint(heap &queue,const int &x_new,const int &y_new,const std::vector<double> &paramvector_init,const int &num_region_new,const int &num_region_old) {
    // Make sure point is within region bounds first
    if (x_new >= roi_new[0].region[num_region_new].leftbound &&
        x_new <= roi_new[0].region[num_region_new].rightbound &&
        y_new >= roi_new[0].region[num_region_new].upperbound &&
        y_new <= roi_new[0].region[num_region_new].lowerbound &&        
        !plot_calcpoints.value[(int)y_new+(int)x_new*plot_calcpoints.height] && 
        roi_new[0].withinregion((int)x_new,(int)y_new,num_region_new)) {
        // Initialize paramvector
        std::vector<double> paramvector(7,0); // [x_new y_new x_old y_old u_old v_old distance]
        
        // Calculate paramvector for point
        OUT outstate_calcpoint = calcpoint(paramvector,x_new,y_new,paramvector_init,num_region_old);
        
        // Make sure parameters are correct before adding them to queue
        if (outstate_calcpoint == SUCCESS && 
            paramvector[6] < DISTANCECUTOFF &&
            (int)ncorr_round(paramvector[2]) >= roi_old[0].region[num_region_old].leftbound &&
            (int)ncorr_round(paramvector[2]) <= roi_old[0].region[num_region_old].rightbound &&
            (int)ncorr_round(paramvector[3]) >= roi_old[0].region[num_region_old].upperbound &&
            (int)ncorr_round(paramvector[3]) <= roi_old[0].region[num_region_old].lowerbound &&
            roi_old[0].mask.value[(int)ncorr_round(paramvector[3])+(int)ncorr_round(paramvector[2])*roi_old[0].mask.height]) {
            // Insert paramvector based on distance
            queue.push(paramvector);
            
            // Valid point
            plot_validpoints.value[y_new+x_new*plot_validpoints.height] = true;        
        }        
        // Calculated point
        plot_calcpoints.value[y_new+x_new*plot_calcpoints.height] = true;
        
        // Increment the waitbar
        waitbar.increment();
    }
}

OUT class_convert::calcpoint(std::vector<double> &paramvector,const int &x_new,const int &y_new,const std::vector<double> &paramvector_init,const int &num_region_old) {
    // Find initial guess -> Refine results with FA-GN -> Return true or false and store output
    
    // Step 1: Get initial guess    
    std::vector<double> defvector_init(2,0); // [x_old y_old]
    defvector_init[0] = paramvector_init[2]+(x_new-paramvector_init[0]);
    defvector_init[1] = paramvector_init[3]+(y_new-paramvector_init[1]);
    
    // Step 2: Perform Iterative refinement
    std::vector<double> defvector(2,0);
    double u_old;
    double v_old;
    double distance;
    OUT outstate_iterative = iterativesearch(defvector,u_old,v_old,distance,x_new,y_new,defvector_init,num_region_old);

    if (outstate_iterative == SUCCESS) {
        // Step 3: Store output and return success
        paramvector[0] = x_new;
        paramvector[1] = y_new;
        paramvector[2] = defvector[0];
        paramvector[3] = defvector[1];
        paramvector[4] = u_old;
        paramvector[5] = v_old;
        paramvector[6] = distance;
        return SUCCESS;
    }    
    return FAILED;
}

OUT class_convert::iterativesearch(std::vector<double> &defvector,double &u_old,double &v_old,double &distance,const int &x_new,const int &y_new,const std::vector<double> &defvector_init,const int &num_region_old) {
    double gradnorm; // norm of gradient vector    
    
    // Start iterations - For first iteration use defvector_init
    OUT outstate_newton = newton(defvector,u_old,v_old,distance,gradnorm,x_new,y_new,defvector_init,num_region_old);
    
    // Initialize counter
    int counter = 1;        
    while (outstate_newton == SUCCESS && gradnorm > GRADNORMCUTOFF && counter < COUNTERCUTOFF) {
        // For rest of iterations use defvector from previous iterations
        outstate_newton = newton(defvector,u_old,v_old,distance,gradnorm,x_new,y_new,defvector,num_region_old);
        ++counter;
    }    
    
    if (outstate_newton == SUCCESS) {
        return SUCCESS;
    }
    return FAILED;
}

OUT class_convert::newton(std::vector<double> &defvector,double &u_old,double &v_old,double &distance,double &gradnorm,const int &x_new,const int &y_new,const std::vector<double> &defvector_init,const int &num_region_old) {        
    // Note that static memory is not safe for multithreading
    static double gradient_buffer[2];
    static double hessian_buffer[4];
    
    // Interpolate values- answers are stored in interpvector. Note that I use interpqbs_convert it saves some
    // calculations by combining the interpolations. In addition, they require
    // gradients to be interpolated as well.
    std::vector<double> interpvector(6,0); // [u v du/dx du/dy dv/dx dv/dy]    
    OUT outstate_interp = interpqbs_convert(interpvector,defvector_init,num_region_old);
    
    if (outstate_interp == SUCCESS) {
        // Calculate Gradient
        gradient_buffer[0] = -2.0*(((double)x_new-(defvector_init[0]+interpvector[0]/((double)spacing+1.0)))*(1+interpvector[2]/((double)spacing+1.0))+((double)y_new-(defvector_init[1]+interpvector[1]/((double)spacing+1.0)))*(interpvector[4]/((double)spacing+1.0)));
        gradient_buffer[1] = -2.0*(((double)x_new-(defvector_init[0]+interpvector[0]/((double)spacing+1.0)))*(interpvector[3]/((double)spacing+1.0))+((double)y_new-(defvector_init[1]+interpvector[1]/((double)spacing+1.0)))*(1+interpvector[5]/((double)spacing+1.0)));
        
        // Calculate Hessian
        hessian_buffer[0] = 2.0*(pow(1+interpvector[2]/((double)spacing+1.0),2)+pow(interpvector[4]/((double)spacing+1.0),2));
        hessian_buffer[1] = 2.0*((interpvector[3]/((double)spacing+1.0))*(1+interpvector[2]/((double)spacing+1.0))+(1+interpvector[5]/((double)spacing+1.0))*(interpvector[4]/((double)spacing+1.0)));
        hessian_buffer[2] = hessian_buffer[1];
        hessian_buffer[3] = 2.0*(pow(interpvector[3]/((double)spacing+1.0),2)+pow(1+interpvector[5]/((double)spacing+1.0),2));
        
        // Determine new current coordinates - make sure hessian is positive definite
        // From :http://www.math.northwestern.edu/~clark/285/2006-07/handouts/pos-def.pdf
        // Make sure det(hess) > 0 and hess(1,1) > 0 
        double det_hess = hessian_buffer[0]*hessian_buffer[3]-hessian_buffer[2]*hessian_buffer[1];
        if (det_hess > 0 && hessian_buffer[0] > 0) { 
            // Solve for new coordinates
            defvector[0] = defvector_init[0]-(1.0/det_hess)*(hessian_buffer[3]*gradient_buffer[0]-hessian_buffer[2]*gradient_buffer[1]);
            defvector[1] = defvector_init[1]-(1.0/det_hess)*(-hessian_buffer[1]*gradient_buffer[0]+hessian_buffer[0]*gradient_buffer[1]);

            // Determine u and v for new coordinates
            outstate_interp = interpqbs_convert(interpvector,defvector,num_region_old);    
            if (outstate_interp == SUCCESS) {
                u_old = interpvector[0];
                v_old = interpvector[1];
                distance = sqrt(pow((double)x_new-(defvector[0]+u_old/((double)spacing+1.0)),2)+pow((double)y_new-(defvector[1]+v_old/((double)spacing+1.0)),2));
                gradnorm = sqrt(pow(gradient_buffer[0],2)+pow(gradient_buffer[1],2));
                return SUCCESS;
            }
        }
    }
    return FAILED;
}

OUT class_convert::interpqbs_convert(std::vector<double> &interpvector,const std::vector<double> &defvector,const int &num_region_old) {
    // Note that static memory is not safe for multithreading
    // interpvector = [u v du/dx du/dy dv/dx dv/dy]
    static double x_vec_buffer[6];
    static double x_vec_dx_buffer[6];
    static double x_vec_star_buffer[6];
    static double y_vec_buffer[6];
    static double y_vec_dy_buffer[6];
    static double y_vec_star_buffer[6];
    static double QK_B_QKT_buffer[36];
    
    double x_tilda = defvector[0];
    double y_tilda = defvector[1];
    
    int x_tilda_floor = (int)floor(x_tilda);
    int y_tilda_floor = (int)floor(y_tilda);
    
    // Get bounds of the desired b-spline coefficients used for interpolation
    int top = y_tilda_floor-roi_old[0].region[num_region_old].upperbound+border_interp-2;
    int left = x_tilda_floor-roi_old[0].region[num_region_old].leftbound+border_interp-2;
    int bottom = y_tilda_floor-roi_old[0].region[num_region_old].upperbound+border_interp+3;
    int right = x_tilda_floor-roi_old[0].region[num_region_old].leftbound+border_interp+3;
        
    if (top >= 0 &&
        left >= 0 && 
        bottom < plot_u_interp_old[num_region_old].height &&
        right < plot_u_interp_old[num_region_old].width) {            
        double x_tilda_delta = x_tilda-(double)x_tilda_floor;
        double y_tilda_delta = y_tilda-(double)y_tilda_floor;

        // Form x_vec, y_vec, x_vec_dx, and y_vec_dy
        // x_vec
        x_vec_buffer[0] = 1.0;
        x_vec_buffer[1] = x_tilda_delta;
        x_vec_buffer[2] = x_vec_buffer[1]*x_tilda_delta;
        x_vec_buffer[3] = x_vec_buffer[2]*x_tilda_delta;
        x_vec_buffer[4] = x_vec_buffer[3]*x_tilda_delta;
        x_vec_buffer[5] = x_vec_buffer[4]*x_tilda_delta;

        // y_vec
        y_vec_buffer[0] = 1.0;
        y_vec_buffer[1] = y_tilda_delta;
        y_vec_buffer[2] = y_vec_buffer[1]*y_tilda_delta;
        y_vec_buffer[3] = y_vec_buffer[2]*y_tilda_delta;
        y_vec_buffer[4] = y_vec_buffer[3]*y_tilda_delta;
        y_vec_buffer[5] = y_vec_buffer[4]*y_tilda_delta;
        
        // x_vec_dx
        x_vec_dx_buffer[0] = 0.0;
        x_vec_dx_buffer[1] = 1.0;
        x_vec_dx_buffer[2] = 2.0*x_vec_buffer[1];
        x_vec_dx_buffer[3] = 3.0*x_vec_buffer[2];
        x_vec_dx_buffer[4] = 4.0*x_vec_buffer[3];
        x_vec_dx_buffer[5] = 5.0*x_vec_buffer[4];

        // y_vec_dx
        y_vec_dy_buffer[0] = 0.0;
        y_vec_dy_buffer[1] = 1.0;
        y_vec_dy_buffer[2] = 2.0*y_vec_buffer[1];
        y_vec_dy_buffer[3] = 3.0*y_vec_buffer[2];
        y_vec_dy_buffer[4] = 4.0*y_vec_buffer[3];
        y_vec_dy_buffer[5] = 5.0*y_vec_buffer[4];
        

        // Compute QK*B_coef*QK^T for u-displacements first
        double b0 = plot_u_interp_old[num_region_old].value[(top)+(left)*plot_u_interp_old[num_region_old].height];
        double b1 = plot_u_interp_old[num_region_old].value[(top+1)+(left)*plot_u_interp_old[num_region_old].height];
        double b2 = plot_u_interp_old[num_region_old].value[(top+2)+(left)*plot_u_interp_old[num_region_old].height];
        double b3 = plot_u_interp_old[num_region_old].value[(top+3)+(left)*plot_u_interp_old[num_region_old].height];
        double b4 = plot_u_interp_old[num_region_old].value[(top+4)+(left)*plot_u_interp_old[num_region_old].height];
        double b5 = plot_u_interp_old[num_region_old].value[(top+5)+(left)*plot_u_interp_old[num_region_old].height];
        double b6 = plot_u_interp_old[num_region_old].value[(top)+(left+1)*plot_u_interp_old[num_region_old].height];
        double b7 = plot_u_interp_old[num_region_old].value[(top+1)+(left+1)*plot_u_interp_old[num_region_old].height];
        double b8 = plot_u_interp_old[num_region_old].value[(top+2)+(left+1)*plot_u_interp_old[num_region_old].height];
        double b9 = plot_u_interp_old[num_region_old].value[(top+3)+(left+1)*plot_u_interp_old[num_region_old].height];
        double b10 = plot_u_interp_old[num_region_old].value[(top+4)+(left+1)*plot_u_interp_old[num_region_old].height];
        double b11 = plot_u_interp_old[num_region_old].value[(top+5)+(left+1)*plot_u_interp_old[num_region_old].height];
        double b12 = plot_u_interp_old[num_region_old].value[(top)+(left+2)*plot_u_interp_old[num_region_old].height];
        double b13 = plot_u_interp_old[num_region_old].value[(top+1)+(left+2)*plot_u_interp_old[num_region_old].height];
        double b14 = plot_u_interp_old[num_region_old].value[(top+2)+(left+2)*plot_u_interp_old[num_region_old].height];
        double b15 = plot_u_interp_old[num_region_old].value[(top+3)+(left+2)*plot_u_interp_old[num_region_old].height];
        double b16 = plot_u_interp_old[num_region_old].value[(top+4)+(left+2)*plot_u_interp_old[num_region_old].height];
        double b17 = plot_u_interp_old[num_region_old].value[(top+5)+(left+2)*plot_u_interp_old[num_region_old].height];
        double b18 = plot_u_interp_old[num_region_old].value[(top)+(left+3)*plot_u_interp_old[num_region_old].height];
        double b19 = plot_u_interp_old[num_region_old].value[(top+1)+(left+3)*plot_u_interp_old[num_region_old].height];
        double b20 = plot_u_interp_old[num_region_old].value[(top+2)+(left+3)*plot_u_interp_old[num_region_old].height];
        double b21 = plot_u_interp_old[num_region_old].value[(top+3)+(left+3)*plot_u_interp_old[num_region_old].height];
        double b22 = plot_u_interp_old[num_region_old].value[(top+4)+(left+3)*plot_u_interp_old[num_region_old].height];
        double b23 = plot_u_interp_old[num_region_old].value[(top+5)+(left+3)*plot_u_interp_old[num_region_old].height];
        double b24 = plot_u_interp_old[num_region_old].value[(top)+(left+4)*plot_u_interp_old[num_region_old].height];
        double b25 = plot_u_interp_old[num_region_old].value[(top+1)+(left+4)*plot_u_interp_old[num_region_old].height];
        double b26 = plot_u_interp_old[num_region_old].value[(top+2)+(left+4)*plot_u_interp_old[num_region_old].height];
        double b27 = plot_u_interp_old[num_region_old].value[(top+3)+(left+4)*plot_u_interp_old[num_region_old].height];
        double b28 = plot_u_interp_old[num_region_old].value[(top+4)+(left+4)*plot_u_interp_old[num_region_old].height];
        double b29 = plot_u_interp_old[num_region_old].value[(top+5)+(left+4)*plot_u_interp_old[num_region_old].height];
        double b30 = plot_u_interp_old[num_region_old].value[(top)+(left+5)*plot_u_interp_old[num_region_old].height];
        double b31 = plot_u_interp_old[num_region_old].value[(top+1)+(left+5)*plot_u_interp_old[num_region_old].height];
        double b32 = plot_u_interp_old[num_region_old].value[(top+2)+(left+5)*plot_u_interp_old[num_region_old].height];
        double b33 = plot_u_interp_old[num_region_old].value[(top+3)+(left+5)*plot_u_interp_old[num_region_old].height];
        double b34 = plot_u_interp_old[num_region_old].value[(top+4)+(left+5)*plot_u_interp_old[num_region_old].height];
        double b35 = plot_u_interp_old[num_region_old].value[(top+5)+(left+5)*plot_u_interp_old[num_region_old].height];

        QK_B_QKT_buffer[0] = 0.00006944444444444444*b0+0.001805555555555556*b1+0.001805555555555556*b10+0.004583333333333333*b12+0.1191666666666667*b13+0.3025*b14+0.1191666666666667*b15+0.004583333333333333*b16+0.001805555555555556*b18+0.04694444444444444*b19+0.004583333333333333*b2+0.1191666666666667*b20+0.04694444444444444*b21+0.001805555555555556*b22+0.00006944444444444444*b24+0.001805555555555556*b25+0.004583333333333333*b26+0.001805555555555556*b27+0.00006944444444444444*b28+0.001805555555555556*b3+0.00006944444444444444*b4+0.001805555555555556*b6+0.04694444444444444*b7+0.1191666666666667*b8+0.04694444444444444*b9;
        QK_B_QKT_buffer[1] = 0.009027777777777778*b10-0.003472222222222222*b1-0.0003472222222222222*b0-0.02291666666666667*b12-0.2291666666666667*b13+0.2291666666666667*b15+0.02291666666666667*b16-0.009027777777777778*b18-0.09027777777777778*b19+0.09027777777777778*b21+0.009027777777777778*b22-0.0003472222222222222*b24-0.003472222222222222*b25+0.003472222222222222*b27+0.0003472222222222222*b28+0.003472222222222222*b3+0.0003472222222222222*b4-0.009027777777777778*b6-0.09027777777777778*b7+0.09027777777777778*b9;
        QK_B_QKT_buffer[2] = 0.0006944444444444444*b0+0.001388888888888889*b1+0.01805555555555556*b10+0.04583333333333333*b12+0.09166666666666667*b13-0.275*b14+0.09166666666666667*b15+0.04583333333333333*b16+0.01805555555555556*b18+0.03611111111111111*b19-0.004166666666666667*b2-0.1083333333333333*b20+0.03611111111111111*b21+0.01805555555555556*b22+0.0006944444444444444*b24+0.001388888888888889*b25-0.004166666666666667*b26+0.001388888888888889*b27+0.0006944444444444444*b28+0.001388888888888889*b3+0.0006944444444444444*b4+0.01805555555555556*b6+0.03611111111111111*b7-0.1083333333333333*b8+0.03611111111111111*b9;
        QK_B_QKT_buffer[3] = 0.001388888888888889*b1-0.0006944444444444444*b0+0.01805555555555556*b10-0.04583333333333333*b12+0.09166666666666667*b13-0.09166666666666667*b15+0.04583333333333333*b16-0.01805555555555556*b18+0.03611111111111111*b19-0.03611111111111111*b21+0.01805555555555556*b22-0.0006944444444444444*b24+0.001388888888888889*b25-0.001388888888888889*b27+0.0006944444444444444*b28-0.001388888888888889*b3+0.0006944444444444444*b4-0.01805555555555556*b6+0.03611111111111111*b7-0.03611111111111111*b9;
        QK_B_QKT_buffer[4] = 0.0003472222222222222*b0-0.001388888888888889*b1+0.009027777777777778*b10+0.02291666666666667*b12-0.09166666666666667*b13+0.1375*b14-0.09166666666666667*b15+0.02291666666666667*b16+0.009027777777777778*b18-0.03611111111111111*b19+0.002083333333333333*b2+0.05416666666666667*b20-0.03611111111111111*b21+0.009027777777777778*b22+0.0003472222222222222*b24-0.001388888888888889*b25+0.002083333333333333*b26-0.001388888888888889*b27+0.0003472222222222222*b28-0.001388888888888889*b3+0.0003472222222222222*b4+0.009027777777777778*b6-0.03611111111111111*b7+0.05416666666666667*b8-0.03611111111111111*b9;
        QK_B_QKT_buffer[5] = 0.0003472222222222222*b1-0.00006944444444444444*b0-0.009027777777777778*b10+0.001805555555555556*b11-0.004583333333333333*b12+0.02291666666666667*b13-0.04583333333333333*b14+0.04583333333333333*b15-0.02291666666666667*b16+0.004583333333333333*b17-0.001805555555555556*b18+0.009027777777777778*b19-0.0006944444444444444*b2-0.01805555555555556*b20+0.01805555555555556*b21-0.009027777777777778*b22+0.001805555555555556*b23-0.00006944444444444444*b24+0.0003472222222222222*b25-0.0006944444444444444*b26+0.0006944444444444444*b27-0.0003472222222222222*b28+0.00006944444444444444*b29+0.0006944444444444444*b3-0.0003472222222222222*b4+0.00006944444444444444*b5-0.001805555555555556*b6+0.009027777777777778*b7-0.01805555555555556*b8+0.01805555555555556*b9;
        QK_B_QKT_buffer[6] = 0.003472222222222222*b18-0.009027777777777778*b1-0.003472222222222222*b10-0.0003472222222222222*b0+0.09027777777777778*b19-0.02291666666666667*b2+0.2291666666666667*b20+0.09027777777777778*b21+0.003472222222222222*b22+0.0003472222222222222*b24+0.009027777777777778*b25+0.02291666666666667*b26+0.009027777777777778*b27+0.0003472222222222222*b28-0.009027777777777778*b3-0.0003472222222222222*b4-0.003472222222222222*b6-0.09027777777777778*b7-0.2291666666666667*b8-0.09027777777777778*b9;
        QK_B_QKT_buffer[7] = 0.001736111111111111*b0+0.01736111111111111*b1-0.01736111111111111*b10-0.01736111111111111*b18-0.1736111111111111*b19+0.1736111111111111*b21+0.01736111111111111*b22-0.001736111111111111*b24-0.01736111111111111*b25+0.01736111111111111*b27+0.001736111111111111*b28-0.01736111111111111*b3-0.001736111111111111*b4+0.01736111111111111*b6+0.1736111111111111*b7-0.1736111111111111*b9;
        QK_B_QKT_buffer[8] = 0.03472222222222222*b18-0.006944444444444444*b1-0.03472222222222222*b10-0.003472222222222222*b0+0.06944444444444444*b19+0.02083333333333333*b2-0.2083333333333333*b20+0.06944444444444444*b21+0.03472222222222222*b22+0.003472222222222222*b24+0.006944444444444444*b25-0.02083333333333333*b26+0.006944444444444444*b27+0.003472222222222222*b28-0.006944444444444444*b3-0.003472222222222222*b4-0.03472222222222222*b6-0.06944444444444444*b7+0.2083333333333333*b8-0.06944444444444444*b9;
        QK_B_QKT_buffer[9] = 0.003472222222222222*b0-0.006944444444444444*b1-0.03472222222222222*b10-0.03472222222222222*b18+0.06944444444444444*b19-0.06944444444444444*b21+0.03472222222222222*b22-0.003472222222222222*b24+0.006944444444444444*b25-0.006944444444444444*b27+0.003472222222222222*b28+0.006944444444444444*b3-0.003472222222222222*b4+0.03472222222222222*b6-0.06944444444444444*b7+0.06944444444444444*b9;
        QK_B_QKT_buffer[10] = 0.006944444444444444*b1-0.001736111111111111*b0-0.01736111111111111*b10+0.01736111111111111*b18-0.06944444444444444*b19-0.01041666666666667*b2+0.1041666666666667*b20-0.06944444444444444*b21+0.01736111111111111*b22+0.001736111111111111*b24-0.006944444444444444*b25+0.01041666666666667*b26-0.006944444444444444*b27+0.001736111111111111*b28+0.006944444444444444*b3-0.001736111111111111*b4-0.01736111111111111*b6+0.06944444444444444*b7-0.1041666666666667*b8+0.06944444444444444*b9;
        QK_B_QKT_buffer[11] = 0.0003472222222222222*b0-0.001736111111111111*b1+0.01736111111111111*b10-0.003472222222222222*b11-0.003472222222222222*b18+0.01736111111111111*b19+0.003472222222222222*b2-0.03472222222222222*b20+0.03472222222222222*b21-0.01736111111111111*b22+0.003472222222222222*b23-0.0003472222222222222*b24+0.001736111111111111*b25-0.003472222222222222*b26+0.003472222222222222*b27-0.001736111111111111*b28+0.0003472222222222222*b29-0.003472222222222222*b3+0.001736111111111111*b4-0.0003472222222222222*b5+0.003472222222222222*b6-0.01736111111111111*b7+0.03472222222222222*b8-0.03472222222222222*b9;
        QK_B_QKT_buffer[12] = 0.0006944444444444444*b0+0.01805555555555556*b1+0.001388888888888889*b10-0.004166666666666667*b12-0.1083333333333333*b13-0.275*b14-0.1083333333333333*b15-0.004166666666666667*b16+0.001388888888888889*b18+0.03611111111111111*b19+0.04583333333333333*b2+0.09166666666666667*b20+0.03611111111111111*b21+0.001388888888888889*b22+0.0006944444444444444*b24+0.01805555555555556*b25+0.04583333333333333*b26+0.01805555555555556*b27+0.0006944444444444444*b28+0.01805555555555556*b3+0.0006944444444444444*b4+0.001388888888888889*b6+0.03611111111111111*b7+0.09166666666666667*b8+0.03611111111111111*b9;
        QK_B_QKT_buffer[13] = 0.006944444444444444*b10-0.03472222222222222*b1-0.003472222222222222*b0+0.02083333333333333*b12+0.2083333333333333*b13-0.2083333333333333*b15-0.02083333333333333*b16-0.006944444444444444*b18-0.06944444444444444*b19+0.06944444444444444*b21+0.006944444444444444*b22-0.003472222222222222*b24-0.03472222222222222*b25+0.03472222222222222*b27+0.003472222222222222*b28+0.03472222222222222*b3+0.003472222222222222*b4-0.006944444444444444*b6-0.06944444444444444*b7+0.06944444444444444*b9;
        QK_B_QKT_buffer[14] = 0.006944444444444444*b0+0.01388888888888889*b1+0.01388888888888889*b10-0.04166666666666667*b12-0.08333333333333333*b13+0.25*b14-0.08333333333333333*b15-0.04166666666666667*b16+0.01388888888888889*b18+0.02777777777777778*b19-0.04166666666666667*b2-0.08333333333333333*b20+0.02777777777777778*b21+0.01388888888888889*b22+0.006944444444444444*b24+0.01388888888888889*b25-0.04166666666666667*b26+0.01388888888888889*b27+0.006944444444444444*b28+0.01388888888888889*b3+0.006944444444444444*b4+0.01388888888888889*b6+0.02777777777777778*b7-0.08333333333333333*b8+0.02777777777777778*b9;
        QK_B_QKT_buffer[15] = 0.01388888888888889*b1-0.006944444444444444*b0+0.01388888888888889*b10+0.04166666666666667*b12-0.08333333333333333*b13+0.08333333333333333*b15-0.04166666666666667*b16-0.01388888888888889*b18+0.02777777777777778*b19-0.02777777777777778*b21+0.01388888888888889*b22-0.006944444444444444*b24+0.01388888888888889*b25-0.01388888888888889*b27+0.006944444444444444*b28-0.01388888888888889*b3+0.006944444444444444*b4-0.01388888888888889*b6+0.02777777777777778*b7-0.02777777777777778*b9;
        QK_B_QKT_buffer[16] = 0.003472222222222222*b0-0.01388888888888889*b1+0.006944444444444444*b10-0.02083333333333333*b12+0.08333333333333333*b13-0.125*b14+0.08333333333333333*b15-0.02083333333333333*b16+0.006944444444444444*b18-0.02777777777777778*b19+0.02083333333333333*b2+0.04166666666666667*b20-0.02777777777777778*b21+0.006944444444444444*b22+0.003472222222222222*b24-0.01388888888888889*b25+0.02083333333333333*b26-0.01388888888888889*b27+0.003472222222222222*b28-0.01388888888888889*b3+0.003472222222222222*b4+0.006944444444444444*b6-0.02777777777777778*b7+0.04166666666666667*b8-0.02777777777777778*b9;
        QK_B_QKT_buffer[17] = 0.003472222222222222*b1-0.0006944444444444444*b0-0.006944444444444444*b10+0.001388888888888889*b11+0.004166666666666667*b12-0.02083333333333333*b13+0.04166666666666667*b14-0.04166666666666667*b15+0.02083333333333333*b16-0.004166666666666667*b17-0.001388888888888889*b18+0.006944444444444444*b19-0.006944444444444444*b2-0.01388888888888889*b20+0.01388888888888889*b21-0.006944444444444444*b22+0.001388888888888889*b23-0.0006944444444444444*b24+0.003472222222222222*b25-0.006944444444444444*b26+0.006944444444444444*b27-0.003472222222222222*b28+0.0006944444444444444*b29+0.006944444444444444*b3-0.003472222222222222*b4+0.0006944444444444444*b5-0.001388888888888889*b6+0.006944444444444444*b7-0.01388888888888889*b8+0.01388888888888889*b9;
        QK_B_QKT_buffer[18] = 0.001388888888888889*b10-0.01805555555555556*b1-0.0006944444444444444*b0-0.001388888888888889*b18-0.03611111111111111*b19-0.04583333333333333*b2-0.09166666666666667*b20-0.03611111111111111*b21-0.001388888888888889*b22+0.0006944444444444444*b24+0.01805555555555556*b25+0.04583333333333333*b26+0.01805555555555556*b27+0.0006944444444444444*b28-0.01805555555555556*b3-0.0006944444444444444*b4+0.001388888888888889*b6+0.03611111111111111*b7+0.09166666666666667*b8+0.03611111111111111*b9;
        QK_B_QKT_buffer[19] = 0.003472222222222222*b0+0.03472222222222222*b1+0.006944444444444444*b10+0.006944444444444444*b18+0.06944444444444444*b19-0.06944444444444444*b21-0.006944444444444444*b22-0.003472222222222222*b24-0.03472222222222222*b25+0.03472222222222222*b27+0.003472222222222222*b28-0.03472222222222222*b3-0.003472222222222222*b4-0.006944444444444444*b6-0.06944444444444444*b7+0.06944444444444444*b9;
        QK_B_QKT_buffer[20] = 0.01388888888888889*b10-0.01388888888888889*b1-0.006944444444444444*b0-0.01388888888888889*b18-0.02777777777777778*b19+0.04166666666666667*b2+0.08333333333333333*b20-0.02777777777777778*b21-0.01388888888888889*b22+0.006944444444444444*b24+0.01388888888888889*b25-0.04166666666666667*b26+0.01388888888888889*b27+0.006944444444444444*b28-0.01388888888888889*b3-0.006944444444444444*b4+0.01388888888888889*b6+0.02777777777777778*b7-0.08333333333333333*b8+0.02777777777777778*b9;
        QK_B_QKT_buffer[21] = 0.006944444444444444*b0-0.01388888888888889*b1+0.01388888888888889*b10+0.01388888888888889*b18-0.02777777777777778*b19+0.02777777777777778*b21-0.01388888888888889*b22-0.006944444444444444*b24+0.01388888888888889*b25-0.01388888888888889*b27+0.006944444444444444*b28+0.01388888888888889*b3-0.006944444444444444*b4-0.01388888888888889*b6+0.02777777777777778*b7-0.02777777777777778*b9;
        QK_B_QKT_buffer[22] = 0.01388888888888889*b1-0.003472222222222222*b0+0.006944444444444444*b10-0.006944444444444444*b18+0.02777777777777778*b19-0.02083333333333333*b2-0.04166666666666667*b20+0.02777777777777778*b21-0.006944444444444444*b22+0.003472222222222222*b24-0.01388888888888889*b25+0.02083333333333333*b26-0.01388888888888889*b27+0.003472222222222222*b28+0.01388888888888889*b3-0.003472222222222222*b4+0.006944444444444444*b6-0.02777777777777778*b7+0.04166666666666667*b8-0.02777777777777778*b9;
        QK_B_QKT_buffer[23] = 0.0006944444444444444*b0-0.003472222222222222*b1-0.006944444444444444*b10+0.001388888888888889*b11+0.001388888888888889*b18-0.006944444444444444*b19+0.006944444444444444*b2+0.01388888888888889*b20-0.01388888888888889*b21+0.006944444444444444*b22-0.001388888888888889*b23-0.0006944444444444444*b24+0.003472222222222222*b25-0.006944444444444444*b26+0.006944444444444444*b27-0.003472222222222222*b28+0.0006944444444444444*b29-0.006944444444444444*b3+0.003472222222222222*b4-0.0006944444444444444*b5-0.001388888888888889*b6+0.006944444444444444*b7-0.01388888888888889*b8+0.01388888888888889*b9;
        QK_B_QKT_buffer[24] = 0.0003472222222222222*b0+0.009027777777777778*b1-0.001388888888888889*b10+0.002083333333333333*b12+0.05416666666666667*b13+0.1375*b14+0.05416666666666667*b15+0.002083333333333333*b16-0.001388888888888889*b18-0.03611111111111111*b19+0.02291666666666667*b2-0.09166666666666667*b20-0.03611111111111111*b21-0.001388888888888889*b22+0.0003472222222222222*b24+0.009027777777777778*b25+0.02291666666666667*b26+0.009027777777777778*b27+0.0003472222222222222*b28+0.009027777777777778*b3+0.0003472222222222222*b4-0.001388888888888889*b6-0.03611111111111111*b7-0.09166666666666667*b8-0.03611111111111111*b9;
        QK_B_QKT_buffer[25] = 0.1041666666666667*b15-0.01736111111111111*b1-0.006944444444444444*b10-0.01041666666666667*b12-0.1041666666666667*b13-0.001736111111111111*b0+0.01041666666666667*b16+0.006944444444444444*b18+0.06944444444444444*b19-0.06944444444444444*b21-0.006944444444444444*b22-0.001736111111111111*b24-0.01736111111111111*b25+0.01736111111111111*b27+0.001736111111111111*b28+0.01736111111111111*b3+0.001736111111111111*b4+0.006944444444444444*b6+0.06944444444444444*b7-0.06944444444444444*b9;
        QK_B_QKT_buffer[26] = 0.003472222222222222*b0+0.006944444444444444*b1-0.01388888888888889*b10+0.02083333333333333*b12+0.04166666666666667*b13-0.125*b14+0.04166666666666667*b15+0.02083333333333333*b16-0.01388888888888889*b18-0.02777777777777778*b19-0.02083333333333333*b2+0.08333333333333333*b20-0.02777777777777778*b21-0.01388888888888889*b22+0.003472222222222222*b24+0.006944444444444444*b25-0.02083333333333333*b26+0.006944444444444444*b27+0.003472222222222222*b28+0.006944444444444444*b3+0.003472222222222222*b4-0.01388888888888889*b6-0.02777777777777778*b7+0.08333333333333333*b8-0.02777777777777778*b9;
        QK_B_QKT_buffer[27] = 0.006944444444444444*b1-0.003472222222222222*b0-0.01388888888888889*b10-0.02083333333333333*b12+0.04166666666666667*b13-0.04166666666666667*b15+0.02083333333333333*b16+0.01388888888888889*b18-0.02777777777777778*b19+0.02777777777777778*b21-0.01388888888888889*b22-0.003472222222222222*b24+0.006944444444444444*b25-0.006944444444444444*b27+0.003472222222222222*b28-0.006944444444444444*b3+0.003472222222222222*b4+0.01388888888888889*b6-0.02777777777777778*b7+0.02777777777777778*b9;
        QK_B_QKT_buffer[28] = 0.001736111111111111*b0-0.006944444444444444*b1-0.006944444444444444*b10+0.01041666666666667*b12-0.04166666666666667*b13+0.0625*b14-0.04166666666666667*b15+0.01041666666666667*b16-0.006944444444444444*b18+0.02777777777777778*b19+0.01041666666666667*b2-0.04166666666666667*b20+0.02777777777777778*b21-0.006944444444444444*b22+0.001736111111111111*b24-0.006944444444444444*b25+0.01041666666666667*b26-0.006944444444444444*b27+0.001736111111111111*b28-0.006944444444444444*b3+0.001736111111111111*b4-0.006944444444444444*b6+0.02777777777777778*b7-0.04166666666666667*b8+0.02777777777777778*b9;
        QK_B_QKT_buffer[29] = 0.001736111111111111*b1-0.0003472222222222222*b0+0.006944444444444444*b10-0.001388888888888889*b11-0.002083333333333333*b12+0.01041666666666667*b13-0.02083333333333333*b14+0.02083333333333333*b15-0.01041666666666667*b16+0.002083333333333333*b17+0.001388888888888889*b18-0.006944444444444444*b19-0.003472222222222222*b2+0.01388888888888889*b20-0.01388888888888889*b21+0.006944444444444444*b22-0.001388888888888889*b23-0.0003472222222222222*b24+0.001736111111111111*b25-0.003472222222222222*b26+0.003472222222222222*b27-0.001736111111111111*b28+0.0003472222222222222*b29+0.003472222222222222*b3-0.001736111111111111*b4+0.0003472222222222222*b5+0.001388888888888889*b6-0.006944444444444444*b7+0.01388888888888889*b8-0.01388888888888889*b9;
        QK_B_QKT_buffer[30] = 0.0003472222222222222*b10-0.001805555555555556*b1-0.00006944444444444444*b0-0.0006944444444444444*b12-0.01805555555555556*b13-0.04583333333333333*b14-0.01805555555555556*b15-0.0006944444444444444*b16+0.0006944444444444444*b18+0.01805555555555556*b19-0.004583333333333333*b2+0.04583333333333333*b20+0.01805555555555556*b21+0.0006944444444444444*b22-0.0003472222222222222*b24-0.009027777777777778*b25-0.02291666666666667*b26-0.009027777777777778*b27-0.0003472222222222222*b28-0.001805555555555556*b3+0.00006944444444444444*b30+0.001805555555555556*b31+0.004583333333333333*b32+0.001805555555555556*b33+0.00006944444444444444*b34-0.00006944444444444444*b4+0.0003472222222222222*b6+0.009027777777777778*b7+0.02291666666666667*b8+0.009027777777777778*b9;
        QK_B_QKT_buffer[31] = 0.0003472222222222222*b0+0.003472222222222222*b1+0.001736111111111111*b10+0.003472222222222222*b12+0.03472222222222222*b13-0.03472222222222222*b15-0.003472222222222222*b16-0.003472222222222222*b18-0.03472222222222222*b19+0.03472222222222222*b21+0.003472222222222222*b22+0.001736111111111111*b24+0.01736111111111111*b25-0.01736111111111111*b27-0.001736111111111111*b28-0.003472222222222222*b3-0.0003472222222222222*b30-0.003472222222222222*b31+0.003472222222222222*b33+0.0003472222222222222*b34-0.0003472222222222222*b4-0.001736111111111111*b6-0.01736111111111111*b7+0.01736111111111111*b9;
        QK_B_QKT_buffer[32] = 0.003472222222222222*b10-0.001388888888888889*b1-0.0006944444444444444*b0-0.006944444444444444*b12-0.01388888888888889*b13+0.04166666666666667*b14-0.01388888888888889*b15-0.006944444444444444*b16+0.006944444444444444*b18+0.01388888888888889*b19+0.004166666666666667*b2-0.04166666666666667*b20+0.01388888888888889*b21+0.006944444444444444*b22-0.003472222222222222*b24-0.006944444444444444*b25+0.02083333333333333*b26-0.006944444444444444*b27-0.003472222222222222*b28-0.001388888888888889*b3+0.0006944444444444444*b30+0.001388888888888889*b31-0.004166666666666667*b32+0.001388888888888889*b33+0.0006944444444444444*b34-0.0006944444444444444*b4+0.003472222222222222*b6+0.006944444444444444*b7-0.02083333333333333*b8+0.006944444444444444*b9;
        QK_B_QKT_buffer[33] = 0.0006944444444444444*b0-0.001388888888888889*b1+0.003472222222222222*b10+0.006944444444444444*b12-0.01388888888888889*b13+0.01388888888888889*b15-0.006944444444444444*b16-0.006944444444444444*b18+0.01388888888888889*b19-0.01388888888888889*b21+0.006944444444444444*b22+0.003472222222222222*b24-0.006944444444444444*b25+0.006944444444444444*b27-0.003472222222222222*b28+0.001388888888888889*b3-0.0006944444444444444*b30+0.001388888888888889*b31-0.001388888888888889*b33+0.0006944444444444444*b34-0.0006944444444444444*b4-0.003472222222222222*b6+0.006944444444444444*b7-0.006944444444444444*b9;
        QK_B_QKT_buffer[34] = 0.001388888888888889*b1-0.0003472222222222222*b0+0.001736111111111111*b10-0.003472222222222222*b12+0.01388888888888889*b13-0.02083333333333333*b14+0.01388888888888889*b15-0.003472222222222222*b16+0.003472222222222222*b18-0.01388888888888889*b19-0.002083333333333333*b2+0.02083333333333333*b20-0.01388888888888889*b21+0.003472222222222222*b22-0.001736111111111111*b24+0.006944444444444444*b25-0.01041666666666667*b26+0.006944444444444444*b27-0.001736111111111111*b28+0.001388888888888889*b3+0.0003472222222222222*b30-0.001388888888888889*b31+0.002083333333333333*b32-0.001388888888888889*b33+0.0003472222222222222*b34-0.0003472222222222222*b4+0.001736111111111111*b6-0.006944444444444444*b7+0.01041666666666667*b8-0.006944444444444444*b9;
        QK_B_QKT_buffer[35] = 0.00006944444444444444*b0-0.0003472222222222222*b1-0.001736111111111111*b10+0.0003472222222222222*b11+0.0006944444444444444*b12-0.003472222222222222*b13+0.006944444444444444*b14-0.006944444444444444*b15+0.003472222222222222*b16-0.0006944444444444444*b17-0.0006944444444444444*b18+0.003472222222222222*b19+0.0006944444444444444*b2-0.006944444444444444*b20+0.006944444444444444*b21-0.003472222222222222*b22+0.0006944444444444444*b23+0.0003472222222222222*b24-0.001736111111111111*b25+0.003472222222222222*b26-0.003472222222222222*b27+0.001736111111111111*b28-0.0003472222222222222*b29-0.0006944444444444444*b3-0.00006944444444444444*b30+0.0003472222222222222*b31-0.0006944444444444444*b32+0.0006944444444444444*b33-0.0003472222222222222*b34+0.00006944444444444444*b35+0.0003472222222222222*b4-0.00006944444444444444*b5-0.0003472222222222222*b6+0.001736111111111111*b7-0.003472222222222222*b8+0.003472222222222222*b9;

        // Calculate x_star_vec for u
        x_vec_star_buffer[0] = x_vec_buffer[0]*QK_B_QKT_buffer[0]+x_vec_buffer[1]*QK_B_QKT_buffer[6]+x_vec_buffer[2]*QK_B_QKT_buffer[12]+x_vec_buffer[3]*QK_B_QKT_buffer[18]+x_vec_buffer[4]*QK_B_QKT_buffer[24]+x_vec_buffer[5]*QK_B_QKT_buffer[30];
        x_vec_star_buffer[1] = x_vec_buffer[0]*QK_B_QKT_buffer[1]+x_vec_buffer[1]*QK_B_QKT_buffer[7]+x_vec_buffer[2]*QK_B_QKT_buffer[13]+x_vec_buffer[3]*QK_B_QKT_buffer[19]+x_vec_buffer[4]*QK_B_QKT_buffer[25]+x_vec_buffer[5]*QK_B_QKT_buffer[31];
        x_vec_star_buffer[2] = x_vec_buffer[0]*QK_B_QKT_buffer[2]+x_vec_buffer[1]*QK_B_QKT_buffer[8]+x_vec_buffer[2]*QK_B_QKT_buffer[14]+x_vec_buffer[3]*QK_B_QKT_buffer[20]+x_vec_buffer[4]*QK_B_QKT_buffer[26]+x_vec_buffer[5]*QK_B_QKT_buffer[32];
        x_vec_star_buffer[3] = x_vec_buffer[0]*QK_B_QKT_buffer[3]+x_vec_buffer[1]*QK_B_QKT_buffer[9]+x_vec_buffer[2]*QK_B_QKT_buffer[15]+x_vec_buffer[3]*QK_B_QKT_buffer[21]+x_vec_buffer[4]*QK_B_QKT_buffer[27]+x_vec_buffer[5]*QK_B_QKT_buffer[33];
        x_vec_star_buffer[4] = x_vec_buffer[0]*QK_B_QKT_buffer[4]+x_vec_buffer[1]*QK_B_QKT_buffer[10]+x_vec_buffer[2]*QK_B_QKT_buffer[16]+x_vec_buffer[3]*QK_B_QKT_buffer[22]+x_vec_buffer[4]*QK_B_QKT_buffer[28]+x_vec_buffer[5]*QK_B_QKT_buffer[34];
        x_vec_star_buffer[5] = x_vec_buffer[0]*QK_B_QKT_buffer[5]+x_vec_buffer[1]*QK_B_QKT_buffer[11]+x_vec_buffer[2]*QK_B_QKT_buffer[17]+x_vec_buffer[3]*QK_B_QKT_buffer[23]+x_vec_buffer[4]*QK_B_QKT_buffer[29]+x_vec_buffer[5]*QK_B_QKT_buffer[35];

        // Calculate y_star_vec for u
        y_vec_star_buffer[0] = y_vec_buffer[0]*QK_B_QKT_buffer[0]+y_vec_buffer[1]*QK_B_QKT_buffer[1]+y_vec_buffer[2]*QK_B_QKT_buffer[2]+y_vec_buffer[3]*QK_B_QKT_buffer[3]+y_vec_buffer[4]*QK_B_QKT_buffer[4]+y_vec_buffer[5]*QK_B_QKT_buffer[5];
        y_vec_star_buffer[1] = y_vec_buffer[0]*QK_B_QKT_buffer[6]+y_vec_buffer[1]*QK_B_QKT_buffer[7]+y_vec_buffer[2]*QK_B_QKT_buffer[8]+y_vec_buffer[3]*QK_B_QKT_buffer[9]+y_vec_buffer[4]*QK_B_QKT_buffer[10]+y_vec_buffer[5]*QK_B_QKT_buffer[11];
        y_vec_star_buffer[2] = y_vec_buffer[0]*QK_B_QKT_buffer[12]+y_vec_buffer[1]*QK_B_QKT_buffer[13]+y_vec_buffer[2]*QK_B_QKT_buffer[14]+y_vec_buffer[3]*QK_B_QKT_buffer[15]+y_vec_buffer[4]*QK_B_QKT_buffer[16]+y_vec_buffer[5]*QK_B_QKT_buffer[17];
        y_vec_star_buffer[3] = y_vec_buffer[0]*QK_B_QKT_buffer[18]+y_vec_buffer[1]*QK_B_QKT_buffer[19]+y_vec_buffer[2]*QK_B_QKT_buffer[20]+y_vec_buffer[3]*QK_B_QKT_buffer[21]+y_vec_buffer[4]*QK_B_QKT_buffer[22]+y_vec_buffer[5]*QK_B_QKT_buffer[23];
        y_vec_star_buffer[4] = y_vec_buffer[0]*QK_B_QKT_buffer[24]+y_vec_buffer[1]*QK_B_QKT_buffer[25]+y_vec_buffer[2]*QK_B_QKT_buffer[26]+y_vec_buffer[3]*QK_B_QKT_buffer[27]+y_vec_buffer[4]*QK_B_QKT_buffer[28]+y_vec_buffer[5]*QK_B_QKT_buffer[29];
        y_vec_star_buffer[5] = y_vec_buffer[0]*QK_B_QKT_buffer[30]+y_vec_buffer[1]*QK_B_QKT_buffer[31]+y_vec_buffer[2]*QK_B_QKT_buffer[32]+y_vec_buffer[3]*QK_B_QKT_buffer[33]+y_vec_buffer[4]*QK_B_QKT_buffer[34]+y_vec_buffer[5]*QK_B_QKT_buffer[35];

        interpvector[0] = x_vec_star_buffer[0]*y_vec_buffer[0]+x_vec_star_buffer[1]*y_vec_buffer[1]+x_vec_star_buffer[2]*y_vec_buffer[2]+x_vec_star_buffer[3]*y_vec_buffer[3]+x_vec_star_buffer[4]*y_vec_buffer[4]+x_vec_star_buffer[5]*y_vec_buffer[5];
        interpvector[2] = x_vec_dx_buffer[0]*y_vec_star_buffer[0]+x_vec_dx_buffer[1]*y_vec_star_buffer[1]+x_vec_dx_buffer[2]*y_vec_star_buffer[2]+x_vec_dx_buffer[3]*y_vec_star_buffer[3]+x_vec_dx_buffer[4]*y_vec_star_buffer[4]+x_vec_dx_buffer[5]*y_vec_star_buffer[5];
        interpvector[3] = x_vec_star_buffer[0]*y_vec_dy_buffer[0]+x_vec_star_buffer[1]*y_vec_dy_buffer[1]+x_vec_star_buffer[2]*y_vec_dy_buffer[2]+x_vec_star_buffer[3]*y_vec_dy_buffer[3]+x_vec_star_buffer[4]*y_vec_dy_buffer[4]+x_vec_star_buffer[5]*y_vec_dy_buffer[5];

        // Compute QK*B_coef*QK^T for v-displacements next
        b0 = plot_v_interp_old[num_region_old].value[(top)+(left)*plot_v_interp_old[num_region_old].height];
        b1 = plot_v_interp_old[num_region_old].value[(top+1)+(left)*plot_v_interp_old[num_region_old].height];
        b2 = plot_v_interp_old[num_region_old].value[(top+2)+(left)*plot_v_interp_old[num_region_old].height];
        b3 = plot_v_interp_old[num_region_old].value[(top+3)+(left)*plot_v_interp_old[num_region_old].height];
        b4 = plot_v_interp_old[num_region_old].value[(top+4)+(left)*plot_v_interp_old[num_region_old].height];
        b5 = plot_v_interp_old[num_region_old].value[(top+5)+(left)*plot_v_interp_old[num_region_old].height];
        b6 = plot_v_interp_old[num_region_old].value[(top)+(left+1)*plot_v_interp_old[num_region_old].height];
        b7 = plot_v_interp_old[num_region_old].value[(top+1)+(left+1)*plot_v_interp_old[num_region_old].height];
        b8 = plot_v_interp_old[num_region_old].value[(top+2)+(left+1)*plot_v_interp_old[num_region_old].height];
        b9 = plot_v_interp_old[num_region_old].value[(top+3)+(left+1)*plot_v_interp_old[num_region_old].height];
        b10 = plot_v_interp_old[num_region_old].value[(top+4)+(left+1)*plot_v_interp_old[num_region_old].height];
        b11 = plot_v_interp_old[num_region_old].value[(top+5)+(left+1)*plot_v_interp_old[num_region_old].height];
        b12 = plot_v_interp_old[num_region_old].value[(top)+(left+2)*plot_v_interp_old[num_region_old].height];
        b13 = plot_v_interp_old[num_region_old].value[(top+1)+(left+2)*plot_v_interp_old[num_region_old].height];
        b14 = plot_v_interp_old[num_region_old].value[(top+2)+(left+2)*plot_v_interp_old[num_region_old].height];
        b15 = plot_v_interp_old[num_region_old].value[(top+3)+(left+2)*plot_v_interp_old[num_region_old].height];
        b16 = plot_v_interp_old[num_region_old].value[(top+4)+(left+2)*plot_v_interp_old[num_region_old].height];
        b17 = plot_v_interp_old[num_region_old].value[(top+5)+(left+2)*plot_v_interp_old[num_region_old].height];
        b18 = plot_v_interp_old[num_region_old].value[(top)+(left+3)*plot_v_interp_old[num_region_old].height];
        b19 = plot_v_interp_old[num_region_old].value[(top+1)+(left+3)*plot_v_interp_old[num_region_old].height];
        b20 = plot_v_interp_old[num_region_old].value[(top+2)+(left+3)*plot_v_interp_old[num_region_old].height];
        b21 = plot_v_interp_old[num_region_old].value[(top+3)+(left+3)*plot_v_interp_old[num_region_old].height];
        b22 = plot_v_interp_old[num_region_old].value[(top+4)+(left+3)*plot_v_interp_old[num_region_old].height];
        b23 = plot_v_interp_old[num_region_old].value[(top+5)+(left+3)*plot_v_interp_old[num_region_old].height];
        b24 = plot_v_interp_old[num_region_old].value[(top)+(left+4)*plot_v_interp_old[num_region_old].height];
        b25 = plot_v_interp_old[num_region_old].value[(top+1)+(left+4)*plot_v_interp_old[num_region_old].height];
        b26 = plot_v_interp_old[num_region_old].value[(top+2)+(left+4)*plot_v_interp_old[num_region_old].height];
        b27 = plot_v_interp_old[num_region_old].value[(top+3)+(left+4)*plot_v_interp_old[num_region_old].height];
        b28 = plot_v_interp_old[num_region_old].value[(top+4)+(left+4)*plot_v_interp_old[num_region_old].height];
        b29 = plot_v_interp_old[num_region_old].value[(top+5)+(left+4)*plot_v_interp_old[num_region_old].height];
        b30 = plot_v_interp_old[num_region_old].value[(top)+(left+5)*plot_v_interp_old[num_region_old].height];
        b31 = plot_v_interp_old[num_region_old].value[(top+1)+(left+5)*plot_v_interp_old[num_region_old].height];
        b32 = plot_v_interp_old[num_region_old].value[(top+2)+(left+5)*plot_v_interp_old[num_region_old].height];
        b33 = plot_v_interp_old[num_region_old].value[(top+3)+(left+5)*plot_v_interp_old[num_region_old].height];
        b34 = plot_v_interp_old[num_region_old].value[(top+4)+(left+5)*plot_v_interp_old[num_region_old].height];
        b35 = plot_v_interp_old[num_region_old].value[(top+5)+(left+5)*plot_v_interp_old[num_region_old].height];

        QK_B_QKT_buffer[0] = 0.00006944444444444444*b0+0.001805555555555556*b1+0.001805555555555556*b10+0.004583333333333333*b12+0.1191666666666667*b13+0.3025*b14+0.1191666666666667*b15+0.004583333333333333*b16+0.001805555555555556*b18+0.04694444444444444*b19+0.004583333333333333*b2+0.1191666666666667*b20+0.04694444444444444*b21+0.001805555555555556*b22+0.00006944444444444444*b24+0.001805555555555556*b25+0.004583333333333333*b26+0.001805555555555556*b27+0.00006944444444444444*b28+0.001805555555555556*b3+0.00006944444444444444*b4+0.001805555555555556*b6+0.04694444444444444*b7+0.1191666666666667*b8+0.04694444444444444*b9;
        QK_B_QKT_buffer[1] = 0.009027777777777778*b10-0.003472222222222222*b1-0.0003472222222222222*b0-0.02291666666666667*b12-0.2291666666666667*b13+0.2291666666666667*b15+0.02291666666666667*b16-0.009027777777777778*b18-0.09027777777777778*b19+0.09027777777777778*b21+0.009027777777777778*b22-0.0003472222222222222*b24-0.003472222222222222*b25+0.003472222222222222*b27+0.0003472222222222222*b28+0.003472222222222222*b3+0.0003472222222222222*b4-0.009027777777777778*b6-0.09027777777777778*b7+0.09027777777777778*b9;
        QK_B_QKT_buffer[2] = 0.0006944444444444444*b0+0.001388888888888889*b1+0.01805555555555556*b10+0.04583333333333333*b12+0.09166666666666667*b13-0.275*b14+0.09166666666666667*b15+0.04583333333333333*b16+0.01805555555555556*b18+0.03611111111111111*b19-0.004166666666666667*b2-0.1083333333333333*b20+0.03611111111111111*b21+0.01805555555555556*b22+0.0006944444444444444*b24+0.001388888888888889*b25-0.004166666666666667*b26+0.001388888888888889*b27+0.0006944444444444444*b28+0.001388888888888889*b3+0.0006944444444444444*b4+0.01805555555555556*b6+0.03611111111111111*b7-0.1083333333333333*b8+0.03611111111111111*b9;
        QK_B_QKT_buffer[3] = 0.001388888888888889*b1-0.0006944444444444444*b0+0.01805555555555556*b10-0.04583333333333333*b12+0.09166666666666667*b13-0.09166666666666667*b15+0.04583333333333333*b16-0.01805555555555556*b18+0.03611111111111111*b19-0.03611111111111111*b21+0.01805555555555556*b22-0.0006944444444444444*b24+0.001388888888888889*b25-0.001388888888888889*b27+0.0006944444444444444*b28-0.001388888888888889*b3+0.0006944444444444444*b4-0.01805555555555556*b6+0.03611111111111111*b7-0.03611111111111111*b9;
        QK_B_QKT_buffer[4] = 0.0003472222222222222*b0-0.001388888888888889*b1+0.009027777777777778*b10+0.02291666666666667*b12-0.09166666666666667*b13+0.1375*b14-0.09166666666666667*b15+0.02291666666666667*b16+0.009027777777777778*b18-0.03611111111111111*b19+0.002083333333333333*b2+0.05416666666666667*b20-0.03611111111111111*b21+0.009027777777777778*b22+0.0003472222222222222*b24-0.001388888888888889*b25+0.002083333333333333*b26-0.001388888888888889*b27+0.0003472222222222222*b28-0.001388888888888889*b3+0.0003472222222222222*b4+0.009027777777777778*b6-0.03611111111111111*b7+0.05416666666666667*b8-0.03611111111111111*b9;
        QK_B_QKT_buffer[5] = 0.0003472222222222222*b1-0.00006944444444444444*b0-0.009027777777777778*b10+0.001805555555555556*b11-0.004583333333333333*b12+0.02291666666666667*b13-0.04583333333333333*b14+0.04583333333333333*b15-0.02291666666666667*b16+0.004583333333333333*b17-0.001805555555555556*b18+0.009027777777777778*b19-0.0006944444444444444*b2-0.01805555555555556*b20+0.01805555555555556*b21-0.009027777777777778*b22+0.001805555555555556*b23-0.00006944444444444444*b24+0.0003472222222222222*b25-0.0006944444444444444*b26+0.0006944444444444444*b27-0.0003472222222222222*b28+0.00006944444444444444*b29+0.0006944444444444444*b3-0.0003472222222222222*b4+0.00006944444444444444*b5-0.001805555555555556*b6+0.009027777777777778*b7-0.01805555555555556*b8+0.01805555555555556*b9;
        QK_B_QKT_buffer[6] = 0.003472222222222222*b18-0.009027777777777778*b1-0.003472222222222222*b10-0.0003472222222222222*b0+0.09027777777777778*b19-0.02291666666666667*b2+0.2291666666666667*b20+0.09027777777777778*b21+0.003472222222222222*b22+0.0003472222222222222*b24+0.009027777777777778*b25+0.02291666666666667*b26+0.009027777777777778*b27+0.0003472222222222222*b28-0.009027777777777778*b3-0.0003472222222222222*b4-0.003472222222222222*b6-0.09027777777777778*b7-0.2291666666666667*b8-0.09027777777777778*b9;
        QK_B_QKT_buffer[7] = 0.001736111111111111*b0+0.01736111111111111*b1-0.01736111111111111*b10-0.01736111111111111*b18-0.1736111111111111*b19+0.1736111111111111*b21+0.01736111111111111*b22-0.001736111111111111*b24-0.01736111111111111*b25+0.01736111111111111*b27+0.001736111111111111*b28-0.01736111111111111*b3-0.001736111111111111*b4+0.01736111111111111*b6+0.1736111111111111*b7-0.1736111111111111*b9;
        QK_B_QKT_buffer[8] = 0.03472222222222222*b18-0.006944444444444444*b1-0.03472222222222222*b10-0.003472222222222222*b0+0.06944444444444444*b19+0.02083333333333333*b2-0.2083333333333333*b20+0.06944444444444444*b21+0.03472222222222222*b22+0.003472222222222222*b24+0.006944444444444444*b25-0.02083333333333333*b26+0.006944444444444444*b27+0.003472222222222222*b28-0.006944444444444444*b3-0.003472222222222222*b4-0.03472222222222222*b6-0.06944444444444444*b7+0.2083333333333333*b8-0.06944444444444444*b9;
        QK_B_QKT_buffer[9] = 0.003472222222222222*b0-0.006944444444444444*b1-0.03472222222222222*b10-0.03472222222222222*b18+0.06944444444444444*b19-0.06944444444444444*b21+0.03472222222222222*b22-0.003472222222222222*b24+0.006944444444444444*b25-0.006944444444444444*b27+0.003472222222222222*b28+0.006944444444444444*b3-0.003472222222222222*b4+0.03472222222222222*b6-0.06944444444444444*b7+0.06944444444444444*b9;
        QK_B_QKT_buffer[10] = 0.006944444444444444*b1-0.001736111111111111*b0-0.01736111111111111*b10+0.01736111111111111*b18-0.06944444444444444*b19-0.01041666666666667*b2+0.1041666666666667*b20-0.06944444444444444*b21+0.01736111111111111*b22+0.001736111111111111*b24-0.006944444444444444*b25+0.01041666666666667*b26-0.006944444444444444*b27+0.001736111111111111*b28+0.006944444444444444*b3-0.001736111111111111*b4-0.01736111111111111*b6+0.06944444444444444*b7-0.1041666666666667*b8+0.06944444444444444*b9;
        QK_B_QKT_buffer[11] = 0.0003472222222222222*b0-0.001736111111111111*b1+0.01736111111111111*b10-0.003472222222222222*b11-0.003472222222222222*b18+0.01736111111111111*b19+0.003472222222222222*b2-0.03472222222222222*b20+0.03472222222222222*b21-0.01736111111111111*b22+0.003472222222222222*b23-0.0003472222222222222*b24+0.001736111111111111*b25-0.003472222222222222*b26+0.003472222222222222*b27-0.001736111111111111*b28+0.0003472222222222222*b29-0.003472222222222222*b3+0.001736111111111111*b4-0.0003472222222222222*b5+0.003472222222222222*b6-0.01736111111111111*b7+0.03472222222222222*b8-0.03472222222222222*b9;
        QK_B_QKT_buffer[12] = 0.0006944444444444444*b0+0.01805555555555556*b1+0.001388888888888889*b10-0.004166666666666667*b12-0.1083333333333333*b13-0.275*b14-0.1083333333333333*b15-0.004166666666666667*b16+0.001388888888888889*b18+0.03611111111111111*b19+0.04583333333333333*b2+0.09166666666666667*b20+0.03611111111111111*b21+0.001388888888888889*b22+0.0006944444444444444*b24+0.01805555555555556*b25+0.04583333333333333*b26+0.01805555555555556*b27+0.0006944444444444444*b28+0.01805555555555556*b3+0.0006944444444444444*b4+0.001388888888888889*b6+0.03611111111111111*b7+0.09166666666666667*b8+0.03611111111111111*b9;
        QK_B_QKT_buffer[13] = 0.006944444444444444*b10-0.03472222222222222*b1-0.003472222222222222*b0+0.02083333333333333*b12+0.2083333333333333*b13-0.2083333333333333*b15-0.02083333333333333*b16-0.006944444444444444*b18-0.06944444444444444*b19+0.06944444444444444*b21+0.006944444444444444*b22-0.003472222222222222*b24-0.03472222222222222*b25+0.03472222222222222*b27+0.003472222222222222*b28+0.03472222222222222*b3+0.003472222222222222*b4-0.006944444444444444*b6-0.06944444444444444*b7+0.06944444444444444*b9;
        QK_B_QKT_buffer[14] = 0.006944444444444444*b0+0.01388888888888889*b1+0.01388888888888889*b10-0.04166666666666667*b12-0.08333333333333333*b13+0.25*b14-0.08333333333333333*b15-0.04166666666666667*b16+0.01388888888888889*b18+0.02777777777777778*b19-0.04166666666666667*b2-0.08333333333333333*b20+0.02777777777777778*b21+0.01388888888888889*b22+0.006944444444444444*b24+0.01388888888888889*b25-0.04166666666666667*b26+0.01388888888888889*b27+0.006944444444444444*b28+0.01388888888888889*b3+0.006944444444444444*b4+0.01388888888888889*b6+0.02777777777777778*b7-0.08333333333333333*b8+0.02777777777777778*b9;
        QK_B_QKT_buffer[15] = 0.01388888888888889*b1-0.006944444444444444*b0+0.01388888888888889*b10+0.04166666666666667*b12-0.08333333333333333*b13+0.08333333333333333*b15-0.04166666666666667*b16-0.01388888888888889*b18+0.02777777777777778*b19-0.02777777777777778*b21+0.01388888888888889*b22-0.006944444444444444*b24+0.01388888888888889*b25-0.01388888888888889*b27+0.006944444444444444*b28-0.01388888888888889*b3+0.006944444444444444*b4-0.01388888888888889*b6+0.02777777777777778*b7-0.02777777777777778*b9;
        QK_B_QKT_buffer[16] = 0.003472222222222222*b0-0.01388888888888889*b1+0.006944444444444444*b10-0.02083333333333333*b12+0.08333333333333333*b13-0.125*b14+0.08333333333333333*b15-0.02083333333333333*b16+0.006944444444444444*b18-0.02777777777777778*b19+0.02083333333333333*b2+0.04166666666666667*b20-0.02777777777777778*b21+0.006944444444444444*b22+0.003472222222222222*b24-0.01388888888888889*b25+0.02083333333333333*b26-0.01388888888888889*b27+0.003472222222222222*b28-0.01388888888888889*b3+0.003472222222222222*b4+0.006944444444444444*b6-0.02777777777777778*b7+0.04166666666666667*b8-0.02777777777777778*b9;
        QK_B_QKT_buffer[17] = 0.003472222222222222*b1-0.0006944444444444444*b0-0.006944444444444444*b10+0.001388888888888889*b11+0.004166666666666667*b12-0.02083333333333333*b13+0.04166666666666667*b14-0.04166666666666667*b15+0.02083333333333333*b16-0.004166666666666667*b17-0.001388888888888889*b18+0.006944444444444444*b19-0.006944444444444444*b2-0.01388888888888889*b20+0.01388888888888889*b21-0.006944444444444444*b22+0.001388888888888889*b23-0.0006944444444444444*b24+0.003472222222222222*b25-0.006944444444444444*b26+0.006944444444444444*b27-0.003472222222222222*b28+0.0006944444444444444*b29+0.006944444444444444*b3-0.003472222222222222*b4+0.0006944444444444444*b5-0.001388888888888889*b6+0.006944444444444444*b7-0.01388888888888889*b8+0.01388888888888889*b9;
        QK_B_QKT_buffer[18] = 0.001388888888888889*b10-0.01805555555555556*b1-0.0006944444444444444*b0-0.001388888888888889*b18-0.03611111111111111*b19-0.04583333333333333*b2-0.09166666666666667*b20-0.03611111111111111*b21-0.001388888888888889*b22+0.0006944444444444444*b24+0.01805555555555556*b25+0.04583333333333333*b26+0.01805555555555556*b27+0.0006944444444444444*b28-0.01805555555555556*b3-0.0006944444444444444*b4+0.001388888888888889*b6+0.03611111111111111*b7+0.09166666666666667*b8+0.03611111111111111*b9;
        QK_B_QKT_buffer[19] = 0.003472222222222222*b0+0.03472222222222222*b1+0.006944444444444444*b10+0.006944444444444444*b18+0.06944444444444444*b19-0.06944444444444444*b21-0.006944444444444444*b22-0.003472222222222222*b24-0.03472222222222222*b25+0.03472222222222222*b27+0.003472222222222222*b28-0.03472222222222222*b3-0.003472222222222222*b4-0.006944444444444444*b6-0.06944444444444444*b7+0.06944444444444444*b9;
        QK_B_QKT_buffer[20] = 0.01388888888888889*b10-0.01388888888888889*b1-0.006944444444444444*b0-0.01388888888888889*b18-0.02777777777777778*b19+0.04166666666666667*b2+0.08333333333333333*b20-0.02777777777777778*b21-0.01388888888888889*b22+0.006944444444444444*b24+0.01388888888888889*b25-0.04166666666666667*b26+0.01388888888888889*b27+0.006944444444444444*b28-0.01388888888888889*b3-0.006944444444444444*b4+0.01388888888888889*b6+0.02777777777777778*b7-0.08333333333333333*b8+0.02777777777777778*b9;
        QK_B_QKT_buffer[21] = 0.006944444444444444*b0-0.01388888888888889*b1+0.01388888888888889*b10+0.01388888888888889*b18-0.02777777777777778*b19+0.02777777777777778*b21-0.01388888888888889*b22-0.006944444444444444*b24+0.01388888888888889*b25-0.01388888888888889*b27+0.006944444444444444*b28+0.01388888888888889*b3-0.006944444444444444*b4-0.01388888888888889*b6+0.02777777777777778*b7-0.02777777777777778*b9;
        QK_B_QKT_buffer[22] = 0.01388888888888889*b1-0.003472222222222222*b0+0.006944444444444444*b10-0.006944444444444444*b18+0.02777777777777778*b19-0.02083333333333333*b2-0.04166666666666667*b20+0.02777777777777778*b21-0.006944444444444444*b22+0.003472222222222222*b24-0.01388888888888889*b25+0.02083333333333333*b26-0.01388888888888889*b27+0.003472222222222222*b28+0.01388888888888889*b3-0.003472222222222222*b4+0.006944444444444444*b6-0.02777777777777778*b7+0.04166666666666667*b8-0.02777777777777778*b9;
        QK_B_QKT_buffer[23] = 0.0006944444444444444*b0-0.003472222222222222*b1-0.006944444444444444*b10+0.001388888888888889*b11+0.001388888888888889*b18-0.006944444444444444*b19+0.006944444444444444*b2+0.01388888888888889*b20-0.01388888888888889*b21+0.006944444444444444*b22-0.001388888888888889*b23-0.0006944444444444444*b24+0.003472222222222222*b25-0.006944444444444444*b26+0.006944444444444444*b27-0.003472222222222222*b28+0.0006944444444444444*b29-0.006944444444444444*b3+0.003472222222222222*b4-0.0006944444444444444*b5-0.001388888888888889*b6+0.006944444444444444*b7-0.01388888888888889*b8+0.01388888888888889*b9;
        QK_B_QKT_buffer[24] = 0.0003472222222222222*b0+0.009027777777777778*b1-0.001388888888888889*b10+0.002083333333333333*b12+0.05416666666666667*b13+0.1375*b14+0.05416666666666667*b15+0.002083333333333333*b16-0.001388888888888889*b18-0.03611111111111111*b19+0.02291666666666667*b2-0.09166666666666667*b20-0.03611111111111111*b21-0.001388888888888889*b22+0.0003472222222222222*b24+0.009027777777777778*b25+0.02291666666666667*b26+0.009027777777777778*b27+0.0003472222222222222*b28+0.009027777777777778*b3+0.0003472222222222222*b4-0.001388888888888889*b6-0.03611111111111111*b7-0.09166666666666667*b8-0.03611111111111111*b9;
        QK_B_QKT_buffer[25] = 0.1041666666666667*b15-0.01736111111111111*b1-0.006944444444444444*b10-0.01041666666666667*b12-0.1041666666666667*b13-0.001736111111111111*b0+0.01041666666666667*b16+0.006944444444444444*b18+0.06944444444444444*b19-0.06944444444444444*b21-0.006944444444444444*b22-0.001736111111111111*b24-0.01736111111111111*b25+0.01736111111111111*b27+0.001736111111111111*b28+0.01736111111111111*b3+0.001736111111111111*b4+0.006944444444444444*b6+0.06944444444444444*b7-0.06944444444444444*b9;
        QK_B_QKT_buffer[26] = 0.003472222222222222*b0+0.006944444444444444*b1-0.01388888888888889*b10+0.02083333333333333*b12+0.04166666666666667*b13-0.125*b14+0.04166666666666667*b15+0.02083333333333333*b16-0.01388888888888889*b18-0.02777777777777778*b19-0.02083333333333333*b2+0.08333333333333333*b20-0.02777777777777778*b21-0.01388888888888889*b22+0.003472222222222222*b24+0.006944444444444444*b25-0.02083333333333333*b26+0.006944444444444444*b27+0.003472222222222222*b28+0.006944444444444444*b3+0.003472222222222222*b4-0.01388888888888889*b6-0.02777777777777778*b7+0.08333333333333333*b8-0.02777777777777778*b9;
        QK_B_QKT_buffer[27] = 0.006944444444444444*b1-0.003472222222222222*b0-0.01388888888888889*b10-0.02083333333333333*b12+0.04166666666666667*b13-0.04166666666666667*b15+0.02083333333333333*b16+0.01388888888888889*b18-0.02777777777777778*b19+0.02777777777777778*b21-0.01388888888888889*b22-0.003472222222222222*b24+0.006944444444444444*b25-0.006944444444444444*b27+0.003472222222222222*b28-0.006944444444444444*b3+0.003472222222222222*b4+0.01388888888888889*b6-0.02777777777777778*b7+0.02777777777777778*b9;
        QK_B_QKT_buffer[28] = 0.001736111111111111*b0-0.006944444444444444*b1-0.006944444444444444*b10+0.01041666666666667*b12-0.04166666666666667*b13+0.0625*b14-0.04166666666666667*b15+0.01041666666666667*b16-0.006944444444444444*b18+0.02777777777777778*b19+0.01041666666666667*b2-0.04166666666666667*b20+0.02777777777777778*b21-0.006944444444444444*b22+0.001736111111111111*b24-0.006944444444444444*b25+0.01041666666666667*b26-0.006944444444444444*b27+0.001736111111111111*b28-0.006944444444444444*b3+0.001736111111111111*b4-0.006944444444444444*b6+0.02777777777777778*b7-0.04166666666666667*b8+0.02777777777777778*b9;
        QK_B_QKT_buffer[29] = 0.001736111111111111*b1-0.0003472222222222222*b0+0.006944444444444444*b10-0.001388888888888889*b11-0.002083333333333333*b12+0.01041666666666667*b13-0.02083333333333333*b14+0.02083333333333333*b15-0.01041666666666667*b16+0.002083333333333333*b17+0.001388888888888889*b18-0.006944444444444444*b19-0.003472222222222222*b2+0.01388888888888889*b20-0.01388888888888889*b21+0.006944444444444444*b22-0.001388888888888889*b23-0.0003472222222222222*b24+0.001736111111111111*b25-0.003472222222222222*b26+0.003472222222222222*b27-0.001736111111111111*b28+0.0003472222222222222*b29+0.003472222222222222*b3-0.001736111111111111*b4+0.0003472222222222222*b5+0.001388888888888889*b6-0.006944444444444444*b7+0.01388888888888889*b8-0.01388888888888889*b9;
        QK_B_QKT_buffer[30] = 0.0003472222222222222*b10-0.001805555555555556*b1-0.00006944444444444444*b0-0.0006944444444444444*b12-0.01805555555555556*b13-0.04583333333333333*b14-0.01805555555555556*b15-0.0006944444444444444*b16+0.0006944444444444444*b18+0.01805555555555556*b19-0.004583333333333333*b2+0.04583333333333333*b20+0.01805555555555556*b21+0.0006944444444444444*b22-0.0003472222222222222*b24-0.009027777777777778*b25-0.02291666666666667*b26-0.009027777777777778*b27-0.0003472222222222222*b28-0.001805555555555556*b3+0.00006944444444444444*b30+0.001805555555555556*b31+0.004583333333333333*b32+0.001805555555555556*b33+0.00006944444444444444*b34-0.00006944444444444444*b4+0.0003472222222222222*b6+0.009027777777777778*b7+0.02291666666666667*b8+0.009027777777777778*b9;
        QK_B_QKT_buffer[31] = 0.0003472222222222222*b0+0.003472222222222222*b1+0.001736111111111111*b10+0.003472222222222222*b12+0.03472222222222222*b13-0.03472222222222222*b15-0.003472222222222222*b16-0.003472222222222222*b18-0.03472222222222222*b19+0.03472222222222222*b21+0.003472222222222222*b22+0.001736111111111111*b24+0.01736111111111111*b25-0.01736111111111111*b27-0.001736111111111111*b28-0.003472222222222222*b3-0.0003472222222222222*b30-0.003472222222222222*b31+0.003472222222222222*b33+0.0003472222222222222*b34-0.0003472222222222222*b4-0.001736111111111111*b6-0.01736111111111111*b7+0.01736111111111111*b9;
        QK_B_QKT_buffer[32] = 0.003472222222222222*b10-0.001388888888888889*b1-0.0006944444444444444*b0-0.006944444444444444*b12-0.01388888888888889*b13+0.04166666666666667*b14-0.01388888888888889*b15-0.006944444444444444*b16+0.006944444444444444*b18+0.01388888888888889*b19+0.004166666666666667*b2-0.04166666666666667*b20+0.01388888888888889*b21+0.006944444444444444*b22-0.003472222222222222*b24-0.006944444444444444*b25+0.02083333333333333*b26-0.006944444444444444*b27-0.003472222222222222*b28-0.001388888888888889*b3+0.0006944444444444444*b30+0.001388888888888889*b31-0.004166666666666667*b32+0.001388888888888889*b33+0.0006944444444444444*b34-0.0006944444444444444*b4+0.003472222222222222*b6+0.006944444444444444*b7-0.02083333333333333*b8+0.006944444444444444*b9;
        QK_B_QKT_buffer[33] = 0.0006944444444444444*b0-0.001388888888888889*b1+0.003472222222222222*b10+0.006944444444444444*b12-0.01388888888888889*b13+0.01388888888888889*b15-0.006944444444444444*b16-0.006944444444444444*b18+0.01388888888888889*b19-0.01388888888888889*b21+0.006944444444444444*b22+0.003472222222222222*b24-0.006944444444444444*b25+0.006944444444444444*b27-0.003472222222222222*b28+0.001388888888888889*b3-0.0006944444444444444*b30+0.001388888888888889*b31-0.001388888888888889*b33+0.0006944444444444444*b34-0.0006944444444444444*b4-0.003472222222222222*b6+0.006944444444444444*b7-0.006944444444444444*b9;
        QK_B_QKT_buffer[34] = 0.001388888888888889*b1-0.0003472222222222222*b0+0.001736111111111111*b10-0.003472222222222222*b12+0.01388888888888889*b13-0.02083333333333333*b14+0.01388888888888889*b15-0.003472222222222222*b16+0.003472222222222222*b18-0.01388888888888889*b19-0.002083333333333333*b2+0.02083333333333333*b20-0.01388888888888889*b21+0.003472222222222222*b22-0.001736111111111111*b24+0.006944444444444444*b25-0.01041666666666667*b26+0.006944444444444444*b27-0.001736111111111111*b28+0.001388888888888889*b3+0.0003472222222222222*b30-0.001388888888888889*b31+0.002083333333333333*b32-0.001388888888888889*b33+0.0003472222222222222*b34-0.0003472222222222222*b4+0.001736111111111111*b6-0.006944444444444444*b7+0.01041666666666667*b8-0.006944444444444444*b9;
        QK_B_QKT_buffer[35] = 0.00006944444444444444*b0-0.0003472222222222222*b1-0.001736111111111111*b10+0.0003472222222222222*b11+0.0006944444444444444*b12-0.003472222222222222*b13+0.006944444444444444*b14-0.006944444444444444*b15+0.003472222222222222*b16-0.0006944444444444444*b17-0.0006944444444444444*b18+0.003472222222222222*b19+0.0006944444444444444*b2-0.006944444444444444*b20+0.006944444444444444*b21-0.003472222222222222*b22+0.0006944444444444444*b23+0.0003472222222222222*b24-0.001736111111111111*b25+0.003472222222222222*b26-0.003472222222222222*b27+0.001736111111111111*b28-0.0003472222222222222*b29-0.0006944444444444444*b3-0.00006944444444444444*b30+0.0003472222222222222*b31-0.0006944444444444444*b32+0.0006944444444444444*b33-0.0003472222222222222*b34+0.00006944444444444444*b35+0.0003472222222222222*b4-0.00006944444444444444*b5-0.0003472222222222222*b6+0.001736111111111111*b7-0.003472222222222222*b8+0.003472222222222222*b9;

        // Calculate x_star_vec for v
        x_vec_star_buffer[0] = x_vec_buffer[0]*QK_B_QKT_buffer[0]+x_vec_buffer[1]*QK_B_QKT_buffer[6]+x_vec_buffer[2]*QK_B_QKT_buffer[12]+x_vec_buffer[3]*QK_B_QKT_buffer[18]+x_vec_buffer[4]*QK_B_QKT_buffer[24]+x_vec_buffer[5]*QK_B_QKT_buffer[30];
        x_vec_star_buffer[1] = x_vec_buffer[0]*QK_B_QKT_buffer[1]+x_vec_buffer[1]*QK_B_QKT_buffer[7]+x_vec_buffer[2]*QK_B_QKT_buffer[13]+x_vec_buffer[3]*QK_B_QKT_buffer[19]+x_vec_buffer[4]*QK_B_QKT_buffer[25]+x_vec_buffer[5]*QK_B_QKT_buffer[31];
        x_vec_star_buffer[2] = x_vec_buffer[0]*QK_B_QKT_buffer[2]+x_vec_buffer[1]*QK_B_QKT_buffer[8]+x_vec_buffer[2]*QK_B_QKT_buffer[14]+x_vec_buffer[3]*QK_B_QKT_buffer[20]+x_vec_buffer[4]*QK_B_QKT_buffer[26]+x_vec_buffer[5]*QK_B_QKT_buffer[32];
        x_vec_star_buffer[3] = x_vec_buffer[0]*QK_B_QKT_buffer[3]+x_vec_buffer[1]*QK_B_QKT_buffer[9]+x_vec_buffer[2]*QK_B_QKT_buffer[15]+x_vec_buffer[3]*QK_B_QKT_buffer[21]+x_vec_buffer[4]*QK_B_QKT_buffer[27]+x_vec_buffer[5]*QK_B_QKT_buffer[33];
        x_vec_star_buffer[4] = x_vec_buffer[0]*QK_B_QKT_buffer[4]+x_vec_buffer[1]*QK_B_QKT_buffer[10]+x_vec_buffer[2]*QK_B_QKT_buffer[16]+x_vec_buffer[3]*QK_B_QKT_buffer[22]+x_vec_buffer[4]*QK_B_QKT_buffer[28]+x_vec_buffer[5]*QK_B_QKT_buffer[34];
        x_vec_star_buffer[5] = x_vec_buffer[0]*QK_B_QKT_buffer[5]+x_vec_buffer[1]*QK_B_QKT_buffer[11]+x_vec_buffer[2]*QK_B_QKT_buffer[17]+x_vec_buffer[3]*QK_B_QKT_buffer[23]+x_vec_buffer[4]*QK_B_QKT_buffer[29]+x_vec_buffer[5]*QK_B_QKT_buffer[35];

        // Calculate y_star_vec for v
        y_vec_star_buffer[0] = y_vec_buffer[0]*QK_B_QKT_buffer[0]+y_vec_buffer[1]*QK_B_QKT_buffer[1]+y_vec_buffer[2]*QK_B_QKT_buffer[2]+y_vec_buffer[3]*QK_B_QKT_buffer[3]+y_vec_buffer[4]*QK_B_QKT_buffer[4]+y_vec_buffer[5]*QK_B_QKT_buffer[5];
        y_vec_star_buffer[1] = y_vec_buffer[0]*QK_B_QKT_buffer[6]+y_vec_buffer[1]*QK_B_QKT_buffer[7]+y_vec_buffer[2]*QK_B_QKT_buffer[8]+y_vec_buffer[3]*QK_B_QKT_buffer[9]+y_vec_buffer[4]*QK_B_QKT_buffer[10]+y_vec_buffer[5]*QK_B_QKT_buffer[11];
        y_vec_star_buffer[2] = y_vec_buffer[0]*QK_B_QKT_buffer[12]+y_vec_buffer[1]*QK_B_QKT_buffer[13]+y_vec_buffer[2]*QK_B_QKT_buffer[14]+y_vec_buffer[3]*QK_B_QKT_buffer[15]+y_vec_buffer[4]*QK_B_QKT_buffer[16]+y_vec_buffer[5]*QK_B_QKT_buffer[17];
        y_vec_star_buffer[3] = y_vec_buffer[0]*QK_B_QKT_buffer[18]+y_vec_buffer[1]*QK_B_QKT_buffer[19]+y_vec_buffer[2]*QK_B_QKT_buffer[20]+y_vec_buffer[3]*QK_B_QKT_buffer[21]+y_vec_buffer[4]*QK_B_QKT_buffer[22]+y_vec_buffer[5]*QK_B_QKT_buffer[23];
        y_vec_star_buffer[4] = y_vec_buffer[0]*QK_B_QKT_buffer[24]+y_vec_buffer[1]*QK_B_QKT_buffer[25]+y_vec_buffer[2]*QK_B_QKT_buffer[26]+y_vec_buffer[3]*QK_B_QKT_buffer[27]+y_vec_buffer[4]*QK_B_QKT_buffer[28]+y_vec_buffer[5]*QK_B_QKT_buffer[29];
        y_vec_star_buffer[5] = y_vec_buffer[0]*QK_B_QKT_buffer[30]+y_vec_buffer[1]*QK_B_QKT_buffer[31]+y_vec_buffer[2]*QK_B_QKT_buffer[32]+y_vec_buffer[3]*QK_B_QKT_buffer[33]+y_vec_buffer[4]*QK_B_QKT_buffer[34]+y_vec_buffer[5]*QK_B_QKT_buffer[35];

        interpvector[1] = x_vec_star_buffer[0]*y_vec_buffer[0]+x_vec_star_buffer[1]*y_vec_buffer[1]+x_vec_star_buffer[2]*y_vec_buffer[2]+x_vec_star_buffer[3]*y_vec_buffer[3]+x_vec_star_buffer[4]*y_vec_buffer[4]+x_vec_star_buffer[5]*y_vec_buffer[5];
        interpvector[4] = x_vec_dx_buffer[0]*y_vec_star_buffer[0]+x_vec_dx_buffer[1]*y_vec_star_buffer[1]+x_vec_dx_buffer[2]*y_vec_star_buffer[2]+x_vec_dx_buffer[3]*y_vec_star_buffer[3]+x_vec_dx_buffer[4]*y_vec_star_buffer[4]+x_vec_dx_buffer[5]*y_vec_star_buffer[5];
        interpvector[5] = x_vec_star_buffer[0]*y_vec_dy_buffer[0]+x_vec_star_buffer[1]*y_vec_dy_buffer[1]+x_vec_star_buffer[2]*y_vec_dy_buffer[2]+x_vec_star_buffer[3]*y_vec_dy_buffer[3]+x_vec_star_buffer[4]*y_vec_dy_buffer[4]+x_vec_star_buffer[5]*y_vec_dy_buffer[5];
        return SUCCESS;
    }    
    return FAILED;
}

void mexFunction(int nlhs,mxArray *plhs[ ],int nrhs,const mxArray *prhs[ ]) {
    if (nrhs == 9 && nlhs == 2) {
        // Create convert 
        class_convert convert(plhs,prhs);
        
        // Run analysis
        convert.analysis();
    } else {
        mexErrMsgTxt("Incorrect number of inputs.\n");
    }
}
