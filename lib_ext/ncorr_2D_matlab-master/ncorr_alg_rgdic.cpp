// This function calculates the displacement fields using RG-DIC. If no flags are passed, then a simple single threaded version will be compiled
// Since this is multithreaded, make sure all functions called in the OpenMP region are threadsafe!

#include <mex.h>
#include <math.h>
#include <queue>
#include <vector>
#include <list>
#include <exception>              // Allow exceptions because this function can allocate a very large chunk of memory for the interpolation lookup table
#include <complex>
#include "standard_datatypes.h"
#include "ncorr_datatypes.h"
#include "ncorr_lib.h"
            
#ifdef NCORR_OPENMP
    #include <omp.h>              // openmp header
#endif

// ----------------------------------------------------------------------//
// Seed info ------------------------------------------------------------//
// ----------------------------------------------------------------------//

struct local_struct_seedinfo {
    // Constructor
    local_struct_seedinfo()   {
        num_region = 0;
        num_thread = 0;
        computepoints = 0;
    }
    
    // Properties
    class_double_array paramvector;
    int num_region;
    int num_thread;
    int computepoints;
};

void get_seedinfo(std::vector<local_struct_seedinfo> &seedinfo,const mxArray *prhs) {    
    // Check input    
    if (mxIsClass(prhs,"struct")) {
        for (int i=0; i<(int)mxGetN(prhs); i++) {
            for (int j=0; j<(int)mxGetM(prhs); j++) {
                int total_threads = (int)mxGetM(prhs);
                mxArray *mat_paramvector = mxGetField(prhs,j+i*total_threads,"paramvector");
                mxArray *mat_num_region = mxGetField(prhs,j+i*total_threads,"num_region");        
                mxArray *mat_num_thread = mxGetField(prhs,j+i*total_threads,"num_thread");    
                mxArray *mat_computepoints = mxGetField(prhs,j+i*total_threads,"computepoints");    
                if (mat_paramvector != 0 && mat_num_region != 0 && mat_num_thread != 0 && mat_computepoints != 0) {
                    // Form seedinfo
                    local_struct_seedinfo seedinfo_template;
                    get_double_array(seedinfo_template.paramvector,mat_paramvector);
                    get_integer_scalar(seedinfo_template.num_region,mat_num_region);    
                    get_integer_scalar(seedinfo_template.num_thread,mat_num_thread);    
                    get_integer_scalar(seedinfo_template.computepoints,mat_computepoints);
                    
                    // Test sizes
                    if (seedinfo_template.paramvector.width != 9 || seedinfo_template.paramvector.height != 1) {
                        mexErrMsgTxt("'paramvector' is supposed to be a 1x9 vector.\n");
                    }
                    
                    // Store seedinfo_template
                    seedinfo.push_back(seedinfo_template);
                } else {
                    mexErrMsgTxt("Some fields are missing for seedinfo.\n");
                }
            }
        }
    } else {
        mexErrMsgTxt("Seedinfo must be of class 'struct'.\n");
    }
}

// ----------------------------------------------------------------------//
// Queue ----------------------------------------------------------------//
// ----------------------------------------------------------------------//

struct comp_queue {
    bool operator ()(std::vector<double> const& a, std::vector<double> const& b) const {
        // Eighth element is the correlation coefficient
        return a[8] > b[8];
    }
};

typedef std::priority_queue<std::vector<double>, std::vector<std::vector<double> >,  comp_queue> heap;


// ----------------------------------------------------//
// Main Class -----------------------------------------//
// ----------------------------------------------------//

class class_rgdic {
public:
    // Constructor
    class_rgdic(mxArray *plhs [ ],const mxArray *prhs [ ]);
    
    // Methods
    void analysis();

private:
    // Properties
    // Inputs:
    std::vector<ncorr_class_img> reference;       // ncorr datatype
    std::vector<ncorr_class_img> current;         // ncorr datatype
    std::vector<ncorr_class_roi> roi;             // ncorr datatype
    std::vector<local_struct_seedinfo> seedinfo;  // local datatype
    class_integer_array threaddiagram;            // standard datatype
    int radius;                                   // standard datatype
    int spacing;                                  // standard datatype
    double cutoff_diffnorm;                       // standard datatype
    int cutoff_iteration;                         // standard datatype
    bool subsettrunc;                             // standard datatype
    int num_img;                                  // standard datatype
    int total_imgs;                               // standard datatype
          
    // Outputs:
    class_double_array plot_u;
    class_double_array plot_v;
    class_double_array plot_corrcoef;
    class_logical_array plot_validpoints;
    double *outstate;
    
    // Other variables:
    int total_threads;
    int total_region;
    class_waitbar waitbar;
    int num_thread_waitbar;               // Keeps track of which thread to base update on
    class_logical_array plot_calcpoints;  // Keeps track of points which have already been calculated    

    // Inverse Compositional Buffers: Must form these explicitly
    // Thread independent: These are thread independent because they are read only
    std::vector<double> df_dx_buffer;     // This is the reference gradient
    std::vector<double> df_dy_buffer;     // This is the reference gradient    
    std::vector<double> QK_B_QKT_buffer;  // Very large 36 x size of reference image
    // Thread dependent: Must form buffers explicitly for each thread to make them thread safe, because these are read/write
    std::vector<std::vector<double> > g_buffer; 
    std::vector<std::vector<double> > df_dp_buffer;  // First order "steepest descent images" - these are computed once per iteration   
    std::vector<std::vector<double> > x_vec_buffer;
    std::vector<std::vector<double> > y_vec_buffer; 
    std::vector<std::vector<double> > gradient_buffer;
    std::vector<std::vector<double> > hessian_gn_buffer;  
    
    // Methods
    void precompute();
    void analyzepoint(heap &queue,const int &x,const int &y,const std::vector<double> &paramvector_init,const int &num_region,const int &num_thread);
    OUT calcpoint(std::vector<double> &paramvector,const int &x,const int &y,const std::vector<double> &paramvector_init,const int &num_region,const int &num_thread);
    OUT iterativesearch(std::vector<double> &defvector,double &corrcoef,const std::vector<double> &defvector_init,const int &num_thread);
    OUT newton(std::vector<double> &defvector,double &corrcoef,double &diffnorm,const std::vector<double> &defvector_init,const double &fm,const double &deltaf_inv,const int &num_thread);
};

class_rgdic::class_rgdic(mxArray *plhs[ ],const mxArray *prhs[ ]) {
    // Get inputs --------------------------------------------------// 
    // input 1: reference image    
    get_imgs(reference,prhs[0]);
    // input 2: current image    
    get_imgs(current,prhs[1]);
    // input 3: ROI
    get_rois(roi,prhs[2]);
    // input 4: seedinfo
    get_seedinfo(seedinfo,prhs[3]);
    // input 5: thread diagram
    get_integer_array(threaddiagram,prhs[4]);
    // input 6: radius
    get_integer_scalar(radius,prhs[5]);
    // input 7: spacing
    get_integer_scalar(spacing,prhs[6]);
    // input 8: diffnorm cutoff
    get_double_scalar(cutoff_diffnorm,prhs[7]);
    // input 9: iteration cutoff
    get_integer_scalar(cutoff_iteration,prhs[8]);
    // input 10: subsettrunc
    get_logical_scalar(subsettrunc,prhs[9]);
    // input 11: image number
    get_integer_scalar(num_img,prhs[10]);
    // input 12: total number of images being processed
    get_integer_scalar(total_imgs,prhs[11]);
        
    // Check inputs - These are very basic checks, inconsistent inputs can possibly cause program to crash
    if (reference[0].gs.width == roi[0].mask.width && 
        reference[0].gs.height == roi[0].mask.height &&
        threaddiagram.height == (int)ceil((double)reference[0].gs.height/((double)(spacing+1))) &&
        threaddiagram.width == (int)ceil((double)reference[0].gs.width/((double)(spacing+1))) &&
        radius > 0 &&
        spacing >= 0 &&
        cutoff_diffnorm > 0.0 &&
        cutoff_iteration > 0 && 
        num_img >= 0 &&
        total_imgs > 0) {
        
        // Get total number of regions and threads from seedinfo --------//
        total_threads = (int)mxGetM(prhs[3]);
        total_region = (int)mxGetN(prhs[3]);
        
        // OpenMP Setup -------------------------------------------------//
        #ifdef NCORR_OPENMP
            // Set number of threads
            omp_set_num_threads(total_threads);
        #endif
                
        // Set cirroi ---------------------------------------------------//
        // Must set one for each thread
        roi[0].set_cirroi(radius,total_threads);
                        
        // Calc points --------------------------------------------------//
        // This buffer keeps track of which points have been analyzed already
        // Dimensions are the same as the thread diagram
        plot_calcpoints.alloc(threaddiagram.height,threaddiagram.width);
        
        // Form Inverse Compositional Buffers ---------------------------//     
        try { 
            // Thread independent:
            df_dx_buffer.resize(reference[0].gs.height*reference[0].gs.width,0.0);
            df_dy_buffer.resize(reference[0].gs.height*reference[0].gs.width,0.0);
            QK_B_QKT_buffer.resize(36*(current[0].bcoef.height-5)*(current[0].bcoef.width-5),0.0); // Very large
            
            // Thread dependent:
            g_buffer.resize(total_threads); 
            df_dp_buffer.resize(total_threads); 
            x_vec_buffer.resize(total_threads); 
            y_vec_buffer.resize(total_threads); 
            gradient_buffer.resize(total_threads); 
            hessian_gn_buffer.resize(total_threads); 
            for (int i=0; i<total_threads; i++) {
                g_buffer[i].resize((radius*2+1)*(radius*2+1),0.0); 
                df_dp_buffer[i].resize((radius*2+1)*(radius*2+1)*6,0.0); 
                x_vec_buffer[i].resize(6,0.0); 
                y_vec_buffer[i].resize(6,0.0); 
                gradient_buffer[i].resize(6,0.0); 
                hessian_gn_buffer[i].resize(36,0.0); 
            }
        } catch (std::exception&) {
            // Program ran out of memory, probably while trying to allocate the interpolation look up table.
            // Use mexErrMsgTxt to return although matlab will handle the exception.
            // This is thread safe since this is single threaded up to this point
            mexErrMsgTxt("Out of memory.\n");
        }
        
        // Form/set outputs ---------------------------------------------//
        // output 1: plot_disp
        // Form displacement structure
        mwSize dims[2] = {1,1};
        int numfields = 4;
        const char *fieldnames[] = {"plot_u","plot_v","plot_corrcoef","plot_validpoints"};
        plhs[0] = mxCreateStructArray(2, dims, numfields, fieldnames);

        // Form fields 
        // Same dimensions as thread diagram
        mxArray *mat_plot_u = mxCreateDoubleMatrix(threaddiagram.height, threaddiagram.width, mxREAL); 
        mxArray *mat_plot_v = mxCreateDoubleMatrix(threaddiagram.height, threaddiagram.width, mxREAL); 
        mxArray *mat_plot_corrcoef = mxCreateDoubleMatrix(threaddiagram.height, threaddiagram.width, mxREAL); 
        mxArray *mat_plot_validpoints = mxCreateLogicalMatrix(threaddiagram.height, threaddiagram.width); 

        // Add fields to structure
        // add u:
        mxSetFieldByNumber(plhs[0],0,0,mat_plot_u);
        // add v:
        mxSetFieldByNumber(plhs[0],0,1,mat_plot_v);
        // add c:
        mxSetFieldByNumber(plhs[0],0,2,mat_plot_corrcoef);
        // add plot_validpoints:
        mxSetFieldByNumber(plhs[0],0,3,mat_plot_validpoints);

        // output 2: outstate
        plhs[1] = mxCreateDoubleMatrix(1,1,mxREAL);
        
        // Get outputs --------------------------------------------------//
        // output 1: plot_disp
        // u:
        get_double_array(plot_u,mat_plot_u);
        // v:
        get_double_array(plot_v,mat_plot_v);
        // corrcoef
        get_double_array(plot_corrcoef,mat_plot_corrcoef);
        // plot_validpoints:
        get_logical_array(plot_validpoints,mat_plot_validpoints);
        // output 2: outstate
        outstate = mxGetPr(plhs[1]);
    } else {
        // Thread safe because it is single threaded up to this point
        mexErrMsgTxt("One of the inputs is incorrect.\n");
    }
}

// ----------------------------------------------------//
// Main Class Methods ---------------------------------//
// ----------------------------------------------------//

void class_rgdic::analysis() {     
    // Initialize outstate to cancelled
    *outstate = (double)CANCELLED; 
    
    // Waitbar Setup -----------------------------------------------//
    // Get thread with the largest number of compute points; this will be the thread that the waitbar's progress will be based on    
    int max_computepoints = 0;
    for (int i=0; i<total_threads; i++) {
        // i is the thread we're looking for; scan over seeds    
        int max_computepoints_buffer = 0;    
        for (int j=0; j<total_region; j++) {
            for (int k=0; k<total_threads; k++) {
                if (seedinfo[k+j*total_threads].num_thread == i) {
                    max_computepoints_buffer += seedinfo[k+j*total_threads].computepoints;
                    break;
                }
            }
        }
        if (max_computepoints_buffer > max_computepoints) {
            // Specify which thread number the waitbar's progress corresponds to
            num_thread_waitbar = i;
            max_computepoints = max_computepoints_buffer;
        }
    }    
    // Form waitbar - NOTE: ONLY THREAD 0 WILL UPDATE WAITBAR, BUT WILL UPDATE BASED ON num_thread_waitbar's PROGRESS
    waitbar.start(num_img,total_imgs,max_computepoints); 
    
    // Begin Analysis ----------------------------------------------//
    // Compute Lookup table and compute reference gradients
     precompute();
        
    // Start parallel region for computation
    bool abort = false; // This will be exit criterion for while loops
    
    // Enter parallel region - anything inside here needs to be threadsafe
    #ifdef NCORR_OPENMP
    #pragma omp parallel 
    {       
    #endif
            
        #ifdef NCORR_OPENMP
            // Get thread number 
            int num_thread = omp_get_thread_num();
        #else
            // Set to zero if openmp is not supported
            int num_thread = 0;
        #endif
                
        // Process one region at a time; Order does not matter
        for (int i=0; i<total_region; i++) {        
            // Find seed corresponding to num_thread
            int lind_seed = -1;
            for (int j=0; j<total_threads; j++) {
                if (seedinfo[j+i*total_threads].num_thread == num_thread) {
                    // Seed corresponding to this thread was found
                    lind_seed = j+i*total_threads;
                    break;
                }
            }
                
            // Update cirroi ---------------------------------------//
            // Need to update cirroi - This is threadsafe
            roi[0].update_cirroi(seedinfo[lind_seed].num_region,num_thread);

            // Initialize queue
            heap queue;
            
            // Add seed to queue
            std::vector<double> paramvector_seed(9,0); // [x y u v du/dx du/dy dv/dx dv/dy corrcoef]
            paramvector_seed[0] = seedinfo[lind_seed].paramvector.value[0];
            paramvector_seed[1] = seedinfo[lind_seed].paramvector.value[1];
            paramvector_seed[2] = seedinfo[lind_seed].paramvector.value[2];
            paramvector_seed[3] = seedinfo[lind_seed].paramvector.value[3];
            paramvector_seed[4] = seedinfo[lind_seed].paramvector.value[4];
            paramvector_seed[5] = seedinfo[lind_seed].paramvector.value[5];
            paramvector_seed[6] = seedinfo[lind_seed].paramvector.value[6];
            paramvector_seed[7] = seedinfo[lind_seed].paramvector.value[7];
            paramvector_seed[8] = seedinfo[lind_seed].paramvector.value[8];
            queue.push(paramvector_seed);

            // Inactivate seed point    
            int x_seed_reduced = (int)paramvector_seed[0]/(spacing+1); // x_seed and y_seed are guaranteed divisible by (spacing+1)
            int y_seed_reduced = (int)paramvector_seed[1]/(spacing+1);
            plot_calcpoints.value[y_seed_reduced+x_seed_reduced*plot_calcpoints.height] = true; 
            plot_validpoints.value[y_seed_reduced+x_seed_reduced*plot_validpoints.height] = true;

            // Enter While Loop - Exit when queue is empty
            while (!queue.empty()) {
                // 1) Load point with lowest correlation coefficient from queue 
                // 2) Delete point from queue
                // 3) Add data to plots 
                // 4) Analyze four surrounding points and sort 
                                
                // Step 1: load
                std::vector<double> paramvector_init = queue.top();

                // Step 2: delete 
                queue.pop();

                // Step 3: add data to plots
                int x_init_reduced = (int)paramvector_init[0]/(spacing+1);
                int y_init_reduced = (int)paramvector_init[1]/(spacing+1);
                plot_u.value[y_init_reduced+x_init_reduced*plot_u.height] = paramvector_init[2];
                plot_v.value[y_init_reduced+x_init_reduced*plot_v.height] = paramvector_init[3];
                plot_corrcoef.value[y_init_reduced+x_init_reduced*plot_corrcoef.height] = paramvector_init[8];
                                                    
                // Step 4: analyze four surrounding points - must increment by spacing parameter
                analyzepoint(queue,(int)paramvector_init[0],(int)paramvector_init[1]-(spacing+1),paramvector_init,seedinfo[lind_seed].num_region,num_thread);  
                analyzepoint(queue,(int)paramvector_init[0]+(spacing+1),(int)paramvector_init[1],paramvector_init,seedinfo[lind_seed].num_region,num_thread); 
                analyzepoint(queue,(int)paramvector_init[0],(int)paramvector_init[1]+(spacing+1),paramvector_init,seedinfo[lind_seed].num_region,num_thread); 
                analyzepoint(queue,(int)paramvector_init[0]-(spacing+1),(int)paramvector_init[1],paramvector_init,seedinfo[lind_seed].num_region,num_thread); 
                
                // Update and check waitbar - must only do this with thread 0
                if (num_thread == 0) {
                    if (!waitbar.updateandcheck()) {
                        // Waitbar was cancelled - set abort equal to true instead of 
                        // returning, this will break all threads from their computation 
                        // if they havent completed yet
                        abort = true;
                    }
                }                
                // Check if aborted
                if (abort) {
                    break;
                }
            }
            // Check if aborted
            if (abort) {
                break;
            }
        }    
        
    #ifdef NCORR_OPENMP
    }    
    #endif
       
    // At this point analysis has been completed successfully
    if (!abort) {
        *outstate = (double)SUCCESS; 
    }
}

void class_rgdic::precompute() {
    // Pre compute interpolation coefficients ---------------------------//
    for (int i=0; i<current[0].bcoef.width-5; i++) {
        for (int j=0; j<current[0].bcoef.height-5; j++) {                
            // Get bspline coefficients
            double b0 = current[0].bcoef.value[(j)+(i)*current[0].bcoef.height];
            double b1 = current[0].bcoef.value[(j+1)+(i)*current[0].bcoef.height];
            double b2 = current[0].bcoef.value[(j+2)+(i)*current[0].bcoef.height];
            double b3 = current[0].bcoef.value[(j+3)+(i)*current[0].bcoef.height];
            double b4 = current[0].bcoef.value[(j+4)+(i)*current[0].bcoef.height];
            double b5 = current[0].bcoef.value[(j+5)+(i)*current[0].bcoef.height];
            double b6 = current[0].bcoef.value[(j)+(i+1)*current[0].bcoef.height];
            double b7 = current[0].bcoef.value[(j+1)+(i+1)*current[0].bcoef.height];
            double b8 = current[0].bcoef.value[(j+2)+(i+1)*current[0].bcoef.height];
            double b9 = current[0].bcoef.value[(j+3)+(i+1)*current[0].bcoef.height];
            double b10 = current[0].bcoef.value[(j+4)+(i+1)*current[0].bcoef.height];
            double b11 = current[0].bcoef.value[(j+5)+(i+1)*current[0].bcoef.height];
            double b12 = current[0].bcoef.value[(j)+(i+2)*current[0].bcoef.height];
            double b13 = current[0].bcoef.value[(j+1)+(i+2)*current[0].bcoef.height];
            double b14 = current[0].bcoef.value[(j+2)+(i+2)*current[0].bcoef.height];
            double b15 = current[0].bcoef.value[(j+3)+(i+2)*current[0].bcoef.height];
            double b16 = current[0].bcoef.value[(j+4)+(i+2)*current[0].bcoef.height];
            double b17 = current[0].bcoef.value[(j+5)+(i+2)*current[0].bcoef.height];
            double b18 = current[0].bcoef.value[(j)+(i+3)*current[0].bcoef.height];
            double b19 = current[0].bcoef.value[(j+1)+(i+3)*current[0].bcoef.height];
            double b20 = current[0].bcoef.value[(j+2)+(i+3)*current[0].bcoef.height];
            double b21 = current[0].bcoef.value[(j+3)+(i+3)*current[0].bcoef.height];
            double b22 = current[0].bcoef.value[(j+4)+(i+3)*current[0].bcoef.height];
            double b23 = current[0].bcoef.value[(j+5)+(i+3)*current[0].bcoef.height];
            double b24 = current[0].bcoef.value[(j)+(i+4)*current[0].bcoef.height];
            double b25 = current[0].bcoef.value[(j+1)+(i+4)*current[0].bcoef.height];
            double b26 = current[0].bcoef.value[(j+2)+(i+4)*current[0].bcoef.height];
            double b27 = current[0].bcoef.value[(j+3)+(i+4)*current[0].bcoef.height];
            double b28 = current[0].bcoef.value[(j+4)+(i+4)*current[0].bcoef.height];
            double b29 = current[0].bcoef.value[(j+5)+(i+4)*current[0].bcoef.height];
            double b30 = current[0].bcoef.value[(j)+(i+5)*current[0].bcoef.height];
            double b31 = current[0].bcoef.value[(j+1)+(i+5)*current[0].bcoef.height];
            double b32 = current[0].bcoef.value[(j+2)+(i+5)*current[0].bcoef.height];
            double b33 = current[0].bcoef.value[(j+3)+(i+5)*current[0].bcoef.height];
            double b34 = current[0].bcoef.value[(j+4)+(i+5)*current[0].bcoef.height];
            double b35 = current[0].bcoef.value[(j+5)+(i+5)*current[0].bcoef.height];

            //Compute base index
            int lind_qkbqkt = j*(36)+i*(36*(current[0].bcoef.height-5));

            //Compute QK_B_QKT vector
            QK_B_QKT_buffer[lind_qkbqkt] = 0.00006944444444444444*b0+0.001805555555555556*b1+0.001805555555555556*b10+0.004583333333333333*b12+0.1191666666666667*b13+0.3025*b14+0.1191666666666667*b15+0.004583333333333333*b16+0.001805555555555556*b18+0.04694444444444444*b19+0.004583333333333333*b2+0.1191666666666667*b20+0.04694444444444444*b21+0.001805555555555556*b22+0.00006944444444444444*b24+0.001805555555555556*b25+0.004583333333333333*b26+0.001805555555555556*b27+0.00006944444444444444*b28+0.001805555555555556*b3+0.00006944444444444444*b4+0.001805555555555556*b6+0.04694444444444444*b7+0.1191666666666667*b8+0.04694444444444444*b9;
            QK_B_QKT_buffer[lind_qkbqkt+1] = 0.009027777777777778*b10-0.003472222222222222*b1-0.0003472222222222222*b0-0.02291666666666667*b12-0.2291666666666667*b13+0.2291666666666667*b15+0.02291666666666667*b16-0.009027777777777778*b18-0.09027777777777778*b19+0.09027777777777778*b21+0.009027777777777778*b22-0.0003472222222222222*b24-0.003472222222222222*b25+0.003472222222222222*b27+0.0003472222222222222*b28+0.003472222222222222*b3+0.0003472222222222222*b4-0.009027777777777778*b6-0.09027777777777778*b7+0.09027777777777778*b9;
            QK_B_QKT_buffer[lind_qkbqkt+2] = 0.0006944444444444444*b0+0.001388888888888889*b1+0.01805555555555556*b10+0.04583333333333333*b12+0.09166666666666667*b13-0.275*b14+0.09166666666666667*b15+0.04583333333333333*b16+0.01805555555555556*b18+0.03611111111111111*b19-0.004166666666666667*b2-0.1083333333333333*b20+0.03611111111111111*b21+0.01805555555555556*b22+0.0006944444444444444*b24+0.001388888888888889*b25-0.004166666666666667*b26+0.001388888888888889*b27+0.0006944444444444444*b28+0.001388888888888889*b3+0.0006944444444444444*b4+0.01805555555555556*b6+0.03611111111111111*b7-0.1083333333333333*b8+0.03611111111111111*b9;
            QK_B_QKT_buffer[lind_qkbqkt+3] = 0.001388888888888889*b1-0.0006944444444444444*b0+0.01805555555555556*b10-0.04583333333333333*b12+0.09166666666666667*b13-0.09166666666666667*b15+0.04583333333333333*b16-0.01805555555555556*b18+0.03611111111111111*b19-0.03611111111111111*b21+0.01805555555555556*b22-0.0006944444444444444*b24+0.001388888888888889*b25-0.001388888888888889*b27+0.0006944444444444444*b28-0.001388888888888889*b3+0.0006944444444444444*b4-0.01805555555555556*b6+0.03611111111111111*b7-0.03611111111111111*b9;
            QK_B_QKT_buffer[lind_qkbqkt+4] = 0.0003472222222222222*b0-0.001388888888888889*b1+0.009027777777777778*b10+0.02291666666666667*b12-0.09166666666666667*b13+0.1375*b14-0.09166666666666667*b15+0.02291666666666667*b16+0.009027777777777778*b18-0.03611111111111111*b19+0.002083333333333333*b2+0.05416666666666667*b20-0.03611111111111111*b21+0.009027777777777778*b22+0.0003472222222222222*b24-0.001388888888888889*b25+0.002083333333333333*b26-0.001388888888888889*b27+0.0003472222222222222*b28-0.001388888888888889*b3+0.0003472222222222222*b4+0.009027777777777778*b6-0.03611111111111111*b7+0.05416666666666667*b8-0.03611111111111111*b9;
            QK_B_QKT_buffer[lind_qkbqkt+5] = 0.0003472222222222222*b1-0.00006944444444444444*b0-0.009027777777777778*b10+0.001805555555555556*b11-0.004583333333333333*b12+0.02291666666666667*b13-0.04583333333333333*b14+0.04583333333333333*b15-0.02291666666666667*b16+0.004583333333333333*b17-0.001805555555555556*b18+0.009027777777777778*b19-0.0006944444444444444*b2-0.01805555555555556*b20+0.01805555555555556*b21-0.009027777777777778*b22+0.001805555555555556*b23-0.00006944444444444444*b24+0.0003472222222222222*b25-0.0006944444444444444*b26+0.0006944444444444444*b27-0.0003472222222222222*b28+0.00006944444444444444*b29+0.0006944444444444444*b3-0.0003472222222222222*b4+0.00006944444444444444*b5-0.001805555555555556*b6+0.009027777777777778*b7-0.01805555555555556*b8+0.01805555555555556*b9;
            QK_B_QKT_buffer[lind_qkbqkt+6] = 0.003472222222222222*b18-0.009027777777777778*b1-0.003472222222222222*b10-0.0003472222222222222*b0+0.09027777777777778*b19-0.02291666666666667*b2+0.2291666666666667*b20+0.09027777777777778*b21+0.003472222222222222*b22+0.0003472222222222222*b24+0.009027777777777778*b25+0.02291666666666667*b26+0.009027777777777778*b27+0.0003472222222222222*b28-0.009027777777777778*b3-0.0003472222222222222*b4-0.003472222222222222*b6-0.09027777777777778*b7-0.2291666666666667*b8-0.09027777777777778*b9;
            QK_B_QKT_buffer[lind_qkbqkt+7] = 0.001736111111111111*b0+0.01736111111111111*b1-0.01736111111111111*b10-0.01736111111111111*b18-0.1736111111111111*b19+0.1736111111111111*b21+0.01736111111111111*b22-0.001736111111111111*b24-0.01736111111111111*b25+0.01736111111111111*b27+0.001736111111111111*b28-0.01736111111111111*b3-0.001736111111111111*b4+0.01736111111111111*b6+0.1736111111111111*b7-0.1736111111111111*b9;
            QK_B_QKT_buffer[lind_qkbqkt+8] = 0.03472222222222222*b18-0.006944444444444444*b1-0.03472222222222222*b10-0.003472222222222222*b0+0.06944444444444444*b19+0.02083333333333333*b2-0.2083333333333333*b20+0.06944444444444444*b21+0.03472222222222222*b22+0.003472222222222222*b24+0.006944444444444444*b25-0.02083333333333333*b26+0.006944444444444444*b27+0.003472222222222222*b28-0.006944444444444444*b3-0.003472222222222222*b4-0.03472222222222222*b6-0.06944444444444444*b7+0.2083333333333333*b8-0.06944444444444444*b9;
            QK_B_QKT_buffer[lind_qkbqkt+9] = 0.003472222222222222*b0-0.006944444444444444*b1-0.03472222222222222*b10-0.03472222222222222*b18+0.06944444444444444*b19-0.06944444444444444*b21+0.03472222222222222*b22-0.003472222222222222*b24+0.006944444444444444*b25-0.006944444444444444*b27+0.003472222222222222*b28+0.006944444444444444*b3-0.003472222222222222*b4+0.03472222222222222*b6-0.06944444444444444*b7+0.06944444444444444*b9;
            QK_B_QKT_buffer[lind_qkbqkt+10] = 0.006944444444444444*b1-0.001736111111111111*b0-0.01736111111111111*b10+0.01736111111111111*b18-0.06944444444444444*b19-0.01041666666666667*b2+0.1041666666666667*b20-0.06944444444444444*b21+0.01736111111111111*b22+0.001736111111111111*b24-0.006944444444444444*b25+0.01041666666666667*b26-0.006944444444444444*b27+0.001736111111111111*b28+0.006944444444444444*b3-0.001736111111111111*b4-0.01736111111111111*b6+0.06944444444444444*b7-0.1041666666666667*b8+0.06944444444444444*b9;
            QK_B_QKT_buffer[lind_qkbqkt+11] = 0.0003472222222222222*b0-0.001736111111111111*b1+0.01736111111111111*b10-0.003472222222222222*b11-0.003472222222222222*b18+0.01736111111111111*b19+0.003472222222222222*b2-0.03472222222222222*b20+0.03472222222222222*b21-0.01736111111111111*b22+0.003472222222222222*b23-0.0003472222222222222*b24+0.001736111111111111*b25-0.003472222222222222*b26+0.003472222222222222*b27-0.001736111111111111*b28+0.0003472222222222222*b29-0.003472222222222222*b3+0.001736111111111111*b4-0.0003472222222222222*b5+0.003472222222222222*b6-0.01736111111111111*b7+0.03472222222222222*b8-0.03472222222222222*b9;
            QK_B_QKT_buffer[lind_qkbqkt+12] = 0.0006944444444444444*b0+0.01805555555555556*b1+0.001388888888888889*b10-0.004166666666666667*b12-0.1083333333333333*b13-0.275*b14-0.1083333333333333*b15-0.004166666666666667*b16+0.001388888888888889*b18+0.03611111111111111*b19+0.04583333333333333*b2+0.09166666666666667*b20+0.03611111111111111*b21+0.001388888888888889*b22+0.0006944444444444444*b24+0.01805555555555556*b25+0.04583333333333333*b26+0.01805555555555556*b27+0.0006944444444444444*b28+0.01805555555555556*b3+0.0006944444444444444*b4+0.001388888888888889*b6+0.03611111111111111*b7+0.09166666666666667*b8+0.03611111111111111*b9;
            QK_B_QKT_buffer[lind_qkbqkt+13] = 0.006944444444444444*b10-0.03472222222222222*b1-0.003472222222222222*b0+0.02083333333333333*b12+0.2083333333333333*b13-0.2083333333333333*b15-0.02083333333333333*b16-0.006944444444444444*b18-0.06944444444444444*b19+0.06944444444444444*b21+0.006944444444444444*b22-0.003472222222222222*b24-0.03472222222222222*b25+0.03472222222222222*b27+0.003472222222222222*b28+0.03472222222222222*b3+0.003472222222222222*b4-0.006944444444444444*b6-0.06944444444444444*b7+0.06944444444444444*b9;
            QK_B_QKT_buffer[lind_qkbqkt+14] = 0.006944444444444444*b0+0.01388888888888889*b1+0.01388888888888889*b10-0.04166666666666667*b12-0.08333333333333333*b13+0.25*b14-0.08333333333333333*b15-0.04166666666666667*b16+0.01388888888888889*b18+0.02777777777777778*b19-0.04166666666666667*b2-0.08333333333333333*b20+0.02777777777777778*b21+0.01388888888888889*b22+0.006944444444444444*b24+0.01388888888888889*b25-0.04166666666666667*b26+0.01388888888888889*b27+0.006944444444444444*b28+0.01388888888888889*b3+0.006944444444444444*b4+0.01388888888888889*b6+0.02777777777777778*b7-0.08333333333333333*b8+0.02777777777777778*b9;
            QK_B_QKT_buffer[lind_qkbqkt+15] = 0.01388888888888889*b1-0.006944444444444444*b0+0.01388888888888889*b10+0.04166666666666667*b12-0.08333333333333333*b13+0.08333333333333333*b15-0.04166666666666667*b16-0.01388888888888889*b18+0.02777777777777778*b19-0.02777777777777778*b21+0.01388888888888889*b22-0.006944444444444444*b24+0.01388888888888889*b25-0.01388888888888889*b27+0.006944444444444444*b28-0.01388888888888889*b3+0.006944444444444444*b4-0.01388888888888889*b6+0.02777777777777778*b7-0.02777777777777778*b9;
            QK_B_QKT_buffer[lind_qkbqkt+16] = 0.003472222222222222*b0-0.01388888888888889*b1+0.006944444444444444*b10-0.02083333333333333*b12+0.08333333333333333*b13-0.125*b14+0.08333333333333333*b15-0.02083333333333333*b16+0.006944444444444444*b18-0.02777777777777778*b19+0.02083333333333333*b2+0.04166666666666667*b20-0.02777777777777778*b21+0.006944444444444444*b22+0.003472222222222222*b24-0.01388888888888889*b25+0.02083333333333333*b26-0.01388888888888889*b27+0.003472222222222222*b28-0.01388888888888889*b3+0.003472222222222222*b4+0.006944444444444444*b6-0.02777777777777778*b7+0.04166666666666667*b8-0.02777777777777778*b9;
            QK_B_QKT_buffer[lind_qkbqkt+17] = 0.003472222222222222*b1-0.0006944444444444444*b0-0.006944444444444444*b10+0.001388888888888889*b11+0.004166666666666667*b12-0.02083333333333333*b13+0.04166666666666667*b14-0.04166666666666667*b15+0.02083333333333333*b16-0.004166666666666667*b17-0.001388888888888889*b18+0.006944444444444444*b19-0.006944444444444444*b2-0.01388888888888889*b20+0.01388888888888889*b21-0.006944444444444444*b22+0.001388888888888889*b23-0.0006944444444444444*b24+0.003472222222222222*b25-0.006944444444444444*b26+0.006944444444444444*b27-0.003472222222222222*b28+0.0006944444444444444*b29+0.006944444444444444*b3-0.003472222222222222*b4+0.0006944444444444444*b5-0.001388888888888889*b6+0.006944444444444444*b7-0.01388888888888889*b8+0.01388888888888889*b9;
            QK_B_QKT_buffer[lind_qkbqkt+18] = 0.001388888888888889*b10-0.01805555555555556*b1-0.0006944444444444444*b0-0.001388888888888889*b18-0.03611111111111111*b19-0.04583333333333333*b2-0.09166666666666667*b20-0.03611111111111111*b21-0.001388888888888889*b22+0.0006944444444444444*b24+0.01805555555555556*b25+0.04583333333333333*b26+0.01805555555555556*b27+0.0006944444444444444*b28-0.01805555555555556*b3-0.0006944444444444444*b4+0.001388888888888889*b6+0.03611111111111111*b7+0.09166666666666667*b8+0.03611111111111111*b9;
            QK_B_QKT_buffer[lind_qkbqkt+19] = 0.003472222222222222*b0+0.03472222222222222*b1+0.006944444444444444*b10+0.006944444444444444*b18+0.06944444444444444*b19-0.06944444444444444*b21-0.006944444444444444*b22-0.003472222222222222*b24-0.03472222222222222*b25+0.03472222222222222*b27+0.003472222222222222*b28-0.03472222222222222*b3-0.003472222222222222*b4-0.006944444444444444*b6-0.06944444444444444*b7+0.06944444444444444*b9;
            QK_B_QKT_buffer[lind_qkbqkt+20] = 0.01388888888888889*b10-0.01388888888888889*b1-0.006944444444444444*b0-0.01388888888888889*b18-0.02777777777777778*b19+0.04166666666666667*b2+0.08333333333333333*b20-0.02777777777777778*b21-0.01388888888888889*b22+0.006944444444444444*b24+0.01388888888888889*b25-0.04166666666666667*b26+0.01388888888888889*b27+0.006944444444444444*b28-0.01388888888888889*b3-0.006944444444444444*b4+0.01388888888888889*b6+0.02777777777777778*b7-0.08333333333333333*b8+0.02777777777777778*b9;
            QK_B_QKT_buffer[lind_qkbqkt+21] = 0.006944444444444444*b0-0.01388888888888889*b1+0.01388888888888889*b10+0.01388888888888889*b18-0.02777777777777778*b19+0.02777777777777778*b21-0.01388888888888889*b22-0.006944444444444444*b24+0.01388888888888889*b25-0.01388888888888889*b27+0.006944444444444444*b28+0.01388888888888889*b3-0.006944444444444444*b4-0.01388888888888889*b6+0.02777777777777778*b7-0.02777777777777778*b9;
            QK_B_QKT_buffer[lind_qkbqkt+22] = 0.01388888888888889*b1-0.003472222222222222*b0+0.006944444444444444*b10-0.006944444444444444*b18+0.02777777777777778*b19-0.02083333333333333*b2-0.04166666666666667*b20+0.02777777777777778*b21-0.006944444444444444*b22+0.003472222222222222*b24-0.01388888888888889*b25+0.02083333333333333*b26-0.01388888888888889*b27+0.003472222222222222*b28+0.01388888888888889*b3-0.003472222222222222*b4+0.006944444444444444*b6-0.02777777777777778*b7+0.04166666666666667*b8-0.02777777777777778*b9;
            QK_B_QKT_buffer[lind_qkbqkt+23] = 0.0006944444444444444*b0-0.003472222222222222*b1-0.006944444444444444*b10+0.001388888888888889*b11+0.001388888888888889*b18-0.006944444444444444*b19+0.006944444444444444*b2+0.01388888888888889*b20-0.01388888888888889*b21+0.006944444444444444*b22-0.001388888888888889*b23-0.0006944444444444444*b24+0.003472222222222222*b25-0.006944444444444444*b26+0.006944444444444444*b27-0.003472222222222222*b28+0.0006944444444444444*b29-0.006944444444444444*b3+0.003472222222222222*b4-0.0006944444444444444*b5-0.001388888888888889*b6+0.006944444444444444*b7-0.01388888888888889*b8+0.01388888888888889*b9;
            QK_B_QKT_buffer[lind_qkbqkt+24] = 0.0003472222222222222*b0+0.009027777777777778*b1-0.001388888888888889*b10+0.002083333333333333*b12+0.05416666666666667*b13+0.1375*b14+0.05416666666666667*b15+0.002083333333333333*b16-0.001388888888888889*b18-0.03611111111111111*b19+0.02291666666666667*b2-0.09166666666666667*b20-0.03611111111111111*b21-0.001388888888888889*b22+0.0003472222222222222*b24+0.009027777777777778*b25+0.02291666666666667*b26+0.009027777777777778*b27+0.0003472222222222222*b28+0.009027777777777778*b3+0.0003472222222222222*b4-0.001388888888888889*b6-0.03611111111111111*b7-0.09166666666666667*b8-0.03611111111111111*b9;
            QK_B_QKT_buffer[lind_qkbqkt+25] = 0.1041666666666667*b15-0.01736111111111111*b1-0.006944444444444444*b10-0.01041666666666667*b12-0.1041666666666667*b13-0.001736111111111111*b0+0.01041666666666667*b16+0.006944444444444444*b18+0.06944444444444444*b19-0.06944444444444444*b21-0.006944444444444444*b22-0.001736111111111111*b24-0.01736111111111111*b25+0.01736111111111111*b27+0.001736111111111111*b28+0.01736111111111111*b3+0.001736111111111111*b4+0.006944444444444444*b6+0.06944444444444444*b7-0.06944444444444444*b9;
            QK_B_QKT_buffer[lind_qkbqkt+26] = 0.003472222222222222*b0+0.006944444444444444*b1-0.01388888888888889*b10+0.02083333333333333*b12+0.04166666666666667*b13-0.125*b14+0.04166666666666667*b15+0.02083333333333333*b16-0.01388888888888889*b18-0.02777777777777778*b19-0.02083333333333333*b2+0.08333333333333333*b20-0.02777777777777778*b21-0.01388888888888889*b22+0.003472222222222222*b24+0.006944444444444444*b25-0.02083333333333333*b26+0.006944444444444444*b27+0.003472222222222222*b28+0.006944444444444444*b3+0.003472222222222222*b4-0.01388888888888889*b6-0.02777777777777778*b7+0.08333333333333333*b8-0.02777777777777778*b9;
            QK_B_QKT_buffer[lind_qkbqkt+27] = 0.006944444444444444*b1-0.003472222222222222*b0-0.01388888888888889*b10-0.02083333333333333*b12+0.04166666666666667*b13-0.04166666666666667*b15+0.02083333333333333*b16+0.01388888888888889*b18-0.02777777777777778*b19+0.02777777777777778*b21-0.01388888888888889*b22-0.003472222222222222*b24+0.006944444444444444*b25-0.006944444444444444*b27+0.003472222222222222*b28-0.006944444444444444*b3+0.003472222222222222*b4+0.01388888888888889*b6-0.02777777777777778*b7+0.02777777777777778*b9;
            QK_B_QKT_buffer[lind_qkbqkt+28] = 0.001736111111111111*b0-0.006944444444444444*b1-0.006944444444444444*b10+0.01041666666666667*b12-0.04166666666666667*b13+0.0625*b14-0.04166666666666667*b15+0.01041666666666667*b16-0.006944444444444444*b18+0.02777777777777778*b19+0.01041666666666667*b2-0.04166666666666667*b20+0.02777777777777778*b21-0.006944444444444444*b22+0.001736111111111111*b24-0.006944444444444444*b25+0.01041666666666667*b26-0.006944444444444444*b27+0.001736111111111111*b28-0.006944444444444444*b3+0.001736111111111111*b4-0.006944444444444444*b6+0.02777777777777778*b7-0.04166666666666667*b8+0.02777777777777778*b9;
            QK_B_QKT_buffer[lind_qkbqkt+29] = 0.001736111111111111*b1-0.0003472222222222222*b0+0.006944444444444444*b10-0.001388888888888889*b11-0.002083333333333333*b12+0.01041666666666667*b13-0.02083333333333333*b14+0.02083333333333333*b15-0.01041666666666667*b16+0.002083333333333333*b17+0.001388888888888889*b18-0.006944444444444444*b19-0.003472222222222222*b2+0.01388888888888889*b20-0.01388888888888889*b21+0.006944444444444444*b22-0.001388888888888889*b23-0.0003472222222222222*b24+0.001736111111111111*b25-0.003472222222222222*b26+0.003472222222222222*b27-0.001736111111111111*b28+0.0003472222222222222*b29+0.003472222222222222*b3-0.001736111111111111*b4+0.0003472222222222222*b5+0.001388888888888889*b6-0.006944444444444444*b7+0.01388888888888889*b8-0.01388888888888889*b9;
            QK_B_QKT_buffer[lind_qkbqkt+30] = 0.0003472222222222222*b10-0.001805555555555556*b1-0.00006944444444444444*b0-0.0006944444444444444*b12-0.01805555555555556*b13-0.04583333333333333*b14-0.01805555555555556*b15-0.0006944444444444444*b16+0.0006944444444444444*b18+0.01805555555555556*b19-0.004583333333333333*b2+0.04583333333333333*b20+0.01805555555555556*b21+0.0006944444444444444*b22-0.0003472222222222222*b24-0.009027777777777778*b25-0.02291666666666667*b26-0.009027777777777778*b27-0.0003472222222222222*b28-0.001805555555555556*b3+0.00006944444444444444*b30+0.001805555555555556*b31+0.004583333333333333*b32+0.001805555555555556*b33+0.00006944444444444444*b34-0.00006944444444444444*b4+0.0003472222222222222*b6+0.009027777777777778*b7+0.02291666666666667*b8+0.009027777777777778*b9;
            QK_B_QKT_buffer[lind_qkbqkt+31] = 0.0003472222222222222*b0+0.003472222222222222*b1+0.001736111111111111*b10+0.003472222222222222*b12+0.03472222222222222*b13-0.03472222222222222*b15-0.003472222222222222*b16-0.003472222222222222*b18-0.03472222222222222*b19+0.03472222222222222*b21+0.003472222222222222*b22+0.001736111111111111*b24+0.01736111111111111*b25-0.01736111111111111*b27-0.001736111111111111*b28-0.003472222222222222*b3-0.0003472222222222222*b30-0.003472222222222222*b31+0.003472222222222222*b33+0.0003472222222222222*b34-0.0003472222222222222*b4-0.001736111111111111*b6-0.01736111111111111*b7+0.01736111111111111*b9;
            QK_B_QKT_buffer[lind_qkbqkt+32] = 0.003472222222222222*b10-0.001388888888888889*b1-0.0006944444444444444*b0-0.006944444444444444*b12-0.01388888888888889*b13+0.04166666666666667*b14-0.01388888888888889*b15-0.006944444444444444*b16+0.006944444444444444*b18+0.01388888888888889*b19+0.004166666666666667*b2-0.04166666666666667*b20+0.01388888888888889*b21+0.006944444444444444*b22-0.003472222222222222*b24-0.006944444444444444*b25+0.02083333333333333*b26-0.006944444444444444*b27-0.003472222222222222*b28-0.001388888888888889*b3+0.0006944444444444444*b30+0.001388888888888889*b31-0.004166666666666667*b32+0.001388888888888889*b33+0.0006944444444444444*b34-0.0006944444444444444*b4+0.003472222222222222*b6+0.006944444444444444*b7-0.02083333333333333*b8+0.006944444444444444*b9;
            QK_B_QKT_buffer[lind_qkbqkt+33] = 0.0006944444444444444*b0-0.001388888888888889*b1+0.003472222222222222*b10+0.006944444444444444*b12-0.01388888888888889*b13+0.01388888888888889*b15-0.006944444444444444*b16-0.006944444444444444*b18+0.01388888888888889*b19-0.01388888888888889*b21+0.006944444444444444*b22+0.003472222222222222*b24-0.006944444444444444*b25+0.006944444444444444*b27-0.003472222222222222*b28+0.001388888888888889*b3-0.0006944444444444444*b30+0.001388888888888889*b31-0.001388888888888889*b33+0.0006944444444444444*b34-0.0006944444444444444*b4-0.003472222222222222*b6+0.006944444444444444*b7-0.006944444444444444*b9;
            QK_B_QKT_buffer[lind_qkbqkt+34] = 0.001388888888888889*b1-0.0003472222222222222*b0+0.001736111111111111*b10-0.003472222222222222*b12+0.01388888888888889*b13-0.02083333333333333*b14+0.01388888888888889*b15-0.003472222222222222*b16+0.003472222222222222*b18-0.01388888888888889*b19-0.002083333333333333*b2+0.02083333333333333*b20-0.01388888888888889*b21+0.003472222222222222*b22-0.001736111111111111*b24+0.006944444444444444*b25-0.01041666666666667*b26+0.006944444444444444*b27-0.001736111111111111*b28+0.001388888888888889*b3+0.0003472222222222222*b30-0.001388888888888889*b31+0.002083333333333333*b32-0.001388888888888889*b33+0.0003472222222222222*b34-0.0003472222222222222*b4+0.001736111111111111*b6-0.006944444444444444*b7+0.01041666666666667*b8-0.006944444444444444*b9;
            QK_B_QKT_buffer[lind_qkbqkt+35] = 0.00006944444444444444*b0-0.0003472222222222222*b1-0.001736111111111111*b10+0.0003472222222222222*b11+0.0006944444444444444*b12-0.003472222222222222*b13+0.006944444444444444*b14-0.006944444444444444*b15+0.003472222222222222*b16-0.0006944444444444444*b17-0.0006944444444444444*b18+0.003472222222222222*b19+0.0006944444444444444*b2-0.006944444444444444*b20+0.006944444444444444*b21-0.003472222222222222*b22+0.0006944444444444444*b23+0.0003472222222222222*b24-0.001736111111111111*b25+0.003472222222222222*b26-0.003472222222222222*b27+0.001736111111111111*b28-0.0003472222222222222*b29-0.0006944444444444444*b3-0.00006944444444444444*b30+0.0003472222222222222*b31-0.0006944444444444444*b32+0.0006944444444444444*b33-0.0003472222222222222*b34+0.00006944444444444444*b35+0.0003472222222222222*b4-0.00006944444444444444*b5-0.0003472222222222222*b6+0.001736111111111111*b7-0.003472222222222222*b8+0.003472222222222222*b9;
        }
    }
    
    // Pre compute reference image gradient -----------------------------//
    for (int i=reference[0].border_bcoef-2; i<reference[0].border_bcoef+reference[0].gs.width-2; i++) {
        for (int j=reference[0].border_bcoef-2; j<reference[0].border_bcoef+reference[0].gs.height-2; j++) {                
            // Get bspline coefficients
            double b0 = reference[0].bcoef.value[(j)+(i)*reference[0].bcoef.height];
            double b1 = reference[0].bcoef.value[(j+1)+(i)*reference[0].bcoef.height];
            double b2 = reference[0].bcoef.value[(j+2)+(i)*reference[0].bcoef.height];
            double b3 = reference[0].bcoef.value[(j+3)+(i)*reference[0].bcoef.height];
            double b4 = reference[0].bcoef.value[(j+4)+(i)*reference[0].bcoef.height];
            double b5 = reference[0].bcoef.value[(j+5)+(i)*reference[0].bcoef.height];            
            double b6 = reference[0].bcoef.value[(j)+(i+1)*reference[0].bcoef.height];
            double b7 = reference[0].bcoef.value[(j+1)+(i+1)*reference[0].bcoef.height];
            double b8 = reference[0].bcoef.value[(j+2)+(i+1)*reference[0].bcoef.height];
            double b9 = reference[0].bcoef.value[(j+3)+(i+1)*reference[0].bcoef.height];
            double b11 = reference[0].bcoef.value[(j+4)+(i+1)*reference[0].bcoef.height];
            double b10 = reference[0].bcoef.value[(j+5)+(i+1)*reference[0].bcoef.height];            
            double b12 = reference[0].bcoef.value[(j)+(i+2)*reference[0].bcoef.height];
            double b13 = reference[0].bcoef.value[(j+1)+(i+2)*reference[0].bcoef.height];
            double b14 = reference[0].bcoef.value[(j+2)+(i+2)*reference[0].bcoef.height];            
            double b15 = reference[0].bcoef.value[(j+3)+(i+2)*reference[0].bcoef.height];
            double b16 = reference[0].bcoef.value[(j+4)+(i+2)*reference[0].bcoef.height];
            double b17 = reference[0].bcoef.value[(j+5)+(i+2)*reference[0].bcoef.height];            
            double b18 = reference[0].bcoef.value[(j)+(i+3)*reference[0].bcoef.height];
            double b19 = reference[0].bcoef.value[(j+1)+(i+3)*reference[0].bcoef.height];
            double b20 = reference[0].bcoef.value[(j+2)+(i+3)*reference[0].bcoef.height];
            double b21 = reference[0].bcoef.value[(j+3)+(i+3)*reference[0].bcoef.height];
            double b22 = reference[0].bcoef.value[(j+4)+(i+3)*reference[0].bcoef.height];
            double b23 = reference[0].bcoef.value[(j+5)+(i+3)*reference[0].bcoef.height];            
            double b24 = reference[0].bcoef.value[(j)+(i+4)*reference[0].bcoef.height];
            double b25 = reference[0].bcoef.value[(j+1)+(i+4)*reference[0].bcoef.height];
            double b26 = reference[0].bcoef.value[(j+2)+(i+4)*reference[0].bcoef.height];
            double b27 = reference[0].bcoef.value[(j+3)+(i+4)*reference[0].bcoef.height];
            double b28 = reference[0].bcoef.value[(j+4)+(i+4)*reference[0].bcoef.height];
            double b29 = reference[0].bcoef.value[(j+5)+(i+4)*reference[0].bcoef.height];
            double b30 = reference[0].bcoef.value[(j)+(i+5)*reference[0].bcoef.height];
            double b31 = reference[0].bcoef.value[(j+1)+(i+5)*reference[0].bcoef.height];
            double b32 = reference[0].bcoef.value[(j+2)+(i+5)*reference[0].bcoef.height];
            double b33 = reference[0].bcoef.value[(j+3)+(i+5)*reference[0].bcoef.height];
            double b34 = reference[0].bcoef.value[(j+4)+(i+5)*reference[0].bcoef.height];
            double b35 = reference[0].bcoef.value[(j+5)+(i+5)*reference[0].bcoef.height];
            
            // Compute base index
            int lind_f = (j-(reference[0].border_bcoef-2))+(i-(reference[0].border_bcoef-2))*reference[0].gs.height;

            // Compute Gradients using b-spline coefficients
            // First order
            df_dx_buffer[lind_f] = 0.003472222222222222*b18-0.009027777777777778*b1-0.003472222222222222*b10-0.0003472222222222222*b0+0.09027777777777778*b19-0.02291666666666667*b2+0.2291666666666667*b20+0.09027777777777778*b21+0.003472222222222222*b22+0.0003472222222222222*b24+0.009027777777777778*b25+0.02291666666666667*b26+0.009027777777777778*b27+0.0003472222222222222*b28-0.009027777777777778*b3-0.0003472222222222222*b4-0.003472222222222222*b6-0.09027777777777778*b7-0.2291666666666667*b8-0.09027777777777778*b9;
            df_dy_buffer[lind_f] = 0.009027777777777778*b10-0.003472222222222222*b1-0.0003472222222222222*b0-0.02291666666666667*b12-0.2291666666666667*b13+0.2291666666666667*b15+0.02291666666666667*b16-0.009027777777777778*b18-0.09027777777777778*b19+0.09027777777777778*b21+0.009027777777777778*b22-0.0003472222222222222*b24-0.003472222222222222*b25+0.003472222222222222*b27+0.0003472222222222222*b28+0.003472222222222222*b3+0.0003472222222222222*b4-0.009027777777777778*b6-0.09027777777777778*b7+0.09027777777777778*b9;
        }
    }
}

void class_rgdic::analyzepoint(heap &queue,const int &x,const int &y,const std::vector<double> &paramvector_init,const int &num_region,const int &num_thread) {
    // These are read only, so they are thread safe
    static double cutoff_corrcoef = 2.0;             // Heuristic, but 2.0 is pretty high. Range is [0,4]. Different images can have different corrcoef cutoffs which work well, so set this to a low value
    static double cutoff_disp = (double)spacing+1;   // Heuristic, this prevents large displacement jumps (most likely incorrect data) form being added
    
    // Reduce coordinates first
    int x_reduced = x/(spacing+1);
    int y_reduced = y/(spacing+1);
    
    // Make sure point is within region bounds first
    if (x >= roi[0].region[num_region].leftbound &&
        x <= roi[0].region[num_region].rightbound &&
        y >= roi[0].region[num_region].upperbound &&
        y <= roi[0].region[num_region].lowerbound &&
        !plot_calcpoints.value[y_reduced+x_reduced*plot_calcpoints.height] && 
        roi[0].withinregion(x,y,num_region) && // MUST TEST IF IN THE CORRECT REGION - its possible for two adjacent regions to have the same thread number
        threaddiagram.value[y_reduced+x_reduced*threaddiagram.height] == num_thread) {   
        // Initialize paramvector
        std::vector<double> paramvector(9,0); // [x y u v du/dx du/dy dv/dx dv/dy corrcoef]
        
        // Calculate paramvector for a point
        OUT outstate_calcpoint = calcpoint(paramvector,x,y,paramvector_init,num_region,num_thread);
        
        // Make sure parameters are correct before adding them to queue
        if (outstate_calcpoint == SUCCESS && 
            paramvector[8] < cutoff_corrcoef && 
            fabs(paramvector_init[2] - paramvector[2]) < cutoff_disp && 
            fabs(paramvector_init[3] - paramvector[3]) < cutoff_disp) {        
            // Insert paramvector based on correlation coefficient
            queue.push(paramvector);
            
            // Valid Point
            plot_validpoints.value[y_reduced+x_reduced*plot_validpoints.height] = true;    
        }
        // Calculated point
        plot_calcpoints.value[y_reduced+x_reduced*plot_calcpoints.height] = true;    

        // Increment waitbar - only use num_thread_waitbar to update these counters
        if (num_thread == num_thread_waitbar) {
            // Increment is threadsafe
            waitbar.increment();
        }
    }
}

OUT class_rgdic::calcpoint(std::vector<double> &paramvector,const int &x,const int &y,const std::vector<double> &paramvector_init,const int &num_region,const int &num_thread) {
    // Get cirroi -> Find initial guess -> Refine results with IC-GN -> Return true or false and store output
    // Step 1: Get cirroi
    roi[0].get_cirroi(x,y,num_region,subsettrunc,num_thread);        
        
    // Step 2: Get initial guess - Use displacement and displacement gradients to get initial guess
    std::vector<double> defvector_init(6,0); // [u v du/dx du/dy dv/dx dv/dy]
    // u_init = u+du/dx*x_delta+du/dy*y_delta;
    // v_init = v+dv/dx*x_delta+dv/dy*y_delta;
    defvector_init[0] = paramvector_init[2]+paramvector_init[4]*(x-paramvector_init[0])+paramvector_init[5]*(y-paramvector_init[1]); 
    defvector_init[1] = paramvector_init[3]+paramvector_init[6]*(x-paramvector_init[0])+paramvector_init[7]*(y-paramvector_init[1]); 
    defvector_init[2] = paramvector_init[4];
    defvector_init[3] = paramvector_init[5];
    defvector_init[4] = paramvector_init[6];
    defvector_init[5] = paramvector_init[7];
            
    // Step 3: Get refined results with IC-GN    
    std::vector<double> defvector(6,0); // [u v du/dx du/dy dv/dx dv/dy]
    double corrcoef;
    OUT outstate_iterative = iterativesearch(defvector,corrcoef,defvector_init,num_thread);            
    
    if (outstate_iterative == SUCCESS) {
        // Step 4: Store output and return true
        paramvector[0] = x;
        paramvector[1] = y;
        paramvector[2] = defvector[0];
        paramvector[3] = defvector[1];
        paramvector[4] = defvector[2];
        paramvector[5] = defvector[3];
        paramvector[6] = defvector[4];
        paramvector[7] = defvector[5];    
        paramvector[8] = corrcoef;
        return SUCCESS;
    }    
    return FAILED;        
}

OUT class_rgdic::iterativesearch(std::vector<double> &defvector,double &corrcoef,const std::vector<double> &defvector_init,const int &num_thread) {   
    // Calculate fm
    double fm = 0.0;
    for (int i=0; i<roi[0].cirroi[num_thread].region.noderange.height; i++) {
        for (int j=0; j<roi[0].cirroi[num_thread].region.noderange.value[i]; j+=2) {
            for (int k=roi[0].cirroi[num_thread].region.nodelist.value[i+j*roi[0].cirroi[num_thread].region.nodelist.height]; k<=roi[0].cirroi[num_thread].region.nodelist.value[i+(j+1)*roi[0].cirroi[num_thread].region.nodelist.height]; k++) {
                int lind_ref = k+(i+(roi[0].cirroi[num_thread].x-roi[0].cirroi[num_thread].radius))*reference[0].gs.height;
                fm += reference[0].gs.value[lind_ref];
            }
        }
    }
    fm = fm/((double)roi[0].cirroi[num_thread].region.totalpoints);
        
    // Calculate deltaf_inf
    double deltaf_inv = 0.0;
    for (int i=0; i<roi[0].cirroi[num_thread].region.noderange.height; i++) {
        for (int j=0; j<roi[0].cirroi[num_thread].region.noderange.value[i]; j+=2) {
            for (int k=roi[0].cirroi[num_thread].region.nodelist.value[i+j*roi[0].cirroi[num_thread].region.nodelist.height]; k<=roi[0].cirroi[num_thread].region.nodelist.value[i+(j+1)*roi[0].cirroi[num_thread].region.nodelist.height]; k++) {
                int lind_ref = k+(i+(roi[0].cirroi[num_thread].x-roi[0].cirroi[num_thread].radius))*reference[0].gs.height;
                deltaf_inv += pow(reference[0].gs.value[lind_ref]-fm,2);
            }
        }
    }
    deltaf_inv = sqrt(deltaf_inv);

    // check to make sure deltaf_inv (strictly positive) isn't close to zero; if it is, iterative search fails
    if (deltaf_inv > LAMBDA) {    
        // Finish deltaf_inv
        deltaf_inv = 1.0/deltaf_inv;
        
        // Precompute "Steepest descent images"                           
        for (int i=0; i<roi[0].cirroi[num_thread].region.noderange.height; i++) {
            for (int j=0; j<roi[0].cirroi[num_thread].region.noderange.value[i]; j+=2) {
                for (int k=roi[0].cirroi[num_thread].region.nodelist.value[i+j*roi[0].cirroi[num_thread].region.nodelist.height]; k<=roi[0].cirroi[num_thread].region.nodelist.value[i+(j+1)*roi[0].cirroi[num_thread].region.nodelist.height]; k++) {
                    // Find new coordinates
                    double dx = (double)(i-roi[0].cirroi[num_thread].radius);
                    double dy = (double)(k-roi[0].cirroi[num_thread].y);
                    
                    int y_tilda_floor = k;
                    int x_tilda_floor = i+(roi[0].cirroi[num_thread].x-roi[0].cirroi[num_thread].radius);

                    // Calculate lind_f for gradient and lind_df for the first order "steepest descent images"
                    int lind_f = y_tilda_floor+x_tilda_floor*reference[0].gs.height;
                    int lind_df = ((k-roi[0].cirroi[num_thread].y)+roi[0].cirroi[num_thread].radius)*6+i*(roi[0].cirroi[num_thread].region.nodelist.height*6);

                    // First order
                    df_dp_buffer[num_thread][lind_df] = df_dx_buffer[lind_f];      // u
                    df_dp_buffer[num_thread][lind_df+1] = df_dy_buffer[lind_f];    // v
                    df_dp_buffer[num_thread][lind_df+2] = df_dx_buffer[lind_f]*dx; // dudx
                    df_dp_buffer[num_thread][lind_df+3] = df_dx_buffer[lind_f]*dy; // dudy
                    df_dp_buffer[num_thread][lind_df+4] = df_dy_buffer[lind_f]*dx; // dvdx
                    df_dp_buffer[num_thread][lind_df+5] = df_dy_buffer[lind_f]*dy; // dvdy
                }
            }        
        }
                
        // Precompute GN hessian
        // Initialize to zero first
        std::fill(hessian_gn_buffer[num_thread].begin(),hessian_gn_buffer[num_thread].end(),0.0);
        for (int i=0; i<roi[0].cirroi[num_thread].region.noderange.height; i++) {
            for (int j=0; j<roi[0].cirroi[num_thread].region.noderange.value[i]; j+=2) {
                for (int k=roi[0].cirroi[num_thread].region.nodelist.value[i+j*roi[0].cirroi[num_thread].region.nodelist.height]; k<=roi[0].cirroi[num_thread].region.nodelist.value[i+(j+1)*roi[0].cirroi[num_thread].region.nodelist.height]; k++) {
                    // Parameters
                    int lind_df = ((k-roi[0].cirroi[num_thread].y)+roi[0].cirroi[num_thread].radius)*6+i*(roi[0].cirroi[num_thread].region.nodelist.height*6);

                    // Hessian - only calculate lower half since hessian is symmetric
                    hessian_gn_buffer[num_thread][0] += df_dp_buffer[num_thread][lind_df]*df_dp_buffer[num_thread][lind_df];
                    hessian_gn_buffer[num_thread][1] += df_dp_buffer[num_thread][lind_df]*df_dp_buffer[num_thread][lind_df+1];
                    hessian_gn_buffer[num_thread][2] += df_dp_buffer[num_thread][lind_df]*df_dp_buffer[num_thread][lind_df+2];
                    hessian_gn_buffer[num_thread][3] += df_dp_buffer[num_thread][lind_df]*df_dp_buffer[num_thread][lind_df+3];
                    hessian_gn_buffer[num_thread][4] += df_dp_buffer[num_thread][lind_df]*df_dp_buffer[num_thread][lind_df+4];
                    hessian_gn_buffer[num_thread][5] += df_dp_buffer[num_thread][lind_df]*df_dp_buffer[num_thread][lind_df+5];

                    hessian_gn_buffer[num_thread][7] += df_dp_buffer[num_thread][lind_df+1]*df_dp_buffer[num_thread][lind_df+1];
                    hessian_gn_buffer[num_thread][8] += df_dp_buffer[num_thread][lind_df+1]*df_dp_buffer[num_thread][lind_df+2];
                    hessian_gn_buffer[num_thread][9] += df_dp_buffer[num_thread][lind_df+1]*df_dp_buffer[num_thread][lind_df+3];
                    hessian_gn_buffer[num_thread][10] += df_dp_buffer[num_thread][lind_df+1]*df_dp_buffer[num_thread][lind_df+4];
                    hessian_gn_buffer[num_thread][11] += df_dp_buffer[num_thread][lind_df+1]*df_dp_buffer[num_thread][lind_df+5];

                    hessian_gn_buffer[num_thread][14] += df_dp_buffer[num_thread][lind_df+2]*df_dp_buffer[num_thread][lind_df+2];
                    hessian_gn_buffer[num_thread][15] += df_dp_buffer[num_thread][lind_df+2]*df_dp_buffer[num_thread][lind_df+3];
                    hessian_gn_buffer[num_thread][16] += df_dp_buffer[num_thread][lind_df+2]*df_dp_buffer[num_thread][lind_df+4];
                    hessian_gn_buffer[num_thread][17] += df_dp_buffer[num_thread][lind_df+2]*df_dp_buffer[num_thread][lind_df+5];

                    hessian_gn_buffer[num_thread][21] += df_dp_buffer[num_thread][lind_df+3]*df_dp_buffer[num_thread][lind_df+3];
                    hessian_gn_buffer[num_thread][22] += df_dp_buffer[num_thread][lind_df+3]*df_dp_buffer[num_thread][lind_df+4];
                    hessian_gn_buffer[num_thread][23] += df_dp_buffer[num_thread][lind_df+3]*df_dp_buffer[num_thread][lind_df+5];

                    hessian_gn_buffer[num_thread][28] += df_dp_buffer[num_thread][lind_df+4]*df_dp_buffer[num_thread][lind_df+4];
                    hessian_gn_buffer[num_thread][29] += df_dp_buffer[num_thread][lind_df+4]*df_dp_buffer[num_thread][lind_df+5];

                    hessian_gn_buffer[num_thread][35] += df_dp_buffer[num_thread][lind_df+5]*df_dp_buffer[num_thread][lind_df+5];    
                }
            }    
        }
                        
        // Multiply components of hessian by 2/deltaf^2
        for (int i=0; i<6; i++) {
            for (int j=i; j<6; j++) {
                hessian_gn_buffer[num_thread][j+i*6] *= 2*pow(deltaf_inv,2);
            }
        }

        // Fill other half of hessian
        hessian_gn_buffer[num_thread][6] = hessian_gn_buffer[num_thread][1];    

        hessian_gn_buffer[num_thread][12] = hessian_gn_buffer[num_thread][2];    
        hessian_gn_buffer[num_thread][13] = hessian_gn_buffer[num_thread][8];    
        
        hessian_gn_buffer[num_thread][18] = hessian_gn_buffer[num_thread][3];    
        hessian_gn_buffer[num_thread][19] = hessian_gn_buffer[num_thread][9];    
        hessian_gn_buffer[num_thread][20] = hessian_gn_buffer[num_thread][15];    

        hessian_gn_buffer[num_thread][24] = hessian_gn_buffer[num_thread][4];    
        hessian_gn_buffer[num_thread][25] = hessian_gn_buffer[num_thread][10];    
        hessian_gn_buffer[num_thread][26] = hessian_gn_buffer[num_thread][16];    
        hessian_gn_buffer[num_thread][27] = hessian_gn_buffer[num_thread][22];    
        
        hessian_gn_buffer[num_thread][30] = hessian_gn_buffer[num_thread][5];    
        hessian_gn_buffer[num_thread][31] = hessian_gn_buffer[num_thread][11];    
        hessian_gn_buffer[num_thread][32] = hessian_gn_buffer[num_thread][17];    
        hessian_gn_buffer[num_thread][33] = hessian_gn_buffer[num_thread][23];    
        hessian_gn_buffer[num_thread][34] = hessian_gn_buffer[num_thread][29];        
        
        // Solve for new parameters via cholesky decomp (from Golub and Van Loan)
        // Lower triangle of Hessian overwritten with parameters used in Cholesky decomp
        // If one of the diagonals is close to zero or negative, then the 
        // hessian is not positive definite
        bool positivedef = true;
        cholesky(hessian_gn_buffer[num_thread],positivedef,6);
                
        if (positivedef) {               
            // Start iterations - For first iteration use defvector_init
            double diffnorm;
            OUT outstate_newton = newton(defvector,corrcoef,diffnorm,defvector_init,fm,deltaf_inv,num_thread);     
            
            // Initialize counter
            int counter = 1; 
            while (outstate_newton == SUCCESS && diffnorm >= cutoff_diffnorm && counter <= cutoff_iteration) {
                // For rest of iterations use defvector from previous iterations
                outstate_newton = newton(defvector,corrcoef,diffnorm,defvector,fm,deltaf_inv,num_thread);     
                ++counter;
            }        
            
            if (outstate_newton == SUCCESS) {                                    
                return SUCCESS;
            }
        } 
    } 
    // Some parameters are invalid - either deltag_inv was zero or the hessian wasn't positive definite
    return FAILED;
}

OUT class_rgdic::newton(std::vector<double> &defvector,double &corrcoef,double &diffnorm,const std::vector<double> &defvector_init,const double &fm,const double &deltaf_inv,const int &num_thread) {
    // Will only overwrite queue_new if parameters are valid
    // Interpolate g subset - do this here instead of interp_qbs because QK_B_QKT has been precomputed
    double gm = 0.0;
    for (int i=0; i<roi[0].cirroi[num_thread].region.noderange.height; i++) {
        for (int j=0; j<roi[0].cirroi[num_thread].region.noderange.value[i]; j+=2) {
            for (int k=roi[0].cirroi[num_thread].region.nodelist.value[i+j*roi[0].cirroi[num_thread].region.nodelist.height]; k<=roi[0].cirroi[num_thread].region.nodelist.value[i+(j+1)*roi[0].cirroi[num_thread].region.nodelist.height]; k++) {
                // Find new coordinates
                double dx = (double)(i-roi[0].cirroi[num_thread].radius);
                double dy = (double)(k-roi[0].cirroi[num_thread].y);

                double y_tilda = (double)k + defvector_init[1] + defvector_init[4]*dx + defvector_init[5]*dy;
                double x_tilda = (double)(i+(roi[0].cirroi[num_thread].x-roi[0].cirroi[num_thread].radius)) + defvector_init[0] + defvector_init[2]*dx + defvector_init[3]*dy;

                int x_tilda_floor = (int)floor(x_tilda);
                int y_tilda_floor = (int)floor(y_tilda);

                int lind_g = (int)dy+roi[0].cirroi[num_thread].radius+i*roi[0].cirroi[num_thread].region.nodelist.height;
                
                // Get bounds of the desired b-spline coefficients used for interpolation
                int top = y_tilda_floor+current[0].border_bcoef-2;
                int left = x_tilda_floor+current[0].border_bcoef-2;
                int bottom = y_tilda_floor+current[0].border_bcoef+3;
                int right = x_tilda_floor+current[0].border_bcoef+3;
                
                if (top >= 0 &&
                    left >= 0 && 
                    bottom < current[0].bcoef.height &&
                    right < current[0].bcoef.width) {                    
                    double x_tilda_delta = x_tilda-(double)x_tilda_floor;
                    double y_tilda_delta = y_tilda-(double)y_tilda_floor;

                    // Form x_vec
                    x_vec_buffer[num_thread][1] = x_tilda_delta;
                    x_vec_buffer[num_thread][2] = x_tilda_delta*x_tilda_delta;
                    x_vec_buffer[num_thread][3] = x_tilda_delta*x_tilda_delta*x_tilda_delta;
                    x_vec_buffer[num_thread][4] = x_tilda_delta*x_tilda_delta*x_tilda_delta*x_tilda_delta;
                    x_vec_buffer[num_thread][5] = x_tilda_delta*x_tilda_delta*x_tilda_delta*x_tilda_delta*x_tilda_delta;

                    // Form y_vec
                    y_vec_buffer[num_thread][1] = y_tilda_delta;
                    y_vec_buffer[num_thread][2] = y_tilda_delta*y_tilda_delta;
                    y_vec_buffer[num_thread][3] = y_tilda_delta*y_tilda_delta*y_tilda_delta;
                    y_vec_buffer[num_thread][4] = y_tilda_delta*y_tilda_delta*y_tilda_delta*y_tilda_delta;
                    y_vec_buffer[num_thread][5] = y_tilda_delta*y_tilda_delta*y_tilda_delta*y_tilda_delta*y_tilda_delta;
                                    
                    // Calculate lind_qkbqkt for QK_B_QKT
                    int lind_qkbqkt = (top*36)+(left*36*(current[0].bcoef.height-5));

                    // Get QK_B_QKT coefficients
                    double QK_B_QKT_0 = QK_B_QKT_buffer[lind_qkbqkt];
                    double QK_B_QKT_1 = QK_B_QKT_buffer[lind_qkbqkt+1];
                    double QK_B_QKT_2 = QK_B_QKT_buffer[lind_qkbqkt+2];
                    double QK_B_QKT_3 = QK_B_QKT_buffer[lind_qkbqkt+3];
                    double QK_B_QKT_4 = QK_B_QKT_buffer[lind_qkbqkt+4];
                    double QK_B_QKT_5 = QK_B_QKT_buffer[lind_qkbqkt+5];
                    double QK_B_QKT_6 = QK_B_QKT_buffer[lind_qkbqkt+6];
                    double QK_B_QKT_7 = QK_B_QKT_buffer[lind_qkbqkt+7];
                    double QK_B_QKT_8 = QK_B_QKT_buffer[lind_qkbqkt+8];
                    double QK_B_QKT_9 = QK_B_QKT_buffer[lind_qkbqkt+9];
                    double QK_B_QKT_10 = QK_B_QKT_buffer[lind_qkbqkt+10];
                    double QK_B_QKT_11 = QK_B_QKT_buffer[lind_qkbqkt+11];
                    double QK_B_QKT_12 = QK_B_QKT_buffer[lind_qkbqkt+12];
                    double QK_B_QKT_13 = QK_B_QKT_buffer[lind_qkbqkt+13];
                    double QK_B_QKT_14 = QK_B_QKT_buffer[lind_qkbqkt+14];
                    double QK_B_QKT_15 = QK_B_QKT_buffer[lind_qkbqkt+15];
                    double QK_B_QKT_16 = QK_B_QKT_buffer[lind_qkbqkt+16];
                    double QK_B_QKT_17 = QK_B_QKT_buffer[lind_qkbqkt+17];
                    double QK_B_QKT_18 = QK_B_QKT_buffer[lind_qkbqkt+18];
                    double QK_B_QKT_19 = QK_B_QKT_buffer[lind_qkbqkt+19];
                    double QK_B_QKT_20 = QK_B_QKT_buffer[lind_qkbqkt+20];
                    double QK_B_QKT_21 = QK_B_QKT_buffer[lind_qkbqkt+21];
                    double QK_B_QKT_22 = QK_B_QKT_buffer[lind_qkbqkt+22];
                    double QK_B_QKT_23 = QK_B_QKT_buffer[lind_qkbqkt+23];
                    double QK_B_QKT_24 = QK_B_QKT_buffer[lind_qkbqkt+24];
                    double QK_B_QKT_25 = QK_B_QKT_buffer[lind_qkbqkt+25];
                    double QK_B_QKT_26 = QK_B_QKT_buffer[lind_qkbqkt+26];
                    double QK_B_QKT_27 = QK_B_QKT_buffer[lind_qkbqkt+27];
                    double QK_B_QKT_28 = QK_B_QKT_buffer[lind_qkbqkt+28];
                    double QK_B_QKT_29 = QK_B_QKT_buffer[lind_qkbqkt+29];
                    double QK_B_QKT_30 = QK_B_QKT_buffer[lind_qkbqkt+30];
                    double QK_B_QKT_31 = QK_B_QKT_buffer[lind_qkbqkt+31];
                    double QK_B_QKT_32 = QK_B_QKT_buffer[lind_qkbqkt+32];
                    double QK_B_QKT_33 = QK_B_QKT_buffer[lind_qkbqkt+33];
                    double QK_B_QKT_34 = QK_B_QKT_buffer[lind_qkbqkt+34];
                    double QK_B_QKT_35 = QK_B_QKT_buffer[lind_qkbqkt+35];

                    // Calculate g - main computational bottleneck of the inverse compositional method with biquintic b-splines
                    g_buffer[num_thread][lind_g] = (QK_B_QKT_0+x_vec_buffer[num_thread][1]*QK_B_QKT_6+x_vec_buffer[num_thread][2]*QK_B_QKT_12+x_vec_buffer[num_thread][3]*QK_B_QKT_18+x_vec_buffer[num_thread][4]*QK_B_QKT_24+x_vec_buffer[num_thread][5]*QK_B_QKT_30)+
                                                   (QK_B_QKT_1+x_vec_buffer[num_thread][1]*QK_B_QKT_7+x_vec_buffer[num_thread][2]*QK_B_QKT_13+x_vec_buffer[num_thread][3]*QK_B_QKT_19+x_vec_buffer[num_thread][4]*QK_B_QKT_25+x_vec_buffer[num_thread][5]*QK_B_QKT_31)*y_vec_buffer[num_thread][1]+
                                                   (QK_B_QKT_2+x_vec_buffer[num_thread][1]*QK_B_QKT_8+x_vec_buffer[num_thread][2]*QK_B_QKT_14+x_vec_buffer[num_thread][3]*QK_B_QKT_20+x_vec_buffer[num_thread][4]*QK_B_QKT_26+x_vec_buffer[num_thread][5]*QK_B_QKT_32)*y_vec_buffer[num_thread][2]+
                                                   (QK_B_QKT_3+x_vec_buffer[num_thread][1]*QK_B_QKT_9+x_vec_buffer[num_thread][2]*QK_B_QKT_15+x_vec_buffer[num_thread][3]*QK_B_QKT_21+x_vec_buffer[num_thread][4]*QK_B_QKT_27+x_vec_buffer[num_thread][5]*QK_B_QKT_33)*y_vec_buffer[num_thread][3]+
                                                   (QK_B_QKT_4+x_vec_buffer[num_thread][1]*QK_B_QKT_10+x_vec_buffer[num_thread][2]*QK_B_QKT_16+x_vec_buffer[num_thread][3]*QK_B_QKT_22+x_vec_buffer[num_thread][4]*QK_B_QKT_28+x_vec_buffer[num_thread][5]*QK_B_QKT_34)*y_vec_buffer[num_thread][4]+
                                                   (QK_B_QKT_5+x_vec_buffer[num_thread][1]*QK_B_QKT_11+x_vec_buffer[num_thread][2]*QK_B_QKT_17+x_vec_buffer[num_thread][3]*QK_B_QKT_23+x_vec_buffer[num_thread][4]*QK_B_QKT_29+x_vec_buffer[num_thread][5]*QK_B_QKT_35)*y_vec_buffer[num_thread][5];
                    
                    // Add components to calculate the mean
                    gm += g_buffer[num_thread][lind_g];
                } else {                    
                    // If this condition is satisfied then we are 
                    // interpolating a point beyond the bounds of the 
                    // original image, so just set the values to zero
                    g_buffer[num_thread][lind_g] = 0.0;
                    
                    // Don't add anything to averages
                    continue;
                }
            }
        }        
    }
    // Divide by totalpoints to get real average
    gm /= (double)roi[0].cirroi[num_thread].region.totalpoints;
    
    // Calculate deltag_inv
    double deltag_inv = 0.0;
    for (int i=0; i<roi[0].cirroi[num_thread].region.noderange.height; i++) {
        for (int j=0; j<roi[0].cirroi[num_thread].region.noderange.value[i]; j+=2) {
            for (int k=roi[0].cirroi[num_thread].region.nodelist.value[i+j*roi[0].cirroi[num_thread].region.nodelist.height]; k<=roi[0].cirroi[num_thread].region.nodelist.value[i+(j+1)*roi[0].cirroi[num_thread].region.nodelist.height]; k++) {
                int lind_g = (k-roi[0].cirroi[num_thread].y)+roi[0].cirroi[num_thread].radius+i*roi[0].cirroi[num_thread].region.nodelist.height;
                deltag_inv = deltag_inv+pow(g_buffer[num_thread][lind_g]-gm,2);
            }
        }
    }    
    deltag_inv = sqrt(deltag_inv); // This is deltag; will take inverse after ensuring it is not close to zero

    // check to make sure deltag_inv (strictly positive) isn't close to zero; if it is, exit newton raphson
    if (deltag_inv > LAMBDA) { 
        // Finish deltag_inv
        deltag_inv = 1.0/deltag_inv;        
        
        // Calculate gradient
        // Initialize to zero first
        std::fill(gradient_buffer[num_thread].begin(),gradient_buffer[num_thread].end(),0.0);
        corrcoef = 0.0;
        for (int i=0; i<roi[0].cirroi[num_thread].region.noderange.height; i++) {
            for (int j=0; j < roi[0].cirroi[num_thread].region.noderange.value[i]; j+=2) {
                for (int k=roi[0].cirroi[num_thread].region.nodelist.value[i+j*roi[0].cirroi[num_thread].region.nodelist.height]; k<=roi[0].cirroi[num_thread].region.nodelist.value[i+(j+1)*roi[0].cirroi[num_thread].region.nodelist.height]; k++) {
                    // Parameters
                    int lind_f = k+(i+(roi[0].cirroi[num_thread].x-roi[0].cirroi[num_thread].radius))*reference[0].gs.height;
                    int lind_df = ((k-roi[0].cirroi[num_thread].y)+roi[0].cirroi[num_thread].radius)*6+i*(roi[0].cirroi[num_thread].region.nodelist.height*6);
                    int lind_g = (k-roi[0].cirroi[num_thread].y)+roi[0].cirroi[num_thread].radius+i*roi[0].cirroi[num_thread].region.nodelist.height;
                    
                    // Gradient Parameters
                    double normalized_diff = (reference[0].gs.value[lind_f]-fm)*deltaf_inv-(g_buffer[num_thread][lind_g]-gm)*deltag_inv;
                                       
                    // Gradient
                    gradient_buffer[num_thread][0] += normalized_diff*df_dp_buffer[num_thread][lind_df];
                    gradient_buffer[num_thread][1] += normalized_diff*df_dp_buffer[num_thread][lind_df+1];
                    gradient_buffer[num_thread][2] += normalized_diff*df_dp_buffer[num_thread][lind_df+2];
                    gradient_buffer[num_thread][3] += normalized_diff*df_dp_buffer[num_thread][lind_df+3];
                    gradient_buffer[num_thread][4] += normalized_diff*df_dp_buffer[num_thread][lind_df+4];
                    gradient_buffer[num_thread][5] += normalized_diff*df_dp_buffer[num_thread][lind_df+5]; 
                    
                    // Correlation coefficient
                    corrcoef += pow(normalized_diff,2);
                }
            }    
        }
        
        // Update gradient; multiply by the inverses
        for (int i=0; i<6; i++) {
            gradient_buffer[num_thread][i] *= 2*deltaf_inv;
        }        

        // Find new change in deformation parameters        
        // Ax = b
        // GG'x = b, where G is lower triangular
        // Gy = b -> G'x = y
        // Step 1: solve for y with forward substitution; y is stored in gradient_buffer
        forwardsub(gradient_buffer[num_thread],hessian_gn_buffer[num_thread],6);

        // Step 2: solve for x with back substitution
        backwardsub(gradient_buffer[num_thread],hessian_gn_buffer[num_thread],6);
        
        // Make gradient_buffer negative
        for (int i=0; i<6; i++) { 
            gradient_buffer[num_thread][i] = -gradient_buffer[num_thread][i];
        }

        // At this point the change in deformation parameters is stored in gradient_buffer    
        // Calculate difference norm - this is stored in gradient_buffer at this point
        diffnorm = 0.0;
        for (int i=0; i<6; i++) { 
            diffnorm += gradient_buffer[num_thread][i]*gradient_buffer[num_thread][i];
        }
        diffnorm = sqrt(diffnorm);  

        // Update parameters using inverse composition
        // Transfer parameters because defvector_init is an alias of defvector after the first iteration
        double defvector_init_u = defvector_init[0];
        double defvector_init_v = defvector_init[1];
        double defvector_init_dudx = defvector_init[2];
        double defvector_init_dudy = defvector_init[3];
        double defvector_init_dvdx = defvector_init[4];
        double defvector_init_dvdy = defvector_init[5];
        defvector[0] = defvector_init_u - ((defvector_init_dudx + 1)*(gradient_buffer[num_thread][0] + gradient_buffer[num_thread][0]*gradient_buffer[num_thread][5] - gradient_buffer[num_thread][1]*gradient_buffer[num_thread][3]))/(gradient_buffer[num_thread][2] + gradient_buffer[num_thread][5] + gradient_buffer[num_thread][2]*gradient_buffer[num_thread][5] - gradient_buffer[num_thread][3]*gradient_buffer[num_thread][4] + 1) - (defvector_init_dudy*(gradient_buffer[num_thread][1] - gradient_buffer[num_thread][0]*gradient_buffer[num_thread][4] + gradient_buffer[num_thread][1]*gradient_buffer[num_thread][2]))/(gradient_buffer[num_thread][2] + gradient_buffer[num_thread][5] + gradient_buffer[num_thread][2]*gradient_buffer[num_thread][5] - gradient_buffer[num_thread][3]*gradient_buffer[num_thread][4] + 1); // u
        defvector[1] = defvector_init_v - ((defvector_init_dvdy + 1)*(gradient_buffer[num_thread][1] - gradient_buffer[num_thread][0]*gradient_buffer[num_thread][4] + gradient_buffer[num_thread][1]*gradient_buffer[num_thread][2]))/(gradient_buffer[num_thread][2] + gradient_buffer[num_thread][5] + gradient_buffer[num_thread][2]*gradient_buffer[num_thread][5] - gradient_buffer[num_thread][3]*gradient_buffer[num_thread][4] + 1) - (defvector_init_dvdx*(gradient_buffer[num_thread][0] + gradient_buffer[num_thread][0]*gradient_buffer[num_thread][5] - gradient_buffer[num_thread][1]*gradient_buffer[num_thread][3]))/(gradient_buffer[num_thread][2] + gradient_buffer[num_thread][5] + gradient_buffer[num_thread][2]*gradient_buffer[num_thread][5] - gradient_buffer[num_thread][3]*gradient_buffer[num_thread][4] + 1); // v
        defvector[2] = ((gradient_buffer[num_thread][5] + 1)*(defvector_init_dudx + 1))/(gradient_buffer[num_thread][2] + gradient_buffer[num_thread][5] + gradient_buffer[num_thread][2]*gradient_buffer[num_thread][5] - gradient_buffer[num_thread][3]*gradient_buffer[num_thread][4] + 1) - (gradient_buffer[num_thread][4]*defvector_init_dudy)/(gradient_buffer[num_thread][2] + gradient_buffer[num_thread][5] + gradient_buffer[num_thread][2]*gradient_buffer[num_thread][5] - gradient_buffer[num_thread][3]*gradient_buffer[num_thread][4] + 1) - 1; // du/dx
        defvector[3] = (defvector_init_dudy*(gradient_buffer[num_thread][2] + 1))/(gradient_buffer[num_thread][2] + gradient_buffer[num_thread][5] + gradient_buffer[num_thread][2]*gradient_buffer[num_thread][5] - gradient_buffer[num_thread][3]*gradient_buffer[num_thread][4] + 1) - (gradient_buffer[num_thread][3]*(defvector_init_dudx + 1))/(gradient_buffer[num_thread][2] + gradient_buffer[num_thread][5] + gradient_buffer[num_thread][2]*gradient_buffer[num_thread][5] - gradient_buffer[num_thread][3]*gradient_buffer[num_thread][4] + 1); // du/dy
        defvector[4] = (defvector_init_dvdx*(gradient_buffer[num_thread][5] + 1))/(gradient_buffer[num_thread][2] + gradient_buffer[num_thread][5] + gradient_buffer[num_thread][2]*gradient_buffer[num_thread][5] - gradient_buffer[num_thread][3]*gradient_buffer[num_thread][4] + 1) - (gradient_buffer[num_thread][4]*(defvector_init_dvdy + 1))/(gradient_buffer[num_thread][2] + gradient_buffer[num_thread][5] + gradient_buffer[num_thread][2]*gradient_buffer[num_thread][5] - gradient_buffer[num_thread][3]*gradient_buffer[num_thread][4] + 1); // dv/dx
        defvector[5] = ((gradient_buffer[num_thread][2] + 1)*(defvector_init_dvdy + 1))/(gradient_buffer[num_thread][2] + gradient_buffer[num_thread][5] + gradient_buffer[num_thread][2]*gradient_buffer[num_thread][5] - gradient_buffer[num_thread][3]*gradient_buffer[num_thread][4] + 1) - (gradient_buffer[num_thread][3]*defvector_init_dvdx)/(gradient_buffer[num_thread][2] + gradient_buffer[num_thread][5] + gradient_buffer[num_thread][2]*gradient_buffer[num_thread][5] - gradient_buffer[num_thread][3]*gradient_buffer[num_thread][4] + 1) - 1; // dv/dy

        // Return successful
        return SUCCESS;   
    }
    // Deltag was close to zero - return failed
    return FAILED;
}

void mexFunction(int nlhs,mxArray *plhs[ ],int nrhs,const mxArray *prhs[ ]) {
    if (nrhs == 12 && nlhs == 2) {
        // Create rgdic
        class_rgdic rgdic(plhs,prhs);
        
        // Run analysis
        rgdic.analysis();
    } else {
        // Thread safe because it is single threaded up to this point
        mexErrMsgTxt("Incorrect number of inputs or outputs.\n");
    }
}
