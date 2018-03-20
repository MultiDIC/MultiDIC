// This function calculates the seeds for the RG-DIC routine for a specified region and current image. 
// If no flags are passed, then a simple single threaded version will be compiled
// Since this is multithreaded, make sure all functions called in the OpenMP region are threadsafe!

#include <mex.h>
#include <math.h>
#include <vector>
#include <algorithm>
#include <list>
#include <exception>              // Allow exceptions because this function can allocate a very large chunk of memory for the interpolation lookup table
#include <complex>
#include "standard_datatypes.h"
#include "ncorr_datatypes.h"
#include "ncorr_lib.h"
            
#ifdef NCORR_OPENMP
    #include <omp.h>              // openmp header
#endif

// ----------------------------------------------------//
// Main Class -----------------------------------------//
// ----------------------------------------------------//

class class_calcseeds {
public:
    // Constructor
    class_calcseeds(mxArray *plhs [ ],const mxArray *prhs [ ]);
    
    // Methods
    void analysis();

private:
    // Properties
    // Inputs:
    std::vector<ncorr_class_img> reference;       // ncorr datatype
    std::vector<ncorr_class_img> current;         // ncorr datatype
    std::vector<ncorr_class_roi> roi;             // ncorr datatype
    int num_region;                               // standard datatype
    class_integer_array pos_seed;                 // standard datatype
    int radius;                                   // standard datatype
    double cutoff_diffnorm;                       // standard datatype
    int cutoff_iteration;                         // standard datatype
    bool enabled_stepanalysis;                    // standard datatype
    bool subsettrunc;                             // standard datatype

    // Outputs:
    std::vector<class_double_array> paramvector_si;
    std::vector<class_double_array> num_region_si;
    std::vector<class_double_array> num_thread_si;
    std::vector<class_double_array> computepoints_si;
    std::vector<class_double_array> num_iterations_c;
    std::vector<class_double_array> diffnorm_c;
    double *outstate;
    
    // Other variables:
    std::vector<int> vec_outstate;                 
    int total_threads;
    
    // Inverse Compositional Buffers: Must form these explicitly
    std::vector<std::vector<double> > g_buffer; 
    std::vector<std::vector<double> > df_dp_buffer;    
    std::vector<std::vector<double> > x_vec_buffer;
    std::vector<std::vector<double> > y_vec_buffer; 
    std::vector<std::vector<double> > gradient_buffer;
    std::vector<std::vector<double> > hessian_gn_buffer; 
    std::vector<std::vector<double> > QK_B_QKT_buffer;

    // Methods    
    OUT calcpoint(std::vector<double> &paramvector,double &diffnorm,int &num_iterations,const int &x,const int &y,const int &num_thread);
    OUT initialguess(std::vector<double> &defvector,const int &num_thread);
    OUT ncc(std::vector<int> &disp_ncc,const int &reduction_multigrid,const std::vector<int> &disp_prev,const int &num_thread);
    OUT iterativesearch(std::vector<double> &defvector,double &corrcoef,double &diffnorm,int &num_iterations,const std::vector<double> &defvector_init,const int &num_thread);
    OUT newton(std::vector<double> &defvector,double &corrcoef,double &diffnorm,const std::vector<double> &defvector_init,const double &fm,const double &deltaf_inv,const int &num_thread);
};

class_calcseeds::class_calcseeds(mxArray *plhs[ ],const mxArray *prhs[ ]) {
    // Get inputs --------------------------------------------------// 
    // input 1: reference image    
    get_imgs(reference,prhs[0]);
    // input 2: current image    
    get_imgs(current,prhs[1]);
    // input 3: ROI
    get_rois(roi,prhs[2]);
    // input 4: region number
    get_integer_scalar(num_region,prhs[3]);    
    // input 5: pos_seed
    get_integer_array(pos_seed,prhs[4]);
    // input 6: radius
    get_integer_scalar(radius,prhs[5]);    
    // input 7: diffnorm cutoff
    get_double_scalar(cutoff_diffnorm,prhs[6]);
    // input 8: iteration cutoff
    get_integer_scalar(cutoff_iteration,prhs[7]);
    // input 9: enabled_stepanalysis
    get_logical_scalar(enabled_stepanalysis,prhs[8]);    
    // input 10: subsettrunc
    get_logical_scalar(subsettrunc,prhs[9]);    
        
    // Check inputs - These are very basic checks
    if (reference[0].gs.width == roi[0].mask.width && 
        reference[0].gs.height == roi[0].mask.height &&
        radius > 0 &&
        cutoff_diffnorm > 0.0 &&
        cutoff_iteration > 0) {
        // Get total number of threads from seed_pos --------------------//
        total_threads = (int)mxGetM(prhs[4]);
        
        // OpenMP Setup -------------------------------------------------//
        #ifdef NCORR_OPENMP
            // Set number of threads
            omp_set_num_threads(total_threads);
        #endif
                
        // Set cirroi ---------------------------------------------------//
        // Must set one for each thread
        roi[0].set_cirroi(radius,total_threads);
           
        // Success Vector -----------------------------------------------//
        vec_outstate.resize(total_threads,0); // One for each thread

        // Form Inverse Compositional Buffers ---------------------------//  
        g_buffer.resize(total_threads); 
        df_dp_buffer.resize(total_threads); 
        x_vec_buffer.resize(total_threads); 
        y_vec_buffer.resize(total_threads); 
        gradient_buffer.resize(total_threads); 
        hessian_gn_buffer.resize(total_threads); 
        QK_B_QKT_buffer.resize(total_threads);
        for (int i=0; i<total_threads; i++) {
            g_buffer[i].resize((radius*2+1)*(radius*2+1),0.0); 
            df_dp_buffer[i].resize((radius*2+1)*(radius*2+1)*6,0.0); 
            x_vec_buffer[i].resize(6,0.0); 
            y_vec_buffer[i].resize(6,0.0); 
            gradient_buffer[i].resize(6,0.0); 
            hessian_gn_buffer[i].resize(36,0.0); 
            QK_B_QKT_buffer[i].resize(36,0.0);
        }

        // Form/set outputs ---------------------------------------------//
        // output 1: seedinfo
        mwSize dims1[2] = {static_cast<mwSize>(total_threads),1};
        int numfields1 = 4;
        const char *fieldnames1[] = {"paramvector","num_region","num_thread","computepoints"};
        mxArray *mat_seedinfo = mxCreateStructArray(2, dims1, numfields1, fieldnames1);      
        plhs[0] = mat_seedinfo;
                
        std::vector<mxArray*> mat_paramvector_si(total_threads,NULL);
        std::vector<mxArray*> mat_num_region_si(total_threads,NULL);
        std::vector<mxArray*> mat_num_thread_si(total_threads,NULL);
        std::vector<mxArray*> mat_computepoints_si(total_threads,NULL);
        for (int i=0; i<total_threads; i++) {
            // Form fields
            mat_paramvector_si[i] = mxCreateDoubleMatrix(1, 9, mxREAL); 
            mat_num_region_si[i] = mxCreateDoubleMatrix(1, 1, mxREAL); 
            mat_num_thread_si[i] = mxCreateDoubleMatrix(1, 1, mxREAL); 
            mat_computepoints_si[i] = mxCreateDoubleMatrix(1, 1, mxREAL); 

            // Add fields to structure
            // add paramvector:
            mxSetFieldByNumber(mat_seedinfo,i,0,mat_paramvector_si[i]);
            // add num_region:
            mxSetFieldByNumber(mat_seedinfo,i,1,mat_num_region_si[i]);
            // add num_thread:
            mxSetFieldByNumber(mat_seedinfo,i,2,mat_num_thread_si[i]);
            // add computepoints:
            mxSetFieldByNumber(mat_seedinfo,i,3,mat_computepoints_si[i]);
        }

        // output 2: convergence
        mwSize dims2[2] = {static_cast<mwSize>(total_threads),1};
        int numfields2 = 2;
        const char *fieldnames2[] = {"num_iterations","diffnorm"};
        mxArray *mat_convergence = mxCreateStructArray(2, dims2, numfields2, fieldnames2);      
        plhs[1] = mat_convergence;

        std::vector<mxArray*> mat_num_iterations_c(total_threads,NULL);
        std::vector<mxArray*> mat_diffnorm_c(total_threads,NULL);
        for (int i=0; i<total_threads; i++) {
            // Form fields
            mat_num_iterations_c[i] = mxCreateDoubleMatrix(1, 1, mxREAL); 
            mat_diffnorm_c[i] = mxCreateDoubleMatrix(1, 1, mxREAL); 

            // Add fields to structure
            // add num_iterations:
            mxSetFieldByNumber(mat_convergence,i,0,mat_num_iterations_c[i]);
            // add diffnorm:
            mxSetFieldByNumber(mat_convergence,i,1,mat_diffnorm_c[i]);
        }

        // output 3: outstate
        plhs[2] = mxCreateDoubleMatrix(1,1,mxREAL);
        
        // Get outputs --------------------------------------------------//
        // output 1: seedinfo
        paramvector_si.resize(total_threads);
        num_region_si.resize(total_threads);
        num_thread_si.resize(total_threads);
        computepoints_si.resize(total_threads);
        for (int i=0; i<total_threads; i++) {
            get_double_array(paramvector_si[i],mat_paramvector_si[i]);
            get_double_array(num_region_si[i],mat_num_region_si[i]);
            get_double_array(num_thread_si[i],mat_num_thread_si[i]);
            get_double_array(computepoints_si[i],mat_computepoints_si[i]);
        }
    
        // output 2: convergence
        num_iterations_c.resize(total_threads);
        diffnorm_c.resize(total_threads);
        for (int i=0; i<total_threads; i++) {
            get_double_array(num_iterations_c[i],mat_num_iterations_c[i]);
            get_double_array(diffnorm_c[i],mat_diffnorm_c[i]);
        }
    
        // output 3: outstate
        outstate = mxGetPr(plhs[2]);
    } else {
        // Thread safe because it is single threaded up to this point
        mexErrMsgTxt("One of the inputs is incorrect.\n");
    }
}

// ----------------------------------------------------//
// Main Class Methods ---------------------------------//
// ----------------------------------------------------//

void class_calcseeds::analysis() {     
    // Initialize outstate to success - note that this function isnt cancellable; it can only fail
    *outstate = (double)SUCCESS; 
            
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
                                
        // Update cirroi ---------------------------------------//
        // Need to update cirroi - This is threadsafe
        roi[0].update_cirroi(num_region,num_thread);

        // Form buffers
        std::vector<double> paramvector(9,0);     // [x y u v du/dx du/dy dv/dx dv/dy corrcoef]
        double diffnorm;
        int num_iterations;

        // Analyze each point ----------------------------------//
        int x = pos_seed.value[num_thread];
        int y = pos_seed.value[num_thread+pos_seed.height];
        vec_outstate[num_thread] = calcpoint(paramvector,diffnorm,num_iterations,x,y,num_thread);

        // Store outputs ---------------------------------------//
        // Paramvector = [x y u v du/dx du/dy dv/dx dv/dy corrcoef]
        paramvector_si[num_thread].value[0] = paramvector[0];
        paramvector_si[num_thread].value[1] = paramvector[1];
        paramvector_si[num_thread].value[2] = paramvector[2];
        paramvector_si[num_thread].value[3] = paramvector[3];
        paramvector_si[num_thread].value[4] = paramvector[4];
        paramvector_si[num_thread].value[5] = paramvector[5];
        paramvector_si[num_thread].value[6] = paramvector[6];
        paramvector_si[num_thread].value[7] = paramvector[7];
        paramvector_si[num_thread].value[8] = paramvector[8];
        num_region_si[num_thread].value[0] = num_region;
        num_thread_si[num_thread].value[0] = num_thread;
        // Note that the number of compute points are assigned later

        num_iterations_c[num_thread].value[0] = num_iterations;
        diffnorm_c[num_thread].value[0] = diffnorm;

    #ifdef NCORR_OPENMP
    }    
    #endif
       
    // Check vec_success to make sure all seeds processed correctly
    for (int i=0; i<total_threads; i++) {
        if (!vec_outstate[i]) {
            *outstate = (double)FAILED;
        }
    }
}

OUT class_calcseeds::calcpoint(std::vector<double> &paramvector,double &diffnorm,int &num_iterations,const int &x,const int &y,const int &num_thread) {
    // Form cirroi -> Find initial guess with NCC -> Refine results with IC-GN -> Store outputs and return true or false
    // Step 1: Get cirroi
    roi[0].get_cirroi(x,y,num_region,subsettrunc,num_thread);    

    // Step 2: Get initial guess
    std::vector<double> defvector_init(6,0);
    OUT outstate_initialguess = initialguess(defvector_init,num_thread);

    if (outstate_initialguess == SUCCESS) {        
        // Step 3: Get refined results with IC-GN
        std::vector<double> defvector(6,0);
        double corrcoef;
        OUT outstate_iterative = iterativesearch(defvector,corrcoef,diffnorm,num_iterations,defvector_init,num_thread);   

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
    }    
    return FAILED;
}

OUT class_calcseeds::initialguess(std::vector<double> &defvector,const int &num_thread) {   
    // This function returns the initial guess. It uses multigrid normalized cross correlation
    // Define nccvec_multigrid - this defines the steps multigrid will take.
    // For now only do two steps. If the first step fails then generally the
    // rest will fail; the converse, where if the first succeeds then 
    // generally the rest will succeed. Thus, the first step is very
    // important and only two steps are used.
    std::vector<int> nccvec_multigrid(1,0);
    if (floor((double)roi[0].cirroi[num_thread].radius/20.0) != 0) {
        nccvec_multigrid.resize(2,0);
        nccvec_multigrid[0] = (int)floor((double)roi[0].cirroi[num_thread].radius/20.0);
    }

    // Perform multigrid normalized circular cross correlation, first
    // iteration uses whole (reduced) image so pass an empty vector as last argument,
    // this argument will generally be the coordinates of the center of the 
    // subset with the highest correlation coefficient from the previous
    // iteration.
    std::vector<int> disp_ncc(2,0);
    std::vector<int> disp_prev;        // Having an empty disp_prev lets NCC know to use the whole current image
    OUT outstate_ncc = ncc(disp_ncc,nccvec_multigrid[0],disp_prev,num_thread);    

    if (outstate_ncc == SUCCESS) {
        // Resize disp_prev - This lets NCC know to use a truncated portion
        disp_prev.resize(2,0);
        // Now iteration over reduction factors        
        for (int i=1; i<(int)nccvec_multigrid.size(); i++) {
            // Transfer answer to previdisp
            disp_prev[0] = disp_ncc[0];
            disp_prev[1] = disp_ncc[1];
            outstate_ncc = ncc(disp_ncc,nccvec_multigrid[i],disp_prev,num_thread); 
            // Check to see if NCC is successful for each iteration
            if (outstate_ncc == FAILED) {
                break;
            }
        }

        if (outstate_ncc == SUCCESS) {
            // Assign outputs
            defvector[0] = (double)disp_ncc[0];
            defvector[1] = (double)disp_ncc[1];
            return SUCCESS;
        }
    }
    return FAILED;
}

OUT class_calcseeds::ncc(std::vector<int> &disp_ncc,const int &reduction_multigrid,const std::vector<int> &disp_prev,const int &num_thread) {   
    // Carries out multiscale normalized cross correlation
    // Get reduced current region -> Get reduced reference subset -> perform NCC
    // This routine actually gets the reduced current region and reduced reference subset
    // and stores them, so add exception handling for the allocation of the
    // current region since this can potentially be large if many threads are used. 
    // This routine re-forms the arrays to avoid strided access which can be slow.
    
    // Note that the up_*, down_*, left_*, and right_* prefixes are meant to
    // specify the bounding box for a region.

    // Step 1: Get reduced current "region"
    std::vector<double> cur_reduced;
    int height_cur_reduced = 0;
    int width_cur_reduced = 0;

    // For the first iteration, use the reduced whole current image for the current region
    // For later iterations, use a reduced portion using the guess from the previous iteration.
    int up_cur_multigrid = 0;   // Distance from top side of image to first grid point
    int left_cur_multigrid = 0; // Distance from left side of image to first grid point
    if (disp_prev.size() == 0) {
        // This is the first iteration, use entire (reduced) current image.
        // These coordinates depend on the reduction_multigrid, coordinates
        // of the subset, and the size of the image.

        // Form reduced current image
        height_cur_reduced = (int)ceil(((double)current[0].gs.height)/((double)reduction_multigrid+1.0));
        width_cur_reduced = (int)ceil(((double)current[0].gs.width)/((double)reduction_multigrid+1.0));
                  
        try { 
            cur_reduced.resize(height_cur_reduced*width_cur_reduced); // Could potentially be large if many threads are uses
        } catch (std::exception&) {
            // Just return failed
            return FAILED;
        }
        
        for (int i=0; i<current[0].gs.width; i+=reduction_multigrid+1) {            
            for (int j=0; j<current[0].gs.height; j+=reduction_multigrid+1) {
                cur_reduced[j/(reduction_multigrid+1)+(i/(reduction_multigrid+1))*height_cur_reduced] = current[0].gs.value[j+i*current[0].gs.height];
            }
        }
    } else {
        // This is the next iteration, use area near region found near
        // estimated subset location from previous iteration. 

        // Get reduced images, use coordinates from before
        int x_subset = roi[0].cirroi[num_thread].x + disp_prev[0];
        int y_subset = roi[0].cirroi[num_thread].y + disp_prev[1];

        // Truncation factor - use successively smaller windows for smaller
        // multigrid reductions. This is because smaller grid reductions
        // result in higher resolution data which is slower to compute. 
        double truncfactor = reduction_multigrid + 1.5;
            
        // Determine bounds and offset - use min and max to prevent using
        // bounds outside the current image bounds. 
        int up_cur = (int)ncorr_round(std::max((double)y_subset-truncfactor*(double)roi[0].cirroi[num_thread].radius,0.0));
        int down_cur = (int)ncorr_round(std::min((double)y_subset+truncfactor*(double)roi[0].cirroi[num_thread].radius,(double)current[0].gs.height-1.0));        
        int left_cur = (int)ncorr_round(std::max((double)x_subset-truncfactor*(double)roi[0].cirroi[num_thread].radius,0.0));
        int right_cur = (int)ncorr_round(std::min((double)x_subset+truncfactor*(double)roi[0].cirroi[num_thread].radius,(double)current[0].gs.width-1.0));
        
        // Set bounds to grid - DO NOT REDFINE up* and left* for cur_multigrid
        up_cur_multigrid = y_subset-(int)floor(((double)y_subset-(double)up_cur)/((double)reduction_multigrid+1.0))*(reduction_multigrid+1);
        int down_cur_multigrid = y_subset+(int)floor(((double)down_cur-(double)y_subset)/((double)reduction_multigrid+1.0))*(reduction_multigrid+1);
        left_cur_multigrid = x_subset-(int)floor(((double)x_subset-(double)left_cur)/((double)reduction_multigrid+1.0))*(reduction_multigrid+1);
        int right_cur_multigrid = x_subset+(int)floor(((double)right_cur-(double)x_subset)/((double)reduction_multigrid+1.0))*(reduction_multigrid+1);
             
        // Get reduced current image
        height_cur_reduced = (int)ceil(((double)down_cur_multigrid-up_cur_multigrid+1)/((double)reduction_multigrid+1.0));
        width_cur_reduced = (int)ceil(((double)right_cur_multigrid-left_cur_multigrid+1)/((double)reduction_multigrid+1.0));
        cur_reduced.resize(height_cur_reduced*width_cur_reduced);
        for (int i=left_cur_multigrid; i<=right_cur_multigrid; i+=reduction_multigrid+1) {            
            for (int j=up_cur_multigrid; j<=down_cur_multigrid; j+=reduction_multigrid+1) {
                cur_reduced[(j-up_cur_multigrid)/(reduction_multigrid+1)+((i-left_cur_multigrid)/(reduction_multigrid+1))*height_cur_reduced] = current[0].gs.value[j+i*current[0].gs.height];
            }
        }
    }

    // Step 2: Get reduced reference subset
    std::vector<double> ref_reduced;
    int height_ref_reduced = 0;
    int width_ref_reduced = 0;

    // NOTE: For a region, the leftbound and rightbound represent where the
    // region starts and stops, so they differ a little bit from the upper
    // and lower bounds. To get the bounding box of a cirroi, you can use
    // the upper and lower bounds for up_* and down_*, but for left_* and
    // right_* you must find them directly
    int up_cirroi = roi[0].cirroi[num_thread].region.upperbound;                          // Copy directly
    int down_cirroi = roi[0].cirroi[num_thread].region.lowerbound;                        // Copy directly
    int left_cirroi = roi[0].cirroi[num_thread].x-roi[0].cirroi[num_thread].radius;       // Initialize
    int right_cirroi = roi[0].cirroi[num_thread].x-roi[0].cirroi[num_thread].radius;      // Initialize

    bool firstpoint = false;
    for (int i=0; i<=2*roi[0].cirroi[num_thread].radius; i++) {
        int x = i+roi[0].cirroi[num_thread].x-roi[0].cirroi[num_thread].radius;
        if (roi[0].cirroi[num_thread].region.noderange.value[i] > 0 && !firstpoint) {
            // This is the left bound
            left_cirroi = x;
            firstpoint = true;
        }
        
        if (roi[0].cirroi[num_thread].region.noderange.value[i] > 0 && firstpoint && x > right_cirroi) {
            right_cirroi = x;        
        }
    }

    // Get multigrid bounding box for cirroi
    int up_cirroi_multigrid = roi[0].cirroi[num_thread].y - (int)floor(((double)roi[0].cirroi[num_thread].y-(double)up_cirroi)/((double)reduction_multigrid+1.0))*(reduction_multigrid+1);
    int down_cirroi_multigrid = roi[0].cirroi[num_thread].y + (int)floor(((double)down_cirroi-(double)roi[0].cirroi[num_thread].y)/((double)reduction_multigrid+1.0))*(reduction_multigrid+1);
    int left_cirroi_multigrid = roi[0].cirroi[num_thread].x - (int)floor(((double)roi[0].cirroi[num_thread].x-(double)left_cirroi)/((double)reduction_multigrid+1.0))*(reduction_multigrid+1);
    int right_cirroi_multigrid = roi[0].cirroi[num_thread].x + (int)floor(((double)right_cirroi-(double)roi[0].cirroi[num_thread].x)/((double)reduction_multigrid+1.0))*(reduction_multigrid+1);
    
    // Get space between first grid point and top and left side of the reference image
    int up_ref_multigrid = roi[0].cirroi[num_thread].y-(int)floor((double)roi[0].cirroi[num_thread].y/((double)reduction_multigrid+1.0))*(reduction_multigrid+1);    
    int left_ref_multigrid = roi[0].cirroi[num_thread].x-(int)floor((double)roi[0].cirroi[num_thread].x/((double)reduction_multigrid+1.0))*(reduction_multigrid+1);
       
    // Get bounding box of cirroi for reduced ref subset
    int up_cirroi_reduced_ref = (up_cirroi_multigrid-up_ref_multigrid)/(reduction_multigrid+1);
    int down_cirroi_reduced_ref = (down_cirroi_multigrid-up_ref_multigrid)/(reduction_multigrid+1);
    int left_cirroi_reduced_ref = (left_cirroi_multigrid-left_ref_multigrid)/(reduction_multigrid+1);
    int right_cirroi_reduced_ref = (right_cirroi_multigrid-left_ref_multigrid)/(reduction_multigrid+1);

    // Resize reduced reference subset - this is usually small so don't check for out of memory exception
    height_ref_reduced = down_cirroi_reduced_ref-up_cirroi_reduced_ref+1;
    width_ref_reduced = right_cirroi_reduced_ref-left_cirroi_reduced_ref+1;
    ref_reduced.resize(height_ref_reduced*width_ref_reduced);

    // Now copy elements to finish reduced reference subset
    // Also get info for total elements, fm, etc
    // Also form cirroi_mask_reduced for ease of debugging and implementation
    std::vector<char> cirroi_mask_reduced(height_ref_reduced*width_ref_reduced,0);
    int totalpoints_reduced = 0;
    double fm = 0.0;
    for (int i=left_cirroi_multigrid; i<=right_cirroi_multigrid; i+=reduction_multigrid+1) {
        int x_ref_reduced = (i-left_cirroi_multigrid)/(reduction_multigrid+1);
        int x_tl = i-roi[0].cirroi[num_thread].x+roi[0].cirroi[num_thread].radius;
        for (int j=up_cirroi_multigrid; j<=down_cirroi_multigrid; j+=reduction_multigrid+1) {
            int y_ref_reduced = (j-up_cirroi_multigrid)/(reduction_multigrid+1);
            int y_tl = j-roi[0].cirroi[num_thread].y+roi[0].cirroi[num_thread].radius;
            if (roi[0].cirroi[num_thread].mask.value[y_tl+x_tl*roi[0].cirroi[num_thread].mask.height]) {
                ref_reduced[y_ref_reduced+x_ref_reduced*height_ref_reduced] = reference[0].gs.value[j+i*reference[0].gs.height];
                cirroi_mask_reduced[y_ref_reduced+x_ref_reduced*height_ref_reduced] = 1;
                fm += ref_reduced[y_ref_reduced+x_ref_reduced*height_ref_reduced];
                ++totalpoints_reduced;
            }
        }
    }
    // Finish fm
    fm /= (double)totalpoints_reduced;

    // Now get fsquaredsum and also subtract mean from reference subwindow
    double fsquaredsum = 0;
    for (int i=0; i<width_ref_reduced; i++) {
        for (int j=0; j<height_ref_reduced; j++) {
            if (cirroi_mask_reduced[j+i*height_ref_reduced]) {
                // Subtract mean
                ref_reduced[j+i*height_ref_reduced] -= fm;
                // Sum the squared value
                fsquaredsum += pow(ref_reduced[j+i*height_ref_reduced],2); // Mean is subtracted already
            }
        }
    }

    // Now we have the mean subtracted reference subset and the fsquared sum
    // in eq.5 precomputed. Now compute eq.5 for each point in the reduced
    // current region. Do this in the spatial domain since an FFT call would
    // require a MATLAB call which isn't safe - OR - would require the FFTW
    // library which would add a lot of files.
    int height_nccmatrix = height_cur_reduced-height_ref_reduced+1; 
    int width_nccmatrix = width_cur_reduced-width_ref_reduced+1;     
    
    double max_val_cc = -1;    // Note range is from [-1 1]
    int idx_max_cc_x = -1;     // Use this for error checking
    int idx_max_cc_y = -1;     // Use this for error checking
    for (int i=0; i<width_nccmatrix; i++) {
        for (int j=0; j<height_nccmatrix; j++) {
            // Get gm first
            double gm = 0;
            for (int k=0; k<width_ref_reduced; k++) {
                int x_cur_reduced = i+k;
                for (int l=0; l<height_ref_reduced; l++) {
                    int y_cur_reduced = j+l;
                    if (cirroi_mask_reduced[l+k*height_ref_reduced]) {
                        gm += cur_reduced[y_cur_reduced+x_cur_reduced*height_cur_reduced];
                    }
                }
            }
            // Finish gm 
            gm /= (double)totalpoints_reduced;

            // Get numerator of eq.5 and gsquaredsum
            double numerator = 0;
            double gsquaredsum = 0;
            for (int k=0; k<width_ref_reduced; k++) {
                int x_cur_reduced = i+k;
                for (int l=0; l<height_ref_reduced; l++) {
                    int y_cur_reduced = j+l;
                    if (cirroi_mask_reduced[l+k*height_ref_reduced]) {
                        // Remember ref_reduced already has mean subtracted
                        numerator += ref_reduced[l+k*height_ref_reduced]*(cur_reduced[y_cur_reduced+x_cur_reduced*height_cur_reduced]-gm);
                        gsquaredsum += pow(cur_reduced[y_cur_reduced+x_cur_reduced*height_cur_reduced]-gm,2);
                    }
                }
            }

            // Get denominator and check to make sure its not close to zero
            double denominator = sqrt(fsquaredsum*gsquaredsum);
            if (denominator > LAMBDA) {
                // Get the NCC value
                double val_cc = numerator/denominator;
                
                if (val_cc > max_val_cc) {
                    max_val_cc = val_cc;
                    idx_max_cc_x = i;
                    idx_max_cc_y = j;
                }
            }
        }
    }

    // Debug
    // if (num_thread == 0) {
    //     imshow(&cur_reduced[0],width_cur_reduced,height_cur_reduced);
    //     imshow(&ref_reduced[0],width_ref_reduced,height_ref_reduced);
    // }

    // Check to make sure indices arent -1
    if (idx_max_cc_x != -1 && idx_max_cc_y != -1) {
        // Determine displacements - get reduced space from top and left of cirroi
        int y_space_reduced = (roi[0].cirroi[num_thread].y-up_cirroi_multigrid)/(reduction_multigrid+1);                    // Distance from the topleft to the center
        int x_space_reduced = (roi[0].cirroi[num_thread].x-left_cirroi_multigrid)/(reduction_multigrid+1);                  // Distance from the topleft to the center
    
        int u = ((idx_max_cc_x+x_space_reduced)*(reduction_multigrid+1)-roi[0].cirroi[num_thread].x)+left_cur_multigrid;    // Coordinates are pixels
        int v = ((idx_max_cc_y+y_space_reduced)*(reduction_multigrid+1)-roi[0].cirroi[num_thread].y)+up_cur_multigrid;      // Coordinates are pixels
        
        // Store
        disp_ncc[0] = u;
        disp_ncc[1] = v;

        return SUCCESS;
    }
    return FAILED;
}

OUT class_calcseeds::iterativesearch(std::vector<double> &defvector,double &corrcoef,double &diffnorm,int &num_iterations,const std::vector<double> &defvector_init, const int &num_thread) {   
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
                    
                    // Get bounds of the desired b-spline coefficients used for interpolation
                    // These will be in bounds if border_bcoef is correctly greater than or
                    // equal to 2.
                    int top = y_tilda_floor+reference[0].border_bcoef-2;
                    int left = x_tilda_floor+reference[0].border_bcoef-2;
                    int bottom = y_tilda_floor+reference[0].border_bcoef+3;
                    int right = x_tilda_floor+reference[0].border_bcoef+3;

                    // Get corresponding b-spline coefficients
                    double b0 = reference[0].bcoef.value[(top)+(left)*reference[0].bcoef.height];
                    double b1 = reference[0].bcoef.value[(top+1)+(left)*reference[0].bcoef.height];
                    double b2 = reference[0].bcoef.value[(top+2)+(left)*reference[0].bcoef.height];
                    double b3 = reference[0].bcoef.value[(top+3)+(left)*reference[0].bcoef.height];
                    double b4 = reference[0].bcoef.value[(top+4)+(left)*reference[0].bcoef.height];
                    double b5 = reference[0].bcoef.value[(top+5)+(left)*reference[0].bcoef.height];            
                    double b6 = reference[0].bcoef.value[(top)+(left+1)*reference[0].bcoef.height];
                    double b7 = reference[0].bcoef.value[(top+1)+(left+1)*reference[0].bcoef.height];
                    double b8 = reference[0].bcoef.value[(top+2)+(left+1)*reference[0].bcoef.height];
                    double b9 = reference[0].bcoef.value[(top+3)+(left+1)*reference[0].bcoef.height];
                    double b10 = reference[0].bcoef.value[(top+4)+(left+1)*reference[0].bcoef.height];
                    double b11 = reference[0].bcoef.value[(top+5)+(left+1)*reference[0].bcoef.height];            
                    double b12 = reference[0].bcoef.value[(top)+(left+2)*reference[0].bcoef.height];
                    double b13 = reference[0].bcoef.value[(top+1)+(left+2)*reference[0].bcoef.height];
                    double b14 = reference[0].bcoef.value[(top+2)+(left+2)*reference[0].bcoef.height];            
                    double b15 = reference[0].bcoef.value[(top+3)+(left+2)*reference[0].bcoef.height];
                    double b16 = reference[0].bcoef.value[(top+4)+(left+2)*reference[0].bcoef.height];
                    double b17 = reference[0].bcoef.value[(top+5)+(left+2)*reference[0].bcoef.height];            
                    double b18 = reference[0].bcoef.value[(top)+(left+3)*reference[0].bcoef.height];
                    double b19 = reference[0].bcoef.value[(top+1)+(left+3)*reference[0].bcoef.height];
                    double b20 = reference[0].bcoef.value[(top+2)+(left+3)*reference[0].bcoef.height];
                    double b21 = reference[0].bcoef.value[(top+3)+(left+3)*reference[0].bcoef.height];
                    double b22 = reference[0].bcoef.value[(top+4)+(left+3)*reference[0].bcoef.height];
                    double b23 = reference[0].bcoef.value[(top+5)+(left+3)*reference[0].bcoef.height];            
                    double b24 = reference[0].bcoef.value[(top)+(left+4)*reference[0].bcoef.height];
                    double b25 = reference[0].bcoef.value[(top+1)+(left+4)*reference[0].bcoef.height];
                    double b26 = reference[0].bcoef.value[(top+2)+(left+4)*reference[0].bcoef.height];
                    double b27 = reference[0].bcoef.value[(top+3)+(left+4)*reference[0].bcoef.height];
                    double b28 = reference[0].bcoef.value[(top+4)+(left+4)*reference[0].bcoef.height];
                    double b29 = reference[0].bcoef.value[(top+5)+(left+4)*reference[0].bcoef.height];
                    double b30 = reference[0].bcoef.value[(top)+(left+5)*reference[0].bcoef.height];
                    double b31 = reference[0].bcoef.value[(top+1)+(left+5)*reference[0].bcoef.height];
                    double b32 = reference[0].bcoef.value[(top+2)+(left+5)*reference[0].bcoef.height];
                    double b33 = reference[0].bcoef.value[(top+3)+(left+5)*reference[0].bcoef.height];
                    double b34 = reference[0].bcoef.value[(top+4)+(left+5)*reference[0].bcoef.height];
                    double b35 = reference[0].bcoef.value[(top+5)+(left+5)*reference[0].bcoef.height];
                                        
                    // Compute Gradients using b-spline coefficients
                    // First order
                    double df_dx_buffer = 0.003472222222222222*b18-0.009027777777777778*b1-0.003472222222222222*b10-0.0003472222222222222*b0+0.09027777777777778*b19-0.02291666666666667*b2+0.2291666666666667*b20+0.09027777777777778*b21+0.003472222222222222*b22+0.0003472222222222222*b24+0.009027777777777778*b25+0.02291666666666667*b26+0.009027777777777778*b27+0.0003472222222222222*b28-0.009027777777777778*b3-0.0003472222222222222*b4-0.003472222222222222*b6-0.09027777777777778*b7-0.2291666666666667*b8-0.09027777777777778*b9;
                    double df_dy_buffer = 0.009027777777777778*b10-0.003472222222222222*b1-0.0003472222222222222*b0-0.02291666666666667*b12-0.2291666666666667*b13+0.2291666666666667*b15+0.02291666666666667*b16-0.009027777777777778*b18-0.09027777777777778*b19+0.09027777777777778*b21+0.009027777777777778*b22-0.0003472222222222222*b24-0.003472222222222222*b25+0.003472222222222222*b27+0.0003472222222222222*b28+0.003472222222222222*b3+0.0003472222222222222*b4-0.009027777777777778*b6-0.09027777777777778*b7+0.09027777777777778*b9;

                    // First order
                    int lind_df = ((k-roi[0].cirroi[num_thread].y)+roi[0].cirroi[num_thread].radius)*6+i*(roi[0].cirroi[num_thread].region.nodelist.height*6);
                    df_dp_buffer[num_thread][lind_df] = df_dx_buffer;      // u
                    df_dp_buffer[num_thread][lind_df+1] = df_dy_buffer;    // v
                    df_dp_buffer[num_thread][lind_df+2] = df_dx_buffer*dx; // dudx
                    df_dp_buffer[num_thread][lind_df+3] = df_dx_buffer*dy; // dudy
                    df_dp_buffer[num_thread][lind_df+4] = df_dy_buffer*dx; // dvdx
                    df_dp_buffer[num_thread][lind_df+5] = df_dy_buffer*dy; // dvdy
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
            OUT outstate_newton = newton(defvector,corrcoef,diffnorm,defvector_init,fm,deltaf_inv,num_thread);    
            
            // Initialize counter
            int counter = 1; 
            while (outstate_newton == SUCCESS && diffnorm > cutoff_diffnorm && counter < cutoff_iteration) {
                // For rest of iterations use defvector from previous iterations
                outstate_newton = newton(defvector,corrcoef,diffnorm,defvector,fm,deltaf_inv,num_thread);    
                ++counter;
            }      

            // Store number of iterations
            num_iterations = counter;
            
            if (outstate_newton == SUCCESS) {                                    
                return SUCCESS;
            }
        } 
    } 
    // Some parameters are invalid - either deltag_inv was zero or the hessian wasn't positive definite
    return FAILED;
}

OUT class_calcseeds::newton(std::vector<double> &defvector,double &corrcoef,double &diffnorm,const std::vector<double> &defvector_init,const double &fm,const double &deltaf_inv,const int &num_thread) {
    // Will only overwrite queue_new if parameters are valid
    // Interpolate g subset - Do this here instead of interp_qbs because QK_B_QKT has been precomputed
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
                                 
                    // Get corresponding b-spline coefficients
                    double b0 = current[0].bcoef.value[(top)+(left)*current[0].bcoef.height];
                    double b1 = current[0].bcoef.value[(top+1)+(left)*current[0].bcoef.height];
                    double b2 = current[0].bcoef.value[(top+2)+(left)*current[0].bcoef.height];
                    double b3 = current[0].bcoef.value[(top+3)+(left)*current[0].bcoef.height];
                    double b4 = current[0].bcoef.value[(top+4)+(left)*current[0].bcoef.height];
                    double b5 = current[0].bcoef.value[(top+5)+(left)*current[0].bcoef.height];            
                    double b6 = current[0].bcoef.value[(top)+(left+1)*current[0].bcoef.height];
                    double b7 = current[0].bcoef.value[(top+1)+(left+1)*current[0].bcoef.height];
                    double b8 = current[0].bcoef.value[(top+2)+(left+1)*current[0].bcoef.height];
                    double b9 = current[0].bcoef.value[(top+3)+(left+1)*current[0].bcoef.height];
                    double b10 = current[0].bcoef.value[(top+4)+(left+1)*current[0].bcoef.height];
                    double b11 = current[0].bcoef.value[(top+5)+(left+1)*current[0].bcoef.height];            
                    double b12 = current[0].bcoef.value[(top)+(left+2)*current[0].bcoef.height];
                    double b13 = current[0].bcoef.value[(top+1)+(left+2)*current[0].bcoef.height];
                    double b14 = current[0].bcoef.value[(top+2)+(left+2)*current[0].bcoef.height];            
                    double b15 = current[0].bcoef.value[(top+3)+(left+2)*current[0].bcoef.height];
                    double b16 = current[0].bcoef.value[(top+4)+(left+2)*current[0].bcoef.height];
                    double b17 = current[0].bcoef.value[(top+5)+(left+2)*current[0].bcoef.height];            
                    double b18 = current[0].bcoef.value[(top)+(left+3)*current[0].bcoef.height];
                    double b19 = current[0].bcoef.value[(top+1)+(left+3)*current[0].bcoef.height];
                    double b20 = current[0].bcoef.value[(top+2)+(left+3)*current[0].bcoef.height];
                    double b21 = current[0].bcoef.value[(top+3)+(left+3)*current[0].bcoef.height];
                    double b22 = current[0].bcoef.value[(top+4)+(left+3)*current[0].bcoef.height];
                    double b23 = current[0].bcoef.value[(top+5)+(left+3)*current[0].bcoef.height];            
                    double b24 = current[0].bcoef.value[(top)+(left+4)*current[0].bcoef.height];
                    double b25 = current[0].bcoef.value[(top+1)+(left+4)*current[0].bcoef.height];
                    double b26 = current[0].bcoef.value[(top+2)+(left+4)*current[0].bcoef.height];
                    double b27 = current[0].bcoef.value[(top+3)+(left+4)*current[0].bcoef.height];
                    double b28 = current[0].bcoef.value[(top+4)+(left+4)*current[0].bcoef.height];
                    double b29 = current[0].bcoef.value[(top+5)+(left+4)*current[0].bcoef.height];
                    double b30 = current[0].bcoef.value[(top)+(left+5)*current[0].bcoef.height];
                    double b31 = current[0].bcoef.value[(top+1)+(left+5)*current[0].bcoef.height];
                    double b32 = current[0].bcoef.value[(top+2)+(left+5)*current[0].bcoef.height];
                    double b33 = current[0].bcoef.value[(top+3)+(left+5)*current[0].bcoef.height];
                    double b34 = current[0].bcoef.value[(top+4)+(left+5)*current[0].bcoef.height];
                    double b35 = current[0].bcoef.value[(top+5)+(left+5)*current[0].bcoef.height];

                    // Compute QK*B*QK^T
                    QK_B_QKT_buffer[num_thread][0] = 0.00006944444444444444*b0+0.001805555555555556*b1+0.001805555555555556*b10+0.004583333333333333*b12+0.1191666666666667*b13+0.3025*b14+0.1191666666666667*b15+0.004583333333333333*b16+0.001805555555555556*b18+0.04694444444444444*b19+0.004583333333333333*b2+0.1191666666666667*b20+0.04694444444444444*b21+0.001805555555555556*b22+0.00006944444444444444*b24+0.001805555555555556*b25+0.004583333333333333*b26+0.001805555555555556*b27+0.00006944444444444444*b28+0.001805555555555556*b3+0.00006944444444444444*b4+0.001805555555555556*b6+0.04694444444444444*b7+0.1191666666666667*b8+0.04694444444444444*b9;
                    QK_B_QKT_buffer[num_thread][1] = 0.009027777777777778*b10-0.003472222222222222*b1-0.0003472222222222222*b0-0.02291666666666667*b12-0.2291666666666667*b13+0.2291666666666667*b15+0.02291666666666667*b16-0.009027777777777778*b18-0.09027777777777778*b19+0.09027777777777778*b21+0.009027777777777778*b22-0.0003472222222222222*b24-0.003472222222222222*b25+0.003472222222222222*b27+0.0003472222222222222*b28+0.003472222222222222*b3+0.0003472222222222222*b4-0.009027777777777778*b6-0.09027777777777778*b7+0.09027777777777778*b9;
                    QK_B_QKT_buffer[num_thread][2] = 0.0006944444444444444*b0+0.001388888888888889*b1+0.01805555555555556*b10+0.04583333333333333*b12+0.09166666666666667*b13-0.275*b14+0.09166666666666667*b15+0.04583333333333333*b16+0.01805555555555556*b18+0.03611111111111111*b19-0.004166666666666667*b2-0.1083333333333333*b20+0.03611111111111111*b21+0.01805555555555556*b22+0.0006944444444444444*b24+0.001388888888888889*b25-0.004166666666666667*b26+0.001388888888888889*b27+0.0006944444444444444*b28+0.001388888888888889*b3+0.0006944444444444444*b4+0.01805555555555556*b6+0.03611111111111111*b7-0.1083333333333333*b8+0.03611111111111111*b9;
                    QK_B_QKT_buffer[num_thread][3] = 0.001388888888888889*b1-0.0006944444444444444*b0+0.01805555555555556*b10-0.04583333333333333*b12+0.09166666666666667*b13-0.09166666666666667*b15+0.04583333333333333*b16-0.01805555555555556*b18+0.03611111111111111*b19-0.03611111111111111*b21+0.01805555555555556*b22-0.0006944444444444444*b24+0.001388888888888889*b25-0.001388888888888889*b27+0.0006944444444444444*b28-0.001388888888888889*b3+0.0006944444444444444*b4-0.01805555555555556*b6+0.03611111111111111*b7-0.03611111111111111*b9;
                    QK_B_QKT_buffer[num_thread][4] = 0.0003472222222222222*b0-0.001388888888888889*b1+0.009027777777777778*b10+0.02291666666666667*b12-0.09166666666666667*b13+0.1375*b14-0.09166666666666667*b15+0.02291666666666667*b16+0.009027777777777778*b18-0.03611111111111111*b19+0.002083333333333333*b2+0.05416666666666667*b20-0.03611111111111111*b21+0.009027777777777778*b22+0.0003472222222222222*b24-0.001388888888888889*b25+0.002083333333333333*b26-0.001388888888888889*b27+0.0003472222222222222*b28-0.001388888888888889*b3+0.0003472222222222222*b4+0.009027777777777778*b6-0.03611111111111111*b7+0.05416666666666667*b8-0.03611111111111111*b9;
                    QK_B_QKT_buffer[num_thread][5] = 0.0003472222222222222*b1-0.00006944444444444444*b0-0.009027777777777778*b10+0.001805555555555556*b11-0.004583333333333333*b12+0.02291666666666667*b13-0.04583333333333333*b14+0.04583333333333333*b15-0.02291666666666667*b16+0.004583333333333333*b17-0.001805555555555556*b18+0.009027777777777778*b19-0.0006944444444444444*b2-0.01805555555555556*b20+0.01805555555555556*b21-0.009027777777777778*b22+0.001805555555555556*b23-0.00006944444444444444*b24+0.0003472222222222222*b25-0.0006944444444444444*b26+0.0006944444444444444*b27-0.0003472222222222222*b28+0.00006944444444444444*b29+0.0006944444444444444*b3-0.0003472222222222222*b4+0.00006944444444444444*b5-0.001805555555555556*b6+0.009027777777777778*b7-0.01805555555555556*b8+0.01805555555555556*b9;
                    QK_B_QKT_buffer[num_thread][6] = 0.003472222222222222*b18-0.009027777777777778*b1-0.003472222222222222*b10-0.0003472222222222222*b0+0.09027777777777778*b19-0.02291666666666667*b2+0.2291666666666667*b20+0.09027777777777778*b21+0.003472222222222222*b22+0.0003472222222222222*b24+0.009027777777777778*b25+0.02291666666666667*b26+0.009027777777777778*b27+0.0003472222222222222*b28-0.009027777777777778*b3-0.0003472222222222222*b4-0.003472222222222222*b6-0.09027777777777778*b7-0.2291666666666667*b8-0.09027777777777778*b9;
                    QK_B_QKT_buffer[num_thread][7] = 0.001736111111111111*b0+0.01736111111111111*b1-0.01736111111111111*b10-0.01736111111111111*b18-0.1736111111111111*b19+0.1736111111111111*b21+0.01736111111111111*b22-0.001736111111111111*b24-0.01736111111111111*b25+0.01736111111111111*b27+0.001736111111111111*b28-0.01736111111111111*b3-0.001736111111111111*b4+0.01736111111111111*b6+0.1736111111111111*b7-0.1736111111111111*b9;
                    QK_B_QKT_buffer[num_thread][8] = 0.03472222222222222*b18-0.006944444444444444*b1-0.03472222222222222*b10-0.003472222222222222*b0+0.06944444444444444*b19+0.02083333333333333*b2-0.2083333333333333*b20+0.06944444444444444*b21+0.03472222222222222*b22+0.003472222222222222*b24+0.006944444444444444*b25-0.02083333333333333*b26+0.006944444444444444*b27+0.003472222222222222*b28-0.006944444444444444*b3-0.003472222222222222*b4-0.03472222222222222*b6-0.06944444444444444*b7+0.2083333333333333*b8-0.06944444444444444*b9;
                    QK_B_QKT_buffer[num_thread][9] = 0.003472222222222222*b0-0.006944444444444444*b1-0.03472222222222222*b10-0.03472222222222222*b18+0.06944444444444444*b19-0.06944444444444444*b21+0.03472222222222222*b22-0.003472222222222222*b24+0.006944444444444444*b25-0.006944444444444444*b27+0.003472222222222222*b28+0.006944444444444444*b3-0.003472222222222222*b4+0.03472222222222222*b6-0.06944444444444444*b7+0.06944444444444444*b9;
                    QK_B_QKT_buffer[num_thread][10] = 0.006944444444444444*b1-0.001736111111111111*b0-0.01736111111111111*b10+0.01736111111111111*b18-0.06944444444444444*b19-0.01041666666666667*b2+0.1041666666666667*b20-0.06944444444444444*b21+0.01736111111111111*b22+0.001736111111111111*b24-0.006944444444444444*b25+0.01041666666666667*b26-0.006944444444444444*b27+0.001736111111111111*b28+0.006944444444444444*b3-0.001736111111111111*b4-0.01736111111111111*b6+0.06944444444444444*b7-0.1041666666666667*b8+0.06944444444444444*b9;
                    QK_B_QKT_buffer[num_thread][11] = 0.0003472222222222222*b0-0.001736111111111111*b1+0.01736111111111111*b10-0.003472222222222222*b11-0.003472222222222222*b18+0.01736111111111111*b19+0.003472222222222222*b2-0.03472222222222222*b20+0.03472222222222222*b21-0.01736111111111111*b22+0.003472222222222222*b23-0.0003472222222222222*b24+0.001736111111111111*b25-0.003472222222222222*b26+0.003472222222222222*b27-0.001736111111111111*b28+0.0003472222222222222*b29-0.003472222222222222*b3+0.001736111111111111*b4-0.0003472222222222222*b5+0.003472222222222222*b6-0.01736111111111111*b7+0.03472222222222222*b8-0.03472222222222222*b9;
                    QK_B_QKT_buffer[num_thread][12] = 0.0006944444444444444*b0+0.01805555555555556*b1+0.001388888888888889*b10-0.004166666666666667*b12-0.1083333333333333*b13-0.275*b14-0.1083333333333333*b15-0.004166666666666667*b16+0.001388888888888889*b18+0.03611111111111111*b19+0.04583333333333333*b2+0.09166666666666667*b20+0.03611111111111111*b21+0.001388888888888889*b22+0.0006944444444444444*b24+0.01805555555555556*b25+0.04583333333333333*b26+0.01805555555555556*b27+0.0006944444444444444*b28+0.01805555555555556*b3+0.0006944444444444444*b4+0.001388888888888889*b6+0.03611111111111111*b7+0.09166666666666667*b8+0.03611111111111111*b9;
                    QK_B_QKT_buffer[num_thread][13] = 0.006944444444444444*b10-0.03472222222222222*b1-0.003472222222222222*b0+0.02083333333333333*b12+0.2083333333333333*b13-0.2083333333333333*b15-0.02083333333333333*b16-0.006944444444444444*b18-0.06944444444444444*b19+0.06944444444444444*b21+0.006944444444444444*b22-0.003472222222222222*b24-0.03472222222222222*b25+0.03472222222222222*b27+0.003472222222222222*b28+0.03472222222222222*b3+0.003472222222222222*b4-0.006944444444444444*b6-0.06944444444444444*b7+0.06944444444444444*b9;
                    QK_B_QKT_buffer[num_thread][14] = 0.006944444444444444*b0+0.01388888888888889*b1+0.01388888888888889*b10-0.04166666666666667*b12-0.08333333333333333*b13+0.25*b14-0.08333333333333333*b15-0.04166666666666667*b16+0.01388888888888889*b18+0.02777777777777778*b19-0.04166666666666667*b2-0.08333333333333333*b20+0.02777777777777778*b21+0.01388888888888889*b22+0.006944444444444444*b24+0.01388888888888889*b25-0.04166666666666667*b26+0.01388888888888889*b27+0.006944444444444444*b28+0.01388888888888889*b3+0.006944444444444444*b4+0.01388888888888889*b6+0.02777777777777778*b7-0.08333333333333333*b8+0.02777777777777778*b9;
                    QK_B_QKT_buffer[num_thread][15] = 0.01388888888888889*b1-0.006944444444444444*b0+0.01388888888888889*b10+0.04166666666666667*b12-0.08333333333333333*b13+0.08333333333333333*b15-0.04166666666666667*b16-0.01388888888888889*b18+0.02777777777777778*b19-0.02777777777777778*b21+0.01388888888888889*b22-0.006944444444444444*b24+0.01388888888888889*b25-0.01388888888888889*b27+0.006944444444444444*b28-0.01388888888888889*b3+0.006944444444444444*b4-0.01388888888888889*b6+0.02777777777777778*b7-0.02777777777777778*b9;
                    QK_B_QKT_buffer[num_thread][16] = 0.003472222222222222*b0-0.01388888888888889*b1+0.006944444444444444*b10-0.02083333333333333*b12+0.08333333333333333*b13-0.125*b14+0.08333333333333333*b15-0.02083333333333333*b16+0.006944444444444444*b18-0.02777777777777778*b19+0.02083333333333333*b2+0.04166666666666667*b20-0.02777777777777778*b21+0.006944444444444444*b22+0.003472222222222222*b24-0.01388888888888889*b25+0.02083333333333333*b26-0.01388888888888889*b27+0.003472222222222222*b28-0.01388888888888889*b3+0.003472222222222222*b4+0.006944444444444444*b6-0.02777777777777778*b7+0.04166666666666667*b8-0.02777777777777778*b9;
                    QK_B_QKT_buffer[num_thread][17] = 0.003472222222222222*b1-0.0006944444444444444*b0-0.006944444444444444*b10+0.001388888888888889*b11+0.004166666666666667*b12-0.02083333333333333*b13+0.04166666666666667*b14-0.04166666666666667*b15+0.02083333333333333*b16-0.004166666666666667*b17-0.001388888888888889*b18+0.006944444444444444*b19-0.006944444444444444*b2-0.01388888888888889*b20+0.01388888888888889*b21-0.006944444444444444*b22+0.001388888888888889*b23-0.0006944444444444444*b24+0.003472222222222222*b25-0.006944444444444444*b26+0.006944444444444444*b27-0.003472222222222222*b28+0.0006944444444444444*b29+0.006944444444444444*b3-0.003472222222222222*b4+0.0006944444444444444*b5-0.001388888888888889*b6+0.006944444444444444*b7-0.01388888888888889*b8+0.01388888888888889*b9;
                    QK_B_QKT_buffer[num_thread][18] = 0.001388888888888889*b10-0.01805555555555556*b1-0.0006944444444444444*b0-0.001388888888888889*b18-0.03611111111111111*b19-0.04583333333333333*b2-0.09166666666666667*b20-0.03611111111111111*b21-0.001388888888888889*b22+0.0006944444444444444*b24+0.01805555555555556*b25+0.04583333333333333*b26+0.01805555555555556*b27+0.0006944444444444444*b28-0.01805555555555556*b3-0.0006944444444444444*b4+0.001388888888888889*b6+0.03611111111111111*b7+0.09166666666666667*b8+0.03611111111111111*b9;
                    QK_B_QKT_buffer[num_thread][19] = 0.003472222222222222*b0+0.03472222222222222*b1+0.006944444444444444*b10+0.006944444444444444*b18+0.06944444444444444*b19-0.06944444444444444*b21-0.006944444444444444*b22-0.003472222222222222*b24-0.03472222222222222*b25+0.03472222222222222*b27+0.003472222222222222*b28-0.03472222222222222*b3-0.003472222222222222*b4-0.006944444444444444*b6-0.06944444444444444*b7+0.06944444444444444*b9;
                    QK_B_QKT_buffer[num_thread][20] = 0.01388888888888889*b10-0.01388888888888889*b1-0.006944444444444444*b0-0.01388888888888889*b18-0.02777777777777778*b19+0.04166666666666667*b2+0.08333333333333333*b20-0.02777777777777778*b21-0.01388888888888889*b22+0.006944444444444444*b24+0.01388888888888889*b25-0.04166666666666667*b26+0.01388888888888889*b27+0.006944444444444444*b28-0.01388888888888889*b3-0.006944444444444444*b4+0.01388888888888889*b6+0.02777777777777778*b7-0.08333333333333333*b8+0.02777777777777778*b9;
                    QK_B_QKT_buffer[num_thread][21] = 0.006944444444444444*b0-0.01388888888888889*b1+0.01388888888888889*b10+0.01388888888888889*b18-0.02777777777777778*b19+0.02777777777777778*b21-0.01388888888888889*b22-0.006944444444444444*b24+0.01388888888888889*b25-0.01388888888888889*b27+0.006944444444444444*b28+0.01388888888888889*b3-0.006944444444444444*b4-0.01388888888888889*b6+0.02777777777777778*b7-0.02777777777777778*b9;
                    QK_B_QKT_buffer[num_thread][22] = 0.01388888888888889*b1-0.003472222222222222*b0+0.006944444444444444*b10-0.006944444444444444*b18+0.02777777777777778*b19-0.02083333333333333*b2-0.04166666666666667*b20+0.02777777777777778*b21-0.006944444444444444*b22+0.003472222222222222*b24-0.01388888888888889*b25+0.02083333333333333*b26-0.01388888888888889*b27+0.003472222222222222*b28+0.01388888888888889*b3-0.003472222222222222*b4+0.006944444444444444*b6-0.02777777777777778*b7+0.04166666666666667*b8-0.02777777777777778*b9;
                    QK_B_QKT_buffer[num_thread][23] = 0.0006944444444444444*b0-0.003472222222222222*b1-0.006944444444444444*b10+0.001388888888888889*b11+0.001388888888888889*b18-0.006944444444444444*b19+0.006944444444444444*b2+0.01388888888888889*b20-0.01388888888888889*b21+0.006944444444444444*b22-0.001388888888888889*b23-0.0006944444444444444*b24+0.003472222222222222*b25-0.006944444444444444*b26+0.006944444444444444*b27-0.003472222222222222*b28+0.0006944444444444444*b29-0.006944444444444444*b3+0.003472222222222222*b4-0.0006944444444444444*b5-0.001388888888888889*b6+0.006944444444444444*b7-0.01388888888888889*b8+0.01388888888888889*b9;
                    QK_B_QKT_buffer[num_thread][24] = 0.0003472222222222222*b0+0.009027777777777778*b1-0.001388888888888889*b10+0.002083333333333333*b12+0.05416666666666667*b13+0.1375*b14+0.05416666666666667*b15+0.002083333333333333*b16-0.001388888888888889*b18-0.03611111111111111*b19+0.02291666666666667*b2-0.09166666666666667*b20-0.03611111111111111*b21-0.001388888888888889*b22+0.0003472222222222222*b24+0.009027777777777778*b25+0.02291666666666667*b26+0.009027777777777778*b27+0.0003472222222222222*b28+0.009027777777777778*b3+0.0003472222222222222*b4-0.001388888888888889*b6-0.03611111111111111*b7-0.09166666666666667*b8-0.03611111111111111*b9;
                    QK_B_QKT_buffer[num_thread][25] = 0.1041666666666667*b15-0.01736111111111111*b1-0.006944444444444444*b10-0.01041666666666667*b12-0.1041666666666667*b13-0.001736111111111111*b0+0.01041666666666667*b16+0.006944444444444444*b18+0.06944444444444444*b19-0.06944444444444444*b21-0.006944444444444444*b22-0.001736111111111111*b24-0.01736111111111111*b25+0.01736111111111111*b27+0.001736111111111111*b28+0.01736111111111111*b3+0.001736111111111111*b4+0.006944444444444444*b6+0.06944444444444444*b7-0.06944444444444444*b9;
                    QK_B_QKT_buffer[num_thread][26] = 0.003472222222222222*b0+0.006944444444444444*b1-0.01388888888888889*b10+0.02083333333333333*b12+0.04166666666666667*b13-0.125*b14+0.04166666666666667*b15+0.02083333333333333*b16-0.01388888888888889*b18-0.02777777777777778*b19-0.02083333333333333*b2+0.08333333333333333*b20-0.02777777777777778*b21-0.01388888888888889*b22+0.003472222222222222*b24+0.006944444444444444*b25-0.02083333333333333*b26+0.006944444444444444*b27+0.003472222222222222*b28+0.006944444444444444*b3+0.003472222222222222*b4-0.01388888888888889*b6-0.02777777777777778*b7+0.08333333333333333*b8-0.02777777777777778*b9;
                    QK_B_QKT_buffer[num_thread][27] = 0.006944444444444444*b1-0.003472222222222222*b0-0.01388888888888889*b10-0.02083333333333333*b12+0.04166666666666667*b13-0.04166666666666667*b15+0.02083333333333333*b16+0.01388888888888889*b18-0.02777777777777778*b19+0.02777777777777778*b21-0.01388888888888889*b22-0.003472222222222222*b24+0.006944444444444444*b25-0.006944444444444444*b27+0.003472222222222222*b28-0.006944444444444444*b3+0.003472222222222222*b4+0.01388888888888889*b6-0.02777777777777778*b7+0.02777777777777778*b9;
                    QK_B_QKT_buffer[num_thread][28] = 0.001736111111111111*b0-0.006944444444444444*b1-0.006944444444444444*b10+0.01041666666666667*b12-0.04166666666666667*b13+0.0625*b14-0.04166666666666667*b15+0.01041666666666667*b16-0.006944444444444444*b18+0.02777777777777778*b19+0.01041666666666667*b2-0.04166666666666667*b20+0.02777777777777778*b21-0.006944444444444444*b22+0.001736111111111111*b24-0.006944444444444444*b25+0.01041666666666667*b26-0.006944444444444444*b27+0.001736111111111111*b28-0.006944444444444444*b3+0.001736111111111111*b4-0.006944444444444444*b6+0.02777777777777778*b7-0.04166666666666667*b8+0.02777777777777778*b9;
                    QK_B_QKT_buffer[num_thread][29] = 0.001736111111111111*b1-0.0003472222222222222*b0+0.006944444444444444*b10-0.001388888888888889*b11-0.002083333333333333*b12+0.01041666666666667*b13-0.02083333333333333*b14+0.02083333333333333*b15-0.01041666666666667*b16+0.002083333333333333*b17+0.001388888888888889*b18-0.006944444444444444*b19-0.003472222222222222*b2+0.01388888888888889*b20-0.01388888888888889*b21+0.006944444444444444*b22-0.001388888888888889*b23-0.0003472222222222222*b24+0.001736111111111111*b25-0.003472222222222222*b26+0.003472222222222222*b27-0.001736111111111111*b28+0.0003472222222222222*b29+0.003472222222222222*b3-0.001736111111111111*b4+0.0003472222222222222*b5+0.001388888888888889*b6-0.006944444444444444*b7+0.01388888888888889*b8-0.01388888888888889*b9;
                    QK_B_QKT_buffer[num_thread][30] = 0.0003472222222222222*b10-0.001805555555555556*b1-0.00006944444444444444*b0-0.0006944444444444444*b12-0.01805555555555556*b13-0.04583333333333333*b14-0.01805555555555556*b15-0.0006944444444444444*b16+0.0006944444444444444*b18+0.01805555555555556*b19-0.004583333333333333*b2+0.04583333333333333*b20+0.01805555555555556*b21+0.0006944444444444444*b22-0.0003472222222222222*b24-0.009027777777777778*b25-0.02291666666666667*b26-0.009027777777777778*b27-0.0003472222222222222*b28-0.001805555555555556*b3+0.00006944444444444444*b30+0.001805555555555556*b31+0.004583333333333333*b32+0.001805555555555556*b33+0.00006944444444444444*b34-0.00006944444444444444*b4+0.0003472222222222222*b6+0.009027777777777778*b7+0.02291666666666667*b8+0.009027777777777778*b9;
                    QK_B_QKT_buffer[num_thread][31] = 0.0003472222222222222*b0+0.003472222222222222*b1+0.001736111111111111*b10+0.003472222222222222*b12+0.03472222222222222*b13-0.03472222222222222*b15-0.003472222222222222*b16-0.003472222222222222*b18-0.03472222222222222*b19+0.03472222222222222*b21+0.003472222222222222*b22+0.001736111111111111*b24+0.01736111111111111*b25-0.01736111111111111*b27-0.001736111111111111*b28-0.003472222222222222*b3-0.0003472222222222222*b30-0.003472222222222222*b31+0.003472222222222222*b33+0.0003472222222222222*b34-0.0003472222222222222*b4-0.001736111111111111*b6-0.01736111111111111*b7+0.01736111111111111*b9;
                    QK_B_QKT_buffer[num_thread][32] = 0.003472222222222222*b10-0.001388888888888889*b1-0.0006944444444444444*b0-0.006944444444444444*b12-0.01388888888888889*b13+0.04166666666666667*b14-0.01388888888888889*b15-0.006944444444444444*b16+0.006944444444444444*b18+0.01388888888888889*b19+0.004166666666666667*b2-0.04166666666666667*b20+0.01388888888888889*b21+0.006944444444444444*b22-0.003472222222222222*b24-0.006944444444444444*b25+0.02083333333333333*b26-0.006944444444444444*b27-0.003472222222222222*b28-0.001388888888888889*b3+0.0006944444444444444*b30+0.001388888888888889*b31-0.004166666666666667*b32+0.001388888888888889*b33+0.0006944444444444444*b34-0.0006944444444444444*b4+0.003472222222222222*b6+0.006944444444444444*b7-0.02083333333333333*b8+0.006944444444444444*b9;
                    QK_B_QKT_buffer[num_thread][33] = 0.0006944444444444444*b0-0.001388888888888889*b1+0.003472222222222222*b10+0.006944444444444444*b12-0.01388888888888889*b13+0.01388888888888889*b15-0.006944444444444444*b16-0.006944444444444444*b18+0.01388888888888889*b19-0.01388888888888889*b21+0.006944444444444444*b22+0.003472222222222222*b24-0.006944444444444444*b25+0.006944444444444444*b27-0.003472222222222222*b28+0.001388888888888889*b3-0.0006944444444444444*b30+0.001388888888888889*b31-0.001388888888888889*b33+0.0006944444444444444*b34-0.0006944444444444444*b4-0.003472222222222222*b6+0.006944444444444444*b7-0.006944444444444444*b9;
                    QK_B_QKT_buffer[num_thread][34] = 0.001388888888888889*b1-0.0003472222222222222*b0+0.001736111111111111*b10-0.003472222222222222*b12+0.01388888888888889*b13-0.02083333333333333*b14+0.01388888888888889*b15-0.003472222222222222*b16+0.003472222222222222*b18-0.01388888888888889*b19-0.002083333333333333*b2+0.02083333333333333*b20-0.01388888888888889*b21+0.003472222222222222*b22-0.001736111111111111*b24+0.006944444444444444*b25-0.01041666666666667*b26+0.006944444444444444*b27-0.001736111111111111*b28+0.001388888888888889*b3+0.0003472222222222222*b30-0.001388888888888889*b31+0.002083333333333333*b32-0.001388888888888889*b33+0.0003472222222222222*b34-0.0003472222222222222*b4+0.001736111111111111*b6-0.006944444444444444*b7+0.01041666666666667*b8-0.006944444444444444*b9;
                    QK_B_QKT_buffer[num_thread][35] = 0.00006944444444444444*b0-0.0003472222222222222*b1-0.001736111111111111*b10+0.0003472222222222222*b11+0.0006944444444444444*b12-0.003472222222222222*b13+0.006944444444444444*b14-0.006944444444444444*b15+0.003472222222222222*b16-0.0006944444444444444*b17-0.0006944444444444444*b18+0.003472222222222222*b19+0.0006944444444444444*b2-0.006944444444444444*b20+0.006944444444444444*b21-0.003472222222222222*b22+0.0006944444444444444*b23+0.0003472222222222222*b24-0.001736111111111111*b25+0.003472222222222222*b26-0.003472222222222222*b27+0.001736111111111111*b28-0.0003472222222222222*b29-0.0006944444444444444*b3-0.00006944444444444444*b30+0.0003472222222222222*b31-0.0006944444444444444*b32+0.0006944444444444444*b33-0.0003472222222222222*b34+0.00006944444444444444*b35+0.0003472222222222222*b4-0.00006944444444444444*b5-0.0003472222222222222*b6+0.001736111111111111*b7-0.003472222222222222*b8+0.003472222222222222*b9;

                    // Interpolate value
                    g_buffer[num_thread][lind_g] = (QK_B_QKT_buffer[num_thread][0]+x_vec_buffer[num_thread][1]*QK_B_QKT_buffer[num_thread][6]+x_vec_buffer[num_thread][2]*QK_B_QKT_buffer[num_thread][12]+x_vec_buffer[num_thread][3]*QK_B_QKT_buffer[num_thread][18]+x_vec_buffer[num_thread][4]*QK_B_QKT_buffer[num_thread][24]+x_vec_buffer[num_thread][5]*QK_B_QKT_buffer[num_thread][30])+
                                                   (QK_B_QKT_buffer[num_thread][1]+x_vec_buffer[num_thread][1]*QK_B_QKT_buffer[num_thread][7]+x_vec_buffer[num_thread][2]*QK_B_QKT_buffer[num_thread][13]+x_vec_buffer[num_thread][3]*QK_B_QKT_buffer[num_thread][19]+x_vec_buffer[num_thread][4]*QK_B_QKT_buffer[num_thread][25]+x_vec_buffer[num_thread][5]*QK_B_QKT_buffer[num_thread][31])*y_vec_buffer[num_thread][1]+
                                                   (QK_B_QKT_buffer[num_thread][2]+x_vec_buffer[num_thread][1]*QK_B_QKT_buffer[num_thread][8]+x_vec_buffer[num_thread][2]*QK_B_QKT_buffer[num_thread][14]+x_vec_buffer[num_thread][3]*QK_B_QKT_buffer[num_thread][20]+x_vec_buffer[num_thread][4]*QK_B_QKT_buffer[num_thread][26]+x_vec_buffer[num_thread][5]*QK_B_QKT_buffer[num_thread][32])*y_vec_buffer[num_thread][2]+
                                                   (QK_B_QKT_buffer[num_thread][3]+x_vec_buffer[num_thread][1]*QK_B_QKT_buffer[num_thread][9]+x_vec_buffer[num_thread][2]*QK_B_QKT_buffer[num_thread][15]+x_vec_buffer[num_thread][3]*QK_B_QKT_buffer[num_thread][21]+x_vec_buffer[num_thread][4]*QK_B_QKT_buffer[num_thread][27]+x_vec_buffer[num_thread][5]*QK_B_QKT_buffer[num_thread][33])*y_vec_buffer[num_thread][3]+
                                                   (QK_B_QKT_buffer[num_thread][4]+x_vec_buffer[num_thread][1]*QK_B_QKT_buffer[num_thread][10]+x_vec_buffer[num_thread][2]*QK_B_QKT_buffer[num_thread][16]+x_vec_buffer[num_thread][3]*QK_B_QKT_buffer[num_thread][22]+x_vec_buffer[num_thread][4]*QK_B_QKT_buffer[num_thread][28]+x_vec_buffer[num_thread][5]*QK_B_QKT_buffer[num_thread][34])*y_vec_buffer[num_thread][4]+
                                                   (QK_B_QKT_buffer[num_thread][5]+x_vec_buffer[num_thread][1]*QK_B_QKT_buffer[num_thread][11]+x_vec_buffer[num_thread][2]*QK_B_QKT_buffer[num_thread][17]+x_vec_buffer[num_thread][3]*QK_B_QKT_buffer[num_thread][23]+x_vec_buffer[num_thread][4]*QK_B_QKT_buffer[num_thread][29]+x_vec_buffer[num_thread][5]*QK_B_QKT_buffer[num_thread][35])*y_vec_buffer[num_thread][5];
                    
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
    
    // Debug
    // if (num_thread == 0) {
    //       imshow(&g_buffer[num_thread][0],radius*2+1,radius*2+1);
    // }

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
            for (int j=0; j<roi[0].cirroi[num_thread].region.noderange.value[i]; j+=2) {
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

        // Return success
        return SUCCESS;   
    }
    // Deltag was close to zero - return failed
    return FAILED;
}

void mexFunction(int nlhs,mxArray *plhs[ ],int nrhs,const mxArray *prhs[ ]) {
    if (nrhs == 10 && nlhs == 3) {
        // Create calcseeds 
        class_calcseeds calcseeds(plhs,prhs);
        
        // Run analysis
        calcseeds.analysis();
    } else {
        // Thread safe because it is single threaded up to this point
        mexErrMsgTxt("Incorrect number of inputs or outputs.\n");
    }
}
