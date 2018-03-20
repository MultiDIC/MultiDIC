// These are function definitions used to obtain ncorr-type inputs to mex files.
// Some functions are not thread safe.

#include <mex.h>
#include <math.h>
#include <algorithm>
#include "standard_datatypes.h"
#include "ncorr_datatypes.h"
#include "ncorr_lib.h"

// ----------------------------------------------------------------------//
// NCORR_CLASS_IMG STUFF ------------------------------------------------//
// ----------------------------------------------------------------------//

// THREADSAFE
ncorr_class_img::ncorr_class_img() {
    // gs uses default constructor for struct_double_array
    max_gs = 0.0;
    // bcoef uses default constructor for struct_double_array
    border_bcoef = 0;
}

// NOT THREADSAFE
void get_imgs(std::vector<ncorr_class_img> &images,const mxArray *prhs) {
    // Check inputs
    if ((mxIsClass(prhs,"ncorr_class_img") || mxIsClass(prhs,"struct")) && mxGetM(prhs) == 1 &&  mxGetN(prhs) >= 1) {    
        for (int i=0; i<(int)mxGetN(prhs); i++) {
            ncorr_class_img img_template;
            mxArray *mat_gs;
            mxArray *mat_max_gs;
            mxArray *mat_bcoef;
            mxArray *mat_border_bcoef;
            mxArray *mat_type;
            if (mxIsClass(prhs,"ncorr_class_img")) {
                mat_gs = mxGetProperty(prhs,i,"gs");
                mat_max_gs = mxGetProperty(prhs,i,"max_gs");
                mat_bcoef = mxGetProperty(prhs,i,"bcoef");
                mat_border_bcoef = mxGetProperty(prhs,i,"border_bcoef");
                mat_type = mxGetProperty(prhs,i,"type");
            } else {
                mat_gs = mxGetField(prhs,i,"gs");
                mat_max_gs = mxGetField(prhs,i,"max_gs");
                mat_bcoef = mxGetField(prhs,i,"bcoef");
                mat_border_bcoef = mxGetField(prhs,i,"border_bcoef");
                mat_type = mxGetField(prhs,i,"type");
            }
            if (mat_gs != 0 && mat_max_gs != 0 && mat_bcoef != 0 && mat_border_bcoef != 0 && mat_type != 0) {
                // Get data
                get_double_array(img_template.gs,mat_gs);
                get_double_scalar(img_template.max_gs,mat_max_gs);
                get_string(img_template.type,mat_type);
                
                // If image is not reduced then get interpolation info
                if (img_template.type.compare(0,7,"reduced") != 0) {
                    get_double_array(img_template.bcoef,mat_bcoef);     
                    get_integer_scalar(img_template.border_bcoef,mat_border_bcoef);
                }
                
                // Store
                images.push_back(img_template);
            } else {
                mexErrMsgTxt("Properties 'gs', 'max_gs', 'bcoef', 'bcoef_border', or 'type' do not exist for image.\n");
            }  
        }
    } else {
         mexErrMsgTxt("Image is not of type 'ncorr_class_img' or 'struct' and size == [1 N].\n");   
    }
}

// ----------------------------------------------------------------------//
// NCORR_CLASS_ROI STUFF ------------------------------------------------//
// ----------------------------------------------------------------------//

// THREADSAFE
ncorr_class_region::ncorr_class_region() {    
    upperbound = 0;
    lowerbound = 0;
    leftbound = 0;
    rightbound = 0;
    totalpoints = 0;
}

// NOT THREADSAFE
void ncorr_class_region::alloc(const int &h,const int &w) {    
    // Allocate space for nodelist and noderange
    nodelist.alloc(h,w);
    noderange.alloc(h,1);
}

// NOT THREADSAFE
void ncorr_class_region::free() {
    // Free nodelist and noderange
    nodelist.free();
    noderange.free();
}

// NOT THREADSAFE
void get_region(std::vector<ncorr_class_region> &region,const mxArray *prhs) {
    // Check input    
    if (mxIsClass(prhs,"struct") && mxGetM(prhs) == 1) {
        for (int i=0; i<(int)mxGetN(prhs); i++) {
            mxArray *mat_nodelist = mxGetField(prhs,i,"nodelist");
            mxArray *mat_noderange = mxGetField(prhs,i,"noderange");
            mxArray *mat_upperbound = mxGetField(prhs,i,"upperbound");
            mxArray *mat_lowerbound = mxGetField(prhs,i,"lowerbound");
            mxArray *mat_leftbound = mxGetField(prhs,i,"leftbound");
            mxArray *mat_rightbound = mxGetField(prhs,i,"rightbound");
            mxArray *mat_totalpoints = mxGetField(prhs,i,"totalpoints");
            if (mat_nodelist != 0 &&
                mat_noderange != 0 &&
                mat_upperbound != 0 &&
                mat_lowerbound != 0 &&
                mat_leftbound != 0 &&
                mat_rightbound != 0 &&
                mat_totalpoints != 0) {
                // Form region
                ncorr_class_region region_template;
                get_integer_array(region_template.nodelist,mat_nodelist);
                get_integer_array(region_template.noderange,mat_noderange);
                get_integer_scalar(region_template.upperbound,mat_upperbound);
                get_integer_scalar(region_template.lowerbound,mat_lowerbound);
                get_integer_scalar(region_template.leftbound,mat_leftbound);
                get_integer_scalar(region_template.rightbound,mat_rightbound);
                get_integer_scalar(region_template.totalpoints,mat_totalpoints);
                
                // Check inputs
                // Check noderange:
                if (region_template.noderange.height != region_template.nodelist.height || region_template.noderange.width != 1) {
                    mexErrMsgTxt("'nodelist' or 'noderange' size is incorrect.\n");
                }
                
                // Store region_template
                region.push_back(region_template);
            } else {
                mexErrMsgTxt("Some fields do not exist for region.\n");
            }
        }
    } else {
        mexErrMsgTxt("Input needs to be of type 'struct'.\n");
    }
}

// Cirroi stuff ---------------------------------------------------------//
// THREADSAFE
struct_cirroi::struct_cirroi() {
    radius = 0;
    x = 0;
    y = 0;
}

// THREADSAFE
struct_info_cirroi::struct_info_cirroi() {
    radius = 0;
    x = 0;
    y = 0;
}

// This is for an ncorr_class_roi input ---------------------------------//
// THREADSAFE
ncorr_class_roi::ncorr_class_roi() {
}

// NOT THREADSAFE
void get_rois(std::vector<ncorr_class_roi> &rois,const mxArray *prhs) {
    // Check if input is an object
    if ((mxIsClass(prhs,"ncorr_class_roi") || mxIsClass(prhs,"struct")) && mxGetM(prhs) == 1 && mxGetN(prhs) >= 1) {   
        for (int i=0; i<(int)mxGetN(prhs); i++) {
            ncorr_class_roi roi_template;
            mxArray *mat_mask;
            mxArray *mat_region;
            if (mxIsClass(prhs,"ncorr_class_roi")) {
                mat_mask = mxGetProperty(prhs,i,"mask");
                mat_region = mxGetProperty(prhs,i,"region");
            } else {
                mat_mask = mxGetField(prhs,i,"mask");
                mat_region = mxGetField(prhs,i,"region");
            }
            if (mat_mask != 0 && mat_region != 0) {
                // Get data            
                get_logical_array(roi_template.mask,mat_mask);
                get_region(roi_template.region,mat_region);
                
                // Store
                rois.push_back(roi_template);                                
            } else {
                mexErrMsgTxt("'region' or 'mask' is not a property of roi.\n");
            }
        }
    } else {
        mexErrMsgTxt("Input needs to be of type 'ncorr_class_roi' or 'struct' and size == [1 N].\n");
    }
}

// NOT THREADSAFE
void ncorr_class_roi::set_cirroi(const int &radius,const int &total_threads) {
    // Used to "set" the cirroi. Must be called before calling update_cirroi, which is followed by get_cirroi
    // This is called only once    
    // Get max nodewidth
    int max_nodewidth = 0;
    for (int i=0; i<(int)region.size(); i++) {
        if (region[i].nodelist.width > max_nodewidth) {
            max_nodewidth = region[i].nodelist.width;
        }    
    }
    
    // Set cirroi and info_cirroi - one for each thread
    cirroi.resize(total_threads);
    info_cirroi.resize(total_threads);
    
    for (int i=0; i<total_threads; i++) {
        info_cirroi[i].region.nodelist.alloc(2*radius+1,max_nodewidth);
        info_cirroi[i].region.noderange.alloc(2*radius+1,1); 
        info_cirroi[i].mask.alloc(2*radius+1,2*radius+1);     
        info_cirroi[i].radius = radius;
        info_cirroi[i].x = 0;
        info_cirroi[i].y = 0;
        
        info_cirroi[i].mask_buffer.alloc(2*radius+1,2*radius+1);
        info_cirroi[i].circletemplate.alloc(2*radius+1,2);     
        info_cirroi[i].queue_buffer.alloc(max_nodewidth,1);         
        
        // Assign space for queue vectors
        info_cirroi[i].queue_nodelist.reserve(max_nodewidth*(2*radius+1));
        info_cirroi[i].queue_nodeindex.reserve(max_nodewidth*(2*radius+1)/2);
        
        // Form circle template
        for (int j=0; j<info_cirroi[i].circletemplate.height; j++) {
            int h = j-info_cirroi[i].radius;
            int top = (int)ceil(-sqrt(pow((double)info_cirroi[i].radius,2)-pow((double)h,2)));
            int bottom = (int)floor(sqrt(pow((double)info_cirroi[i].radius,2)-pow((double)h,2)));
            info_cirroi[i].circletemplate.value[j] = top;
            info_cirroi[i].circletemplate.value[j+info_cirroi[i].circletemplate.height] = bottom;
        }
    }
}

// THREADSAFE
void ncorr_class_roi::update_cirroi(const int &num_region,const int &num_thread) {
    // This function is called once per region before calling get_cirroi.
    if (info_cirroi[num_thread].region.nodelist.value != NULL && info_cirroi[num_thread].region.noderange.value != NULL && info_cirroi[num_thread].mask.value != NULL) {      
        // Resize activelines
        info_cirroi[num_thread].activelines.resize(region[num_region].nodelist.height*region[num_region].nodelist.width/2);
    } else {
        // Cannot use mexErrMsgTxt because it is not thread safe. mexPrintf is only marginally thread safe.
        mexPrintf("Cirroi has not been set yet for thread %d; this will most likely cause a crash, hopefully this message is seen before matlab closes.\n",num_thread);
    }      
}

// THREADSAFE
void ncorr_class_roi::get_cirroi(const int &x,const int &y,const int &num_region,const bool &subsettrunc,const int &num_thread) {
    // Function creates cirroi which can be accessed through the cirroi member    
    // Quick check to see if cirroi has been set yet:
    if (info_cirroi[num_thread].region.nodelist.value != NULL && info_cirroi[num_thread].region.noderange.value != NULL && info_cirroi[num_thread].mask.value != NULL) {
        // Reset Values -------------------------------------------------//
        // Reset mask
        info_cirroi[num_thread].mask.reset();

        // Reset active lines - set true
        std::fill(info_cirroi[num_thread].activelines.begin(),info_cirroi[num_thread].activelines.end(),1);
                
        // Reset noderange
        info_cirroi[num_thread].region.noderange.reset();
        
        // Reset totalpoints
        info_cirroi[num_thread].region.totalpoints = 0;
        
        // Set x --------------------------------------------------------//
        info_cirroi[num_thread].x = x;

        // Update Circle template ---------------------------------------//
        // Subtract old y -> add new y 
        // Steps 1&2:
        for (int i=0; i<info_cirroi[num_thread].circletemplate.height; i++) {
            info_cirroi[num_thread].circletemplate.value[i] += -info_cirroi[num_thread].y+y;
            info_cirroi[num_thread].circletemplate.value[i+info_cirroi[num_thread].circletemplate.height] += -info_cirroi[num_thread].y+y;
        }
        
        // Set y --------------------------------------------------------//
        info_cirroi[num_thread].y = y;
        
        // Initialize parameter which tells us whether 
        // or not the subset has been truncated
        bool circ_untrunc = true;

        // Set ROI and range --------------------------------------------//
        // Find idx_nodelist
        int idx_nodelist = -1;
        int idx_roi_x = info_cirroi[num_thread].x-region[num_region].leftbound;
        if (idx_roi_x >= 0 && idx_roi_x < region[num_region].noderange.height) {
            for (int i=0; i<region[num_region].noderange.value[idx_roi_x]; i+=2) {
                if (info_cirroi[num_thread].y >= region[num_region].nodelist.value[idx_roi_x+i*region[num_region].nodelist.height] && info_cirroi[num_thread].y <= region[num_region].nodelist.value[idx_roi_x+(i+1)*region[num_region].nodelist.height]) {
                    idx_nodelist = i;
                    break;
                }
            }
        }
        if (idx_nodelist == -1) {        
            // Cannot use mexErrMsgTxt because it is not thread safe. mexPrintf is only marginally thread safe.
            mexPrintf("(X,Y) are not in region.\n");
        }

        // Add center nodes to queue
        int node_top;
        int node_bottom;
        if (region[num_region].nodelist.value[idx_roi_x+idx_nodelist*region[num_region].nodelist.height] < info_cirroi[num_thread].circletemplate.value[info_cirroi[num_thread].radius]) {
            node_top = info_cirroi[num_thread].circletemplate.value[info_cirroi[num_thread].radius];
        } else {
            node_top = region[num_region].nodelist.value[idx_roi_x+idx_nodelist*region[num_region].nodelist.height];
            circ_untrunc = false;
        }
        if (region[num_region].nodelist.value[idx_roi_x+(idx_nodelist+1)*region[num_region].nodelist.height] > info_cirroi[num_thread].circletemplate.value[info_cirroi[num_thread].radius+info_cirroi[num_thread].circletemplate.height]) {
            node_bottom = info_cirroi[num_thread].circletemplate.value[info_cirroi[num_thread].radius+info_cirroi[num_thread].circletemplate.height];
        } else {
            node_bottom = region[num_region].nodelist.value[idx_roi_x+(idx_nodelist+1)*region[num_region].nodelist.height];
            circ_untrunc = false;
        }
        
        // Update mask
        for (int i=node_top-(info_cirroi[num_thread].y-info_cirroi[num_thread].radius); i<=node_bottom-(info_cirroi[num_thread].y-info_cirroi[num_thread].radius); i++) {
            info_cirroi[num_thread].mask.value[i+info_cirroi[num_thread].radius*info_cirroi[num_thread].mask.height] = true;
        }
        
        // Update queue
        info_cirroi[num_thread].queue_nodelist.push_back(node_top);
        info_cirroi[num_thread].queue_nodelist.push_back(node_bottom);
        info_cirroi[num_thread].queue_nodeindex.push_back(info_cirroi[num_thread].radius);
        
        // Inactivate nodepair from region
        info_cirroi[num_thread].activelines[idx_roi_x+(idx_nodelist/2)*(region[num_region].nodelist.height)] = 0; //Set to false

        // Enter while loop until queue is empty
        while (info_cirroi[num_thread].queue_nodelist.size() > 0) {
            // STEPS:
            // 1) Load nodepair from queue
            // 2) Update queue
            // 3) Add nodepair to cirroi and sort
            // 4) Compare nodepair to nodepairs left and right and add nodes to queue
            // 5) Update totalpoints

            // Steps 1 & 2: Do top and bottom in reverse order
            int queue_nodelist_bottom_buffer = info_cirroi[num_thread].queue_nodelist.back();
            info_cirroi[num_thread].queue_nodelist.pop_back();
            int queue_nodelist_top_buffer = info_cirroi[num_thread].queue_nodelist.back();
            info_cirroi[num_thread].queue_nodelist.pop_back();
            int queue_nodeindex_buffer = info_cirroi[num_thread].queue_nodeindex.back();
            info_cirroi[num_thread].queue_nodeindex.pop_back();

            // Step 3:          
            if (info_cirroi[num_thread].region.noderange.value[queue_nodeindex_buffer] == 0) {
                info_cirroi[num_thread].region.nodelist.value[queue_nodeindex_buffer] = queue_nodelist_top_buffer;
                info_cirroi[num_thread].region.nodelist.value[queue_nodeindex_buffer+info_cirroi[num_thread].region.nodelist.height] = queue_nodelist_bottom_buffer;
            } else {
                bool inserted = false;
                for (int i=0; i<info_cirroi[num_thread].region.noderange.value[queue_nodeindex_buffer]; i+=2) {
                    if (queue_nodelist_bottom_buffer < info_cirroi[num_thread].region.nodelist.value[queue_nodeindex_buffer+i*info_cirroi[num_thread].region.nodelist.height]) {
                        for (int j=i; j<info_cirroi[num_thread].region.noderange.value[queue_nodeindex_buffer]; j++) {
                            info_cirroi[num_thread].queue_buffer.value[j-i] = info_cirroi[num_thread].region.nodelist.value[queue_nodeindex_buffer+j*info_cirroi[num_thread].region.nodelist.height];
                        }
                        
                        info_cirroi[num_thread].region.nodelist.value[queue_nodeindex_buffer+i*info_cirroi[num_thread].region.nodelist.height] = queue_nodelist_top_buffer;
                        info_cirroi[num_thread].region.nodelist.value[queue_nodeindex_buffer+(i+1)*info_cirroi[num_thread].region.nodelist.height] = queue_nodelist_bottom_buffer;
                        
                        for (int j=i+2; j<info_cirroi[num_thread].region.noderange.value[queue_nodeindex_buffer]+2; j++) {
                            info_cirroi[num_thread].region.nodelist.value[queue_nodeindex_buffer+j*info_cirroi[num_thread].region.nodelist.height] = info_cirroi[num_thread].queue_buffer.value[j-(i+2)];
                        }
                        
                        inserted = true;
                        break;
                    }
                }
                
                if (!inserted){
                    info_cirroi[num_thread].region.nodelist.value[queue_nodeindex_buffer+info_cirroi[num_thread].region.noderange.value[queue_nodeindex_buffer]*info_cirroi[num_thread].region.nodelist.height] = queue_nodelist_top_buffer;
                    info_cirroi[num_thread].region.nodelist.value[queue_nodeindex_buffer+(info_cirroi[num_thread].region.noderange.value[queue_nodeindex_buffer]+1)*info_cirroi[num_thread].region.nodelist.height] = queue_nodelist_bottom_buffer;
                }
            }
            info_cirroi[num_thread].region.noderange.value[queue_nodeindex_buffer] = info_cirroi[num_thread].region.noderange.value[queue_nodeindex_buffer]+2;

            // Step 4:
            int idx_froi = queue_nodeindex_buffer+(info_cirroi[num_thread].x-info_cirroi[num_thread].radius)-region[num_region].leftbound;

            // Check nodepairs LEFT idx_froi (i.e. subtract 1)
            if (queue_nodeindex_buffer-1 >= 0 && idx_froi-1 >= 0) {
                for (int i=0; i<region[num_region].noderange.value[idx_froi-1]; i+=2) {
                    if (info_cirroi[num_thread].activelines[idx_froi-1+(i/2)*(region[num_region].nodelist.height)] == 0) {
                        continue;
                    }

                    if (region[num_region].nodelist.value[(idx_froi-1)+(i+1)*(region[num_region].nodelist.height)] < queue_nodelist_top_buffer) {
                        continue;
                    } else if (region[num_region].nodelist.value[(idx_froi-1)+i*(region[num_region].nodelist.height)] <= queue_nodelist_bottom_buffer && region[num_region].nodelist.value[(idx_froi-1)+(i+1)*(region[num_region].nodelist.height)] >= queue_nodelist_top_buffer) {

                        if (region[num_region].nodelist.value[(idx_froi-1)+i*(region[num_region].nodelist.height)] <  info_cirroi[num_thread].circletemplate.value[queue_nodeindex_buffer-1]) {
                            node_top = info_cirroi[num_thread].circletemplate.value[queue_nodeindex_buffer-1];
                        } else {
                            node_top = region[num_region].nodelist.value[(idx_froi-1)+i*(region[num_region].nodelist.height)];
                            circ_untrunc = false;
                        }    

                        if (region[num_region].nodelist.value[(idx_froi-1)+(i+1)*(region[num_region].nodelist.height)] > info_cirroi[num_thread].circletemplate.value[(queue_nodeindex_buffer-1)+info_cirroi[num_thread].circletemplate.height]) {
                            node_bottom = info_cirroi[num_thread].circletemplate.value[(queue_nodeindex_buffer-1)+info_cirroi[num_thread].circletemplate.height];
                        } else {
                            node_bottom = region[num_region].nodelist.value[(idx_froi-1)+(i+1)*(region[num_region].nodelist.height)];
                            circ_untrunc = false;
                        }    

                        if (node_top > node_bottom || node_top > queue_nodelist_bottom_buffer || node_bottom < queue_nodelist_top_buffer) {
                            continue;
                        }

                        // Add nodes to queue
                        info_cirroi[num_thread].queue_nodelist.push_back(node_top);
                        info_cirroi[num_thread].queue_nodelist.push_back(node_bottom);
                        info_cirroi[num_thread].queue_nodeindex.push_back(queue_nodeindex_buffer-1);

                        // Make node pair inactive
                        info_cirroi[num_thread].activelines[idx_froi-1+(i/2)*(region[num_region].nodelist.height)] = 0;
                        
                        // Update mask
                        for (int j=node_top-(info_cirroi[num_thread].y-info_cirroi[num_thread].radius); j<=node_bottom-(info_cirroi[num_thread].y-info_cirroi[num_thread].radius); j++) {
                            info_cirroi[num_thread].mask.value[j+(queue_nodeindex_buffer-1)*info_cirroi[num_thread].mask.height] = true;
                        }
                    } else {
                        break;
                    }
                }
            }

            // Check nodepairs RIGHT idx_froi (i.e. subtract 1)
            if (queue_nodeindex_buffer+1 <= 2*info_cirroi[num_thread].radius && idx_froi+1 <= region[num_region].rightbound-region[num_region].leftbound) {
                for (int i=0; i<region[num_region].noderange.value[idx_froi+1]; i+=2) {
                    if (info_cirroi[num_thread].activelines[idx_froi+1+(i/2)*(region[num_region].nodelist.height)] == 0) {
                        continue;
                    }

                    if (region[num_region].nodelist.value[(idx_froi+1)+(i+1)*(region[num_region].nodelist.height)] < queue_nodelist_top_buffer) {
                        continue;
                    } else if (region[num_region].nodelist.value[(idx_froi+1)+i*(region[num_region].nodelist.height)] <= queue_nodelist_bottom_buffer && region[num_region].nodelist.value[(idx_froi+1)+(i+1)*(region[num_region].nodelist.height)] >= queue_nodelist_top_buffer) {

                        if (region[num_region].nodelist.value[(idx_froi+1)+i*(region[num_region].nodelist.height)] <  info_cirroi[num_thread].circletemplate.value[queue_nodeindex_buffer+1]) {
                            node_top = info_cirroi[num_thread].circletemplate.value[queue_nodeindex_buffer+1];
                        } else {
                            node_top = region[num_region].nodelist.value[(idx_froi+1)+i*(region[num_region].nodelist.height)];
                            circ_untrunc = false;
                        }    

                        if (region[num_region].nodelist.value[(idx_froi+1)+(i+1)*(region[num_region].nodelist.height)] > info_cirroi[num_thread].circletemplate.value[(queue_nodeindex_buffer+1)+info_cirroi[num_thread].circletemplate.height]) {
                            node_bottom = info_cirroi[num_thread].circletemplate.value[(queue_nodeindex_buffer+1)+info_cirroi[num_thread].circletemplate.height];
                        } else {
                            node_bottom = region[num_region].nodelist.value[(idx_froi+1)+(i+1)*(region[num_region].nodelist.height)];
                            circ_untrunc = false;
                        }    

                        if (node_top > node_bottom || node_top > queue_nodelist_bottom_buffer || node_bottom < queue_nodelist_top_buffer) {
                            continue;
                        }

                        // Add nodes to queue
                        info_cirroi[num_thread].queue_nodelist.push_back(node_top);
                        info_cirroi[num_thread].queue_nodelist.push_back(node_bottom);
                        info_cirroi[num_thread].queue_nodeindex.push_back(queue_nodeindex_buffer+1);

                        // Make node pair inactive
                        info_cirroi[num_thread].activelines[idx_froi+1+(i/2)*(region[num_region].nodelist.height)] = 0;
                        
                        // Update mask
                        for (int j=node_top-(info_cirroi[num_thread].y-info_cirroi[num_thread].radius); j<=node_bottom-(info_cirroi[num_thread].y-info_cirroi[num_thread].radius); j++) {
                            info_cirroi[num_thread].mask.value[j+(queue_nodeindex_buffer+1)*info_cirroi[num_thread].mask.height] = true;
                        }
                    } else {
                        break;
                    }
                }
            }

            // Update totalpoints
            info_cirroi[num_thread].region.totalpoints += queue_nodelist_bottom_buffer-queue_nodelist_top_buffer+1;
        }
        
        // Set bounds ---------------------------------------------------//
        // Left and right bounds - Since nodelist/noderange is always 2*radius+1, set left and right bounds:
        info_cirroi[num_thread].region.leftbound = info_cirroi[num_thread].x-info_cirroi[num_thread].radius;
        info_cirroi[num_thread].region.rightbound = info_cirroi[num_thread].x+info_cirroi[num_thread].radius;
        
        // Find Upper and lower bounds
        info_cirroi[num_thread].region.upperbound = mask.height;
        info_cirroi[num_thread].region.lowerbound = 0;
        for (int i=0; i<info_cirroi[num_thread].region.noderange.height; i++) {
            if (info_cirroi[num_thread].region.noderange.value[i] > 0) {
                if (info_cirroi[num_thread].region.upperbound > info_cirroi[num_thread].region.nodelist.value[i]) {
                    info_cirroi[num_thread].region.upperbound = info_cirroi[num_thread].region.nodelist.value[i];
                }                
                if (info_cirroi[num_thread].region.lowerbound < info_cirroi[num_thread].region.nodelist.value[i+(info_cirroi[num_thread].region.noderange.value[i]-1)*info_cirroi[num_thread].region.nodelist.height]) {
                    info_cirroi[num_thread].region.lowerbound = info_cirroi[num_thread].region.nodelist.value[i+(info_cirroi[num_thread].region.noderange.value[i]-1)*info_cirroi[num_thread].region.nodelist.height];
                }
            }
        }
        
        // See if there are any noderanges with zero value, these will
        // not trigger circ_untruc so you must check it here
        for (int i=0; i<info_cirroi[num_thread].region.noderange.height; i++) {
            if (info_cirroi[num_thread].region.noderange.value[i] == 0) {
                circ_untrunc = false;
                break;
            }
        }
        
        // Perform subset truncation ------------------------------------//
        // Check to see if this feature is enabled.
        if (subsettrunc && !circ_untrunc) {
            // Initialize vec_boundary
            std::vector<std::vector<int> > vec_boundary;
            
            // Find the top left point
            std::vector<int> vec_point_topleft;
            vec_point_topleft.resize(2,-1);
            for (int i=0; i<info_cirroi[num_thread].mask.width; i++) {
                for (int j=0; j<info_cirroi[num_thread].mask.height; j++) {
                    if (info_cirroi[num_thread].mask.value[j+i*info_cirroi[num_thread].mask.height]) {
                        vec_point_topleft[0] = i;
                        vec_point_topleft[1] = j;
                        break;
                    }
                }
                if (vec_point_topleft[0] != -1 && vec_point_topleft[1] != -1) {
                    break;
                }
            }
                                          
            // Find boundary
            int direc = 7;
            form_boundary(vec_boundary,vec_point_topleft,info_cirroi[num_thread].mask,direc);
            
            // Find closest point -> Find line equation -> Find 
            // number of points to the left and right of the 
            // line; the side with the least amount of points 
            // will be discarded
            int idx_min = 0; // Start with first boundary point
            double val_min = pow((double)(vec_boundary[0][0]-info_cirroi[num_thread].radius),2.0) + pow((double)(vec_boundary[0][1]-info_cirroi[num_thread].radius),2.0);
            // Start at 1 since first value has been taken
            for (int i=1; i<(int)vec_boundary.size(); i++) {
                double val_min_buffer = pow((double)(vec_boundary[i][0]-info_cirroi[num_thread].radius),2.0) + pow((double)(vec_boundary[i][1]-info_cirroi[num_thread].radius),2.0);
                if (val_min_buffer < val_min) {
                    val_min = val_min_buffer;
                    idx_min = i;
                }
            }
        
            // Make sure closest point is not on the boundary, which can
            // happen since the boundary isnt perfectly circular
            if ((int)ceil(sqrt(val_min)) < info_cirroi[num_thread].radius) {
                // Get two points on either side of the closest point. Use a 
                // nonlinear solver with a single step to determine the two points
                std::vector<double> p0(2,0);
                std::vector<double> p1(2,0);

                // Get indices of points 1 behind and 1 ahead   
                int idx_space = 3;
                int idx_plus = mod_pos(idx_min+idx_space,(int)vec_boundary.size());  
                int idx_minus = mod_pos(idx_min-idx_space,(int)vec_boundary.size());

                // Initialize points to calculate derivatives with
                double x_plus_f = 0;
                double x_min_f = 0; 
                double x_minus_f = 0;
                double y_plus_f = 0;
                double y_min_f = 0; 
                double y_minus_f = 0;
                // Set filt length
                int length_filt = 2;
                for (int i=-length_filt; i<=length_filt; i++) {
                    x_plus_f += (double)vec_boundary[mod_pos(idx_plus+i,(int)vec_boundary.size())][0];
                    x_min_f += (double)vec_boundary[mod_pos(idx_min+i,(int)vec_boundary.size())][0];
                    x_minus_f += (double)vec_boundary[mod_pos(idx_minus+i,(int)vec_boundary.size())][0];

                    y_plus_f += (double)vec_boundary[mod_pos(idx_plus+i,(int)vec_boundary.size())][1];
                    y_min_f += (double)vec_boundary[mod_pos(idx_min+i,(int)vec_boundary.size())][1];
                    y_minus_f += (double)vec_boundary[mod_pos(idx_minus+i,(int)vec_boundary.size())][1];
                }            
                // Divide by length
                x_plus_f /= 2*length_filt+1;
                x_min_f /= 2*length_filt+1;
                x_minus_f /= 2*length_filt+1;
                y_plus_f /= 2*length_filt+1;
                y_min_f /= 2*length_filt+1; 
                y_minus_f /= 2*length_filt+1;   

                // Get derivatives
                double dx_di_f = (x_plus_f - x_minus_f)/(2*idx_space);
                double d2x_di2_f = (x_plus_f - 2*x_min_f + x_minus_f)/pow((double)idx_space,2);
                double dy_di_f = (y_plus_f - y_minus_f)/(2*idx_space);
                double d2y_di2_f = (y_plus_f - 2*y_min_f + y_minus_f)/pow((double)idx_space,2);

                // Do one iteration to find approximation of index
                double deltai = ((-x_min_f*dx_di_f+(double)info_cirroi[num_thread].radius*dx_di_f) + (-y_min_f*dy_di_f+(double)info_cirroi[num_thread].radius*dy_di_f))/((pow(dx_di_f,2)+x_min_f*d2x_di2_f-(double)info_cirroi[num_thread].radius*d2x_di2_f)+(pow(dy_di_f,2)+y_min_f*d2y_di2_f-(double)info_cirroi[num_thread].radius*d2y_di2_f));

                // Test deltai, make sure its between subsequent spacings; If
                // not then set the points minus and plus points  
                if (fabs(deltai) < idx_space) {    
                    // Calculate two points based on dx_di_i and dy_di_i
                    double x_i = x_min_f + dx_di_f*deltai + (1.0/2.0)*d2x_di2_f*pow(deltai,2); 
                    double y_i = y_min_f + dy_di_f*deltai + (1.0/2.0)*d2y_di2_f*pow(deltai,2); 
                    double dx_di_i = dx_di_f + d2x_di2_f*deltai;
                    double dy_di_i = dy_di_f + d2y_di2_f*deltai;

                    // Set stepsize - magnitude doesnt matter
                    double stepsize = 0.5;

                    // Get points
                    p0[0] = x_i - dx_di_i*stepsize;
                    p0[1] = y_i - dy_di_i*stepsize;
                    p1[0] = x_i + dx_di_i*stepsize;
                    p1[1] = y_i + dy_di_i*stepsize;
                } else {
                    // Just set points to averaged points on either side
                    p0[0] = x_minus_f;
                    p0[1] = y_minus_f;
                    p1[0] = x_plus_f;
                    p1[1] = y_plus_f;                
                }

                // Find which side to clear: p_subset is a point 
                // defined to be part of the subset, so the side
                // opposite of p_subset needs to be cleared. If the center 
                // of the subset does not lie on boundary, then find
                // displacement from closest point to the p0. Add
                // this displacement to the center, and then
                // determine which side the center is on.
                // Displacement from closest point to p0:
                std::vector<double> p_subset(2,0);

                if (vec_boundary[idx_min][0] == info_cirroi[num_thread].radius && vec_boundary[idx_min][1] == info_cirroi[num_thread].radius) {
                    // Center is the closest point. Find valid points around 
                    // the boundary and get the centroid of these points.
                    int width_win = 1; // This determines window of points collected
                    int counter = 0; // Counts number of points added so a proper average can be taken
                    for (int i=-width_win; i<=width_win; i++) {
                        int x_mask = vec_boundary[idx_min][0]+i;
                        for (int j=-width_win;j <= width_win; j++ ) {
                            int y_mask = vec_boundary[idx_min][1]+j;

                            // Make sure points are within mask and
                            // that the mask is valid here.
                            if (x_mask >= 0 && x_mask <= 2*info_cirroi[num_thread].radius && y_mask >= 0 && y_mask <= 2*info_cirroi[num_thread].radius && 
                                info_cirroi[num_thread].mask.value[y_mask+x_mask*info_cirroi[num_thread].mask.height]) {
                                p_subset[0] += x_mask;
                                p_subset[1] += y_mask;
                                ++counter;
                            }
                        }
                    }                           
                    // Divide by counter to get the average
                    p_subset[0] /= counter;
                    p_subset[1] /= counter;
                } else {
                    // Set the centerpoint 
                    p_subset[0] = info_cirroi[num_thread].radius;
                    p_subset[1] = info_cirroi[num_thread].radius;                
                }

                // Normalize point against line since it can shift
                p_subset[0] += p0[0] - vec_boundary[idx_min][0];
                p_subset[1] += p0[1] - vec_boundary[idx_min][1];

                // Now determine which side p_subset lies and
                // clear the other side
                int sign_clear = -sign((p1[0]-p0[0])*(p_subset[1]-p0[1])-(p_subset[0]-p0[0])*(p1[1]-p0[1]));  

                // Make a copy of the mask - this is because truncation can potentially
                // cause the mask to be empty. In this case simply leave the subset
                // untruncated
                for (int i=0; i<info_cirroi[num_thread].mask_buffer.width; i++) { 
                    for (int j=0; j<info_cirroi[num_thread].mask_buffer.height; j++) { 
                        info_cirroi[num_thread].mask_buffer.value[j+i*info_cirroi[num_thread].mask_buffer.height] = info_cirroi[num_thread].mask.value[j+i*info_cirroi[num_thread].mask.height];
                    }
                }

                for (int i=0; i<info_cirroi[num_thread].region.noderange.height; i++) {  
                    int x = i;
                    for (int j=0; j<info_cirroi[num_thread].region.noderange.value[i]; j+=2) {
                        for (int k=info_cirroi[num_thread].region.nodelist.value[i+j*info_cirroi[num_thread].region.nodelist.height]; k<=info_cirroi[num_thread].region.nodelist.value[i+(j+1)*info_cirroi[num_thread].region.nodelist.height]; k++) {
                            int y = k-(info_cirroi[num_thread].y-info_cirroi[num_thread].radius);
                            if (sign((p1[0]-p0[0])*((double)y-p0[1])-((double)x-p0[0])*(p1[1]-p0[1])) == sign_clear) {
                                // Clear points on this side
                                info_cirroi[num_thread].mask_buffer.value[y+x*info_cirroi[num_thread].mask_buffer.height] = false;
                            }
                        }
                    }
                }

                // Get region again from the mask, then select the
                // biggest region. This is required because subset truncation
                // does not preserve the cirroi's contiguous nature.
                // Setting the last parameter to true for formregion
                // will make the length of the noderange the same as
                // the width of the mask.
                
                // Use a vec_struct_region instead of ncorr_struct_region because it is thread safe
                std::vector<vec_struct_region> region_cirroi_buffer;
                bool removed = false;
                form_regions(region_cirroi_buffer,removed,info_cirroi[num_thread].mask_buffer,0,true);   

                // Its possible for region to be empty. Must check it
                if (region_cirroi_buffer.size() > 0) {      
                    // Find region with the highest number of nodes and keep this region
                    int idx_max = 0;
                    int val_max = region_cirroi_buffer[0].totalpoints;
                    for (int i=1; i<(int)region_cirroi_buffer.size(); i++) {
                        if (region_cirroi_buffer[i].totalpoints > val_max) {
                            // Store
                            idx_max = i;
                            val_max = region_cirroi_buffer[i].totalpoints;
                        }
                    }   

                    // Set mask -----------------------------------------//
                    info_cirroi[num_thread].mask.reset();
                    for (int i=0; i<region_cirroi_buffer[idx_max].height_nodelist; i++) {  
                        for (int j=0; j<region_cirroi_buffer[idx_max].noderange[i]; j+=2) {
                            for (int k=region_cirroi_buffer[idx_max].nodelist[i+j*region_cirroi_buffer[idx_max].height_nodelist]; k<=region_cirroi_buffer[idx_max].nodelist[i+(j+1)*region_cirroi_buffer[idx_max].height_nodelist]; k++) {
                               info_cirroi[num_thread].mask.value[k+i*info_cirroi[num_thread].mask.height] = true;
                            }
                        }
                    }                          

                    // Set region ---------------------------------------//
                    // Convert back to regular coordinates first
                    for (int i=0; i<region_cirroi_buffer[idx_max].height_nodelist; i++) {  
                        for (int j=0; j<region_cirroi_buffer[idx_max].noderange[i]; j++) {
                            region_cirroi_buffer[idx_max].nodelist[i+j*region_cirroi_buffer[idx_max].height_nodelist] += (info_cirroi[num_thread].y-info_cirroi[num_thread].radius);
                        }
                    }
                    region_cirroi_buffer[idx_max].leftbound += (info_cirroi[num_thread].x-info_cirroi[num_thread].radius);
                    region_cirroi_buffer[idx_max].rightbound += (info_cirroi[num_thread].x-info_cirroi[num_thread].radius);
                    region_cirroi_buffer[idx_max].upperbound += (info_cirroi[num_thread].y-info_cirroi[num_thread].radius);
                    region_cirroi_buffer[idx_max].lowerbound += (info_cirroi[num_thread].y-info_cirroi[num_thread].radius);    

                    // Store in info_cirroi[num_thread]
                    for (int i=0; i<region_cirroi_buffer[idx_max].height_nodelist; i++) {  
                        info_cirroi[num_thread].region.noderange.value[i] = region_cirroi_buffer[idx_max].noderange[i];
                        for (int j=0; j<region_cirroi_buffer[idx_max].noderange[i]; j++) {
                            info_cirroi[num_thread].region.nodelist.value[i+j*info_cirroi[num_thread].region.nodelist.height] = region_cirroi_buffer[idx_max].nodelist[i+j*region_cirroi_buffer[idx_max].height_nodelist];
                        }
                    }
                    info_cirroi[num_thread].region.leftbound = region_cirroi_buffer[idx_max].leftbound;
                    info_cirroi[num_thread].region.rightbound = region_cirroi_buffer[idx_max].rightbound;
                    info_cirroi[num_thread].region.upperbound = region_cirroi_buffer[idx_max].upperbound;
                    info_cirroi[num_thread].region.lowerbound = region_cirroi_buffer[idx_max].lowerbound;
                    info_cirroi[num_thread].region.totalpoints = region_cirroi_buffer[idx_max].totalpoints;
                }
            } 
        } 

        // Set output
        cirroi[num_thread].region = info_cirroi[num_thread].region;
        cirroi[num_thread].mask = info_cirroi[num_thread].mask;
        cirroi[num_thread].radius = info_cirroi[num_thread].radius;
        cirroi[num_thread].x = info_cirroi[num_thread].x;
        cirroi[num_thread].y = info_cirroi[num_thread].y;
    } else {
        // Cannot use mexErrMsgTxt because it is not thread safe. mexPrintf is only marginally thread safe.
        mexPrintf("Cirroi has not been set yet for thread %d; this will most likely cause a crash, hopefully this message is seen before matlab closes.\n",num_thread);
    }        
}

// THREADSAFE
bool ncorr_class_roi::withinregion(const int &x,const int &y,const int &num_region) {
    // Returns true if x and y are located in the region specified by num_region
    int idx_roi_x = x-region[num_region].leftbound;
    if (idx_roi_x >= 0 && idx_roi_x < region[num_region].noderange.height) {
        for (int i=0; i < region[num_region].noderange.value[idx_roi_x]; i+=2) {
            if (y >= region[num_region].nodelist.value[idx_roi_x+i*region[num_region].nodelist.height] && y <= region[num_region].nodelist.value[idx_roi_x+(i+1)*region[num_region].nodelist.height]) {
                return true;
            }
        }
    }
    return false;
}

// Inherited Types ------------------------------------------------------//

// NOT THREADSAFE
ncorr_class_inverseregion::ncorr_class_inverseregion(ncorr_class_region &region, const int &border_extrap){
    // Note that border_extrap must be greater than or equal to 1.
    // Form inverseregion given a region input    
    // Allocate space for region - note that upperbound of width will be width of the region + 2
    alloc(region.nodelist.height+2*border_extrap,region.nodelist.width+2);

    // Set bounds - relative to new mat_plot_extrap array
    upperbound = 0;
    lowerbound = region.lowerbound-region.upperbound+2*border_extrap;
    leftbound = 0;
    rightbound = region.rightbound-region.leftbound+2*border_extrap;

    // Find new nodes ---------------------------------------------------//
    // Left nodes: border_extrap space doesnt interact with previous nodes, so just add them
    for (int i=0; i<border_extrap; i++) {
        // Nodelist
        nodelist.value[i] = 0; // Top node
        nodelist.value[i+nodelist.height] = lowerbound; // Bottom Node

        // Noderange
        noderange.value[i] = 2; 

        // totalpoints
        totalpoints += (nodelist.value[i+nodelist.height]-nodelist.value[i])+1;
    }

    // Middle nodes: middle nodes will interact with region nodes
    for (int i=border_extrap; i<nodelist.height-border_extrap; i++) {
        if (region.noderange.value[i-border_extrap] > 1) {
            // First node pair
            nodelist.value[i] = 0;
            nodelist.value[i+nodelist.height] = region.nodelist.value[i-border_extrap]-region.upperbound+border_extrap-1;

            // totalpoints
            totalpoints += (nodelist.value[i+nodelist.height]-nodelist.value[i]+1);

            // Middle node pairs
            for (int j=0; j<region.noderange.value[i-border_extrap]-2; j+=2) {
                // Nodelist
                nodelist.value[i+(j+2)*nodelist.height] = region.nodelist.value[(i-border_extrap)+(j+1)*region.nodelist.height]-region.upperbound+border_extrap+1;
                nodelist.value[i+(j+3)*nodelist.height] = region.nodelist.value[(i-border_extrap)+(j+2)*region.nodelist.height]-region.upperbound+border_extrap-1;

                // totalpoints
                totalpoints += (nodelist.value[i+(j+3)*nodelist.height]-nodelist.value[i+(j+2)*nodelist.height])+1;
            }

            // Last node pair
            nodelist.value[i+region.noderange.value[i-border_extrap]*nodelist.height] = region.nodelist.value[(i-border_extrap)+(region.noderange.value[i-border_extrap]-1)*region.nodelist.height]-region.upperbound+border_extrap+1;
            nodelist.value[i+(region.noderange.value[i-border_extrap]+1)*nodelist.height] = lowerbound;

            // totalpoints
            totalpoints += (nodelist.value[i+(region.noderange.value[i-border_extrap]+1)*nodelist.height]-nodelist.value[i+region.noderange.value[i-border_extrap]*nodelist.height])+1;

            // Noderange
            noderange.value[i] = region.noderange.value[i-border_extrap]+2;
        } else {
            // Nodelist
            nodelist.value[i] = 0;
            nodelist.value[i+nodelist.height] = lowerbound;

            // Noderange
            noderange.value[i] = 2;

            // totalpoints
            totalpoints += (nodelist.value[i+nodelist.height]-nodelist.value[i])+1;
        }
    }

    // Right nodes: border_extrap space doesnt interact with previous nodes, so just add them
    for (int i=nodelist.height-border_extrap; i<nodelist.height; i++) {
        // Nodelist
        nodelist.value[i] = 0; // Top node
        nodelist.value[i+nodelist.height] = lowerbound; // Bottom Node

        // Noderange
        noderange.value[i] = 2;

        // totalpoints
        totalpoints += (nodelist.value[i+nodelist.height]-nodelist.value[i])+1;
    }       
}
