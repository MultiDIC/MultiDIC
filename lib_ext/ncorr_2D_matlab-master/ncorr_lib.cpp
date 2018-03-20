// These are functions used by ncorr. Any function call to the mex API is 
// not threadsafe and is typically noted above the function definition.

#include <mex.h>
#include <math.h>
#include <sstream> 
#include <algorithm> 
#include "standard_datatypes.h"
#include "ncorr_datatypes.h"
#include "ncorr_lib.h"

// ----------------------------------------------------------------------//
// WAITBAR STUFF --------------------------------------------------------//
// ----------------------------------------------------------------------//

// THREAD SAFE 
class_waitbar::class_waitbar() {
    computepoints = 0;
    counter_total = 0;
    counter_sub = 0;
    lhs_waitbar_create[0] = NULL; 
    lhs_getappdata[0] = NULL;       
    rhs_getappdata[0] = NULL;      
    rhs_getappdata[1] = NULL; 
    rhs_waitbar_update[0] = NULL; 
    rhs_waitbar_update[1] = NULL; 
    mat_waitbar_update_1 = NULL;
    mat_getappdata_1 = NULL;
    fraction_waitbar_update = NULL;
}

// NOT THREAD SAFE 
class_waitbar::~class_waitbar() {
    // Close waitbar
    if (lhs_waitbar_create[0] != NULL) {
        mexCallMATLAB(0,NULL,1,lhs_waitbar_create,"delete");
    }        
    
    // Free variables
    if (mat_waitbar_update_1 != NULL) {
        mxDestroyArray(mat_waitbar_update_1);  
    }
    if (mat_getappdata_1 != NULL) {
        mxDestroyArray(mat_getappdata_1);     
    }    
}

// NOT THREAD SAFE 
void class_waitbar::start(const int &num_img,const int &total_imgs,const int &computepoints_buf) {
    // Initialize computepoints
    computepoints = computepoints_buf;
    
    // Create waitbar object
    mxArray *rhs_waitbar_create[4];

    // Create progress string
    std::stringstream stream_progress;
    stream_progress << "Processing image " << num_img+1 << " of " << total_imgs << ".";

    mxArray *mat_wb_0 = mxCreateDoubleMatrix(1, 1, mxREAL); // Make sure to deallocate this - this also initialized to zero
    mxArray *mat_wb_1 = mxCreateString((stream_progress.str()).c_str()); // Make sure to deallocate this
    mxArray *mat_wb_2 = mxCreateString("Name"); // Make sure to deallocate this
    mxArray *mat_wb_3 = mxCreateString("Progress..."); // Make sure to deallocate this

    rhs_waitbar_create[0] = mat_wb_0; 
    rhs_waitbar_create[1] = mat_wb_1; 
    rhs_waitbar_create[2] = mat_wb_2; 
    rhs_waitbar_create[3] = mat_wb_3; 

    // Create Waitbar - handle is stored in lhs_waitbar_create[1]
    mexCallMATLAB(1,lhs_waitbar_create,4,rhs_waitbar_create,"waitbar");
    
    // Set closerequestfcn
    mxArray *rhs_crf[3];    
    mxArray *mat_crf_1 = mxCreateString("CloseRequestFcn"); // Make sure to deallocate this
    mxArray *mat_crf_2 = mxCreateString("setappdata(gcbf,'canceling',true)"); // Make sure to deallocate this
    
    rhs_crf[0] = lhs_waitbar_create[0];
    rhs_crf[1] = mat_crf_1; 
    rhs_crf[2] = mat_crf_2; 
    
    // Set closerequestfcn 
    mexCallMATLAB(0,NULL,3,rhs_crf,"set");
    
    // Set appdata
    mxArray *rhs_setappdata[3];
    
    mxArray *mat_setappdata_1 = mxCreateString("canceling"); // Make sure to deallocate this
    mxArray *mat_setappdata_2 =  mxCreateLogicalMatrix(1,1); // Make sure to deallocate this - this also initialized to false

    rhs_setappdata[0] = lhs_waitbar_create[0];
    rhs_setappdata[1] = mat_setappdata_1;
    rhs_setappdata[2] = mat_setappdata_2;

    // Setappdata
    mexCallMATLAB(0,NULL,3,rhs_setappdata,"setappdata");
    
    // Initialize variables for getappdata
    mat_getappdata_1 = mxCreateString("canceling"); // Make sure to deallocate this
    rhs_getappdata[0] = lhs_waitbar_create[0];
    rhs_getappdata[1] = mat_getappdata_1;
    
    // Form rhs for waitbar update
    mat_waitbar_update_1 =  mxCreateDoubleMatrix(1, 1, mxREAL); // Make sure to deallocate this - this also initialized to zero
    fraction_waitbar_update = mxGetPr(mat_waitbar_update_1);
    rhs_waitbar_update[0] = mat_waitbar_update_1;
    rhs_waitbar_update[1] = lhs_waitbar_create[0];        
    
    // Free initialization variables
    mxDestroyArray(mat_wb_0);
    mxDestroyArray(mat_wb_1);
    mxDestroyArray(mat_wb_2);
    mxDestroyArray(mat_wb_3);
    mxDestroyArray(mat_crf_1);
    mxDestroyArray(mat_crf_2);
    mxDestroyArray(mat_setappdata_1);
    mxDestroyArray(mat_setappdata_2);
}

// THREAD SAFE 
void class_waitbar::increment() {
    ++counter_sub;
    ++counter_total; 
}

// NOT THREAD SAFE 
bool class_waitbar::updateandcheck() {
    // Do not update waitbar every iteration- instead, do it only when counter_sub exceeds WAITBAR_UPDATE
    if (counter_sub > WAITBAR_UPDATE) {
        // Reset counter_sub to zero
        counter_sub = 0;
            
        // See if waitbar has been cancelled
        if (!getappdata()) {
            // Update waitbar - note that fraction_waitbar_update is the pr for rhs_waitbar_update
            *fraction_waitbar_update = std::min(1.0,(double)counter_total/(double)computepoints);
            mexCallMATLAB(0,NULL,2,rhs_waitbar_update,"waitbar");
        } else {
            return false;
        }
    }   
    return true;    
}

// NOT THREAD SAFE 
bool class_waitbar::getappdata() {
    // This function is slow, so only call it every so often
    // Used to see if waitbar has been cancelled
    mexCallMATLAB(1,lhs_getappdata,2,rhs_getappdata,"getappdata");
    
    if (*(mxGetLogicals(lhs_getappdata[0]))) {
        return true;  // If cancelled, return true
    } else {
        return false; // Else, return false and continue with analysis
    }
}

// ----------------------------------------------------------------------//
// FUNCTIONS ------------------------------------------------------------//
// ----------------------------------------------------------------------//

// THREAD SAFE 
double ncorr_round(const double &r) {
// Older compilers dont have round, but newer ones do, so I've changed the 
// name to ncorr_round to prevent a redefinition
    return (r>0.0) ? floor(r+0.5) : ceil(r-0.5);
}

// THREAD SAFE 
int sign(const double &r) {
    if (r > 0) {
        return 1;
    } else if (r < 0) {
        return -1;
    } else {
        // Possibly alter this function to use LAMDA and fabs to test if value is zero
        return 0;
    }
}

// THREAD SAFE 
int mod_pos(const int &i,const int &n) {
    return (i % n + n) % n;
}

// THREAD SAFE 
void form_boundary(std::vector<std::vector<int> > &vec_boundary,const std::vector<int> &point_init,const class_logical_array &mask,int &direc) {        
    // Initialize newpoint
    std::vector<int> newpoint(2,0);
    newpoint[0] = point_init[0];
    newpoint[1] = point_init[1];
    
    // Push newpoint into vec_boundary
    vec_boundary.push_back(newpoint);
    newpoint[0] = -1; // Set this to -1 (an impossible value to initialize while condition
    newpoint[1] = -1;    
    while (newpoint[0] != vec_boundary[0][0] || newpoint[1] != vec_boundary[0][1]) { // Exit when newpoint equals first boundary point - i.e. when the boundary is closed
        // Get new direction first
        direc = (direc+6) % 8;
        
        // Now search
        int i;
        for (i=0; i<9; i++) {
            // Convert direction into increment
            int x_inc;
            int y_inc;
            if (direc == 0) {
                x_inc = 1;
                y_inc = 0;
            } else if (direc == 1) {
                x_inc = 1;
                y_inc = -1;
            } else if (direc == 2) {
                x_inc = 0;
                y_inc = -1;
            } else if (direc == 3) {
                x_inc = -1;
                y_inc = -1;
            } else if (direc == 4) {
                x_inc = -1;
                y_inc = 0;
            } else if (direc == 5) {
                x_inc = -1;
                y_inc = 1;
            } else if (direc == 6) {
                x_inc = 0;
                y_inc = 1;
            } else {
                x_inc = 1;
                y_inc = 1;
            }
            
            // Set newpoint
            newpoint[0] = vec_boundary.back()[0]+x_inc;
            newpoint[1] = vec_boundary.back()[1]+y_inc;
            if (newpoint[0] >= 0 && newpoint[1] >= 0 &&
                newpoint[0] < mask.width && newpoint[1] < mask.height && 
                mask.value[newpoint[1]+ newpoint[0]*mask.height]) {
                vec_boundary.push_back(newpoint);
                break;
            } else {
                direc = (direc+1) % 8;
            }                               
        } 
        
        // This means the direction went past opposite direction of initial 
        // direction - this means this is just a single point, so break;
        if (i == 9) {
            break;
        }
    }    
}

// Structures for form_regions function ----------------------------//

struct local_struct_roi_overall {    
    // Constructor
    local_struct_roi_overall() {    // THREAD SAFE 
        leftbound = 0;
        rightbound = 0;        
    }    
    
    // Properties
    std::vector<std::vector<int> > nodelist;
    std::vector<std::vector<bool> > activelines;
    int leftbound;
    int rightbound;
};

struct local_struct_roi_separate {
public:
    // Constructor
    local_struct_roi_separate() {    // THREAD SAFE 
        leftbound = 0;
        rightbound = 0;
        totalpoints = 0;
    }        
    
    // Properties
    std::vector<std::vector<int> > nodelist;
    int leftbound;
    int rightbound;
    int totalpoints;
};

struct local_struct_queue_roi_separate {
public:
    // Constructor
    local_struct_queue_roi_separate() {    // THREAD SAFE 
        node_top = 0;
        node_bottom = 0;
        index = 0;
    }    
    
    // Properties
    int node_top;
    int node_bottom;
    int index;
};

// THREAD SAFE 
vec_struct_region::vec_struct_region() {
    height_nodelist = 0;
    width_nodelist = 0;
    upperbound = 0;
    lowerbound = 0;
    leftbound = 0;
    rightbound = 0;
    totalpoints = 0;
}

// THREAD SAFE 
void form_regions(std::vector<vec_struct_region> &region,bool &removed,const class_logical_array &mask,const int &cutoff,const bool &preservelength) {
    // Initialize removed to false
    removed = false;
    
    // Form overall ROI -----------------------------//    
    local_struct_roi_overall roi_overall;
    
    // Find bounds - if preservelength is set, then set the left and right bounds to the edges of the mask    
    bool firstpoint = true;
    if (preservelength) {
        roi_overall.leftbound = 0;
        roi_overall.rightbound = mask.width-1;
        for (int i=0; i<mask.width; i++) {
            for (int j=0; j<mask.height; j++) {
                if (mask.value[j+i*mask.height]) {
                    firstpoint = false;
                    break;
                }
            }
            if (!firstpoint) {
                break;
            }
        }
    } else {
        roi_overall.leftbound = mask.width-1;
        roi_overall.rightbound = 0;
        for (int i=0; i<mask.width; i++) {
            for (int j=0; j<mask.height; j++) {
                if (firstpoint && mask.value[j+i*mask.height]) {
                    // This is the first element encountered - set this to the leftbound
                    roi_overall.leftbound = i;
                    firstpoint = false;
                }
                if (mask.value[j+i*mask.height]) {
                    // This is the next update for the rightbound
                    roi_overall.rightbound = i;
                    break; // Skip to next line
                }
            }
        }
    }

    // See if algorithm is still trying to find first point
    if (firstpoint) {
        // Mask is empty, return.
        return;
    }
    
    // Get node pairs - form nodelist and activelines based on bounds
    roi_overall.nodelist.resize(roi_overall.rightbound-roi_overall.leftbound+1);
    roi_overall.activelines.resize(roi_overall.rightbound-roi_overall.leftbound+1);
    for (int i=roi_overall.leftbound; i<=roi_overall.rightbound; i++) {
        bool start = false;
        int node_top;
        int node_bottom;
        for (int j=0; j<mask.height; j++) { 
            if (!start && mask.value[j+i*mask.height]) {
                // This finds the top node
                start = true;
                node_top = j;
            }

            if (start && (!mask.value[j+i*mask.height] || j == mask.height-1)) {
                // This finds the bottom node
                start = false;
                if (j == mask.height-1 && mask.value[j+i*mask.height]) {
                    node_bottom = j; // Dont subtract one because edge of image has been reached
                } else {
                    node_bottom = j-1; // Subtract one to step back
                }
                // Now add node pair to roi_overall
                roi_overall.nodelist[i-roi_overall.leftbound].push_back(node_top); 
                roi_overall.nodelist[i-roi_overall.leftbound].push_back(node_bottom); 
                roi_overall.activelines[i-roi_overall.leftbound].push_back(true); 
            }
        }
    }
    
    // Separate Regions -----------------------------//   
    // Regions are made 4-way contiguous
    std::vector<local_struct_roi_separate> roi_separate;
        
    // Form queue
    std::vector<local_struct_queue_roi_separate> queue_roi_separate;    
    // Scan over nodelist - keep track of columns with col
    int col = -1;
    while (col < (int)roi_overall.nodelist.size()-1) {
        // Increment column every iteration
        ++col;
        
        if (roi_overall.nodelist[col].size() == 0) {
            // If nodelist is empty then continue
            continue;
        } else {
            int idx_node = 0;
            bool activenodes = false;
            for (int j=0; j<(int)roi_overall.activelines[col].size(); j++) {
                if (roi_overall.activelines[col][j]) {
                    idx_node = j*2;
                    roi_overall.activelines[col][j] = false; // Inactivate node pair
                    activenodes = true;
                    break;
                }
            }
            if (!activenodes) {
                // No active nodes are left - continue to next scanline
                continue;
            }

            // Initialize roi_separate_buffer
            local_struct_roi_separate roi_separate_buffer; 
            roi_separate_buffer.nodelist.resize(roi_overall.rightbound-roi_overall.leftbound+1);
            roi_separate_buffer.totalpoints = 0;
            // For bounds, set them according to preservelength
            if (preservelength) {
                roi_separate_buffer.leftbound = roi_overall.leftbound;
                roi_separate_buffer.rightbound = roi_overall.rightbound;
            } else {
                // Set left bound here; right bound will get updated during computation
                roi_separate_buffer.leftbound = col+roi_overall.leftbound;
                roi_separate_buffer.rightbound = 0;
            }

            // Initialize queue_roi_separate              
            local_struct_queue_roi_separate queue_template;
            queue_template.node_top = roi_overall.nodelist[col][idx_node];
            queue_template.node_bottom = roi_overall.nodelist[col][idx_node+1];
            queue_template.index = col;
            queue_roi_separate.push_back(queue_template);

            // Enter while loop - Exit when queue_roi_separate is empty
            while (queue_roi_separate.size() > 0) {
                // Take node pair out of queue_roi_separate and compare it to node pairs left and right
                local_struct_queue_roi_separate queue_roi_separate_buffer = queue_roi_separate.back(); // Load nodes
                queue_roi_separate.pop_back(); // Delete nodes

                // Insert node pairs and then sort
                roi_separate_buffer.nodelist[queue_roi_separate_buffer.index].push_back(queue_roi_separate_buffer.node_top);
                roi_separate_buffer.nodelist[queue_roi_separate_buffer.index].push_back(queue_roi_separate_buffer.node_bottom);
                sort(roi_separate_buffer.nodelist[queue_roi_separate_buffer.index].begin(),roi_separate_buffer.nodelist[queue_roi_separate_buffer.index].end());

                // Update totalpoints
                roi_separate_buffer.totalpoints += queue_roi_separate_buffer.node_bottom-queue_roi_separate_buffer.node_top+1;

                // Compare to node pairs LEFT. Any node pairs which interact are added to the queue_roi_separate
                if (queue_roi_separate_buffer.index > 0) {
                    for (int j=0; j<(int)roi_overall.nodelist[queue_roi_separate_buffer.index-1].size(); j+=2) {
                        if (roi_overall.nodelist[queue_roi_separate_buffer.index-1][j] > queue_roi_separate_buffer.node_bottom) {
                            // top node comes after bottom node of buffer
                            break;
                        } else if (!(roi_overall.activelines[queue_roi_separate_buffer.index-1][j/2])) {
                            continue;
                        } else if (roi_overall.nodelist[queue_roi_separate_buffer.index-1][j+1] < queue_roi_separate_buffer.node_top) {
                            // bottom node comes before top node of buffer
                            continue;
                        } else {
                            // Nodes interact
                            queue_template.node_top = roi_overall.nodelist[queue_roi_separate_buffer.index-1][j];
                            queue_template.node_bottom = roi_overall.nodelist[queue_roi_separate_buffer.index-1][j+1];
                            queue_template.index = queue_roi_separate_buffer.index-1;
                            queue_roi_separate.push_back(queue_template);
                            roi_overall.activelines[queue_roi_separate_buffer.index-1][j/2] = 0;
                        } 
                    }
                }

                // Compare to node pairs RIGHT. Any node pairs which interact are added to the queue_roi_separate
                if (queue_roi_separate_buffer.index < (int)roi_overall.nodelist.size()-1) {
                    for (int j=0; j<(int)roi_overall.nodelist[queue_roi_separate_buffer.index+1].size(); j+=2) {
                        if (roi_overall.nodelist[queue_roi_separate_buffer.index+1][j] > queue_roi_separate_buffer.node_bottom) {
                            // top node comes after bottom node of buffer
                            break;
                        } else if (!(roi_overall.activelines[queue_roi_separate_buffer.index+1][j/2])) {
                            continue;
                        } else if (roi_overall.nodelist[queue_roi_separate_buffer.index+1][j+1] < queue_roi_separate_buffer.node_top) {
                            // bottom node comes before top node of buffer
                            continue;
                        } else {
                            // Nodes interact
                            queue_template.node_top = roi_overall.nodelist[queue_roi_separate_buffer.index+1][j];
                            queue_template.node_bottom = roi_overall.nodelist[queue_roi_separate_buffer.index+1][j+1];
                            queue_template.index = queue_roi_separate_buffer.index+1;
                            queue_roi_separate.push_back(queue_template);
                            roi_overall.activelines[queue_roi_separate_buffer.index+1][j/2] = 0;
                        }
                    }
                }

                // Update rightbound if preservelength is disabled
                if (!preservelength && queue_roi_separate_buffer.index > roi_separate_buffer.rightbound) {
                    roi_separate_buffer.rightbound = queue_roi_separate_buffer.index;
                }
            }
            
            // Finalize rightbound if preservelength is disabled
            if (!preservelength) {
                // Add roi_overall.leftbound to roi_separate_buffer.rightbound to account for offset
                roi_separate_buffer.rightbound += roi_overall.leftbound;
            }
            
            // Subtract one from col in order to recheck the scan line to ensure all nodes are gone before proceeding
            --col;

            // Insert roi_separate buffer into roi_separate if the totalpoints is above the cutoff
            if (roi_separate_buffer.totalpoints > cutoff) {      
                roi_separate.push_back(roi_separate_buffer);
            } else {
                removed = true;
            }
        } 
    }

    // Finish regions -------------------------------//
    region.resize(roi_separate.size());
    // Set data - cycle over separate ROIs
    for (int i=0; i<(int)roi_separate.size(); i++) {
        // Find the max number of nodes to determine what width to set nodelist to
        int maxnodes = 0; 
        for (int j=0; j<(int)roi_separate[i].nodelist.size(); j++) {
            if ((int)roi_separate[i].nodelist[j].size() > maxnodes) {
                maxnodes = (int)roi_separate[i].nodelist[j].size();
            }
        }

        // Allocate memory
        region[i].height_nodelist = roi_separate[i].rightbound-roi_separate[i].leftbound+1;
        region[i].width_nodelist = maxnodes;
        region[i].nodelist.resize(region[i].height_nodelist*region[i].width_nodelist);
        region[i].noderange.resize(region[i].height_nodelist);

        // Set bounds
        region[i].leftbound = roi_separate[i].leftbound;
        region[i].rightbound = roi_separate[i].rightbound;
        region[i].totalpoints = roi_separate[i].totalpoints;
        
        // Initialize upper and lower bounds and then calculate them
        region[i].upperbound = mask.height;
        region[i].lowerbound = 0;

        // Set noderange and nodelist
        int idx_x = 0; // Use idx_x for index since roi_separate's nodelist size is based on roi_overall's
        for (int j=0; j<(int)roi_separate[i].nodelist.size(); j++) {
            if (preservelength || roi_separate[i].nodelist[j].size()>0) {
                // Noderange:
                region[i].noderange[idx_x] = (int)roi_separate[i].nodelist[j].size();
                
                // Nodelist and upper and lower bounds:
                for (int k=0; k<(int)roi_separate[i].nodelist[j].size(); k+=2) {
                    region[i].nodelist[idx_x+k*region[i].height_nodelist] = roi_separate[i].nodelist[j][k];
                    region[i].nodelist[idx_x+(k+1)*region[i].height_nodelist] = roi_separate[i].nodelist[j][k+1];

                    if (roi_separate[i].nodelist[j][k] < region[i].upperbound) {
                        region[i].upperbound = roi_separate[i].nodelist[j][k];
                    }

                    if (roi_separate[i].nodelist[j][k+1] > region[i].lowerbound) {
                        region[i].lowerbound = roi_separate[i].nodelist[j][k+1];
                    }                
                }
                ++idx_x;
            }
        }
    }
}

// THREAD SAFE 
void form_union(std::vector<vec_struct_region> &region_union,const std::vector<ncorr_class_region> &region,const class_logical_array &mask,const bool &inplace) {    
// If this is done in-place, then region_union and region are the same. If not, region_union is empty and needs to be allocated space
    // Cycle over regions
    for (int i=0; i<(int)region.size(); i++) {
        // If this is not done in-place, you need to figure out what size the region should be and then allocate space for it
        if (!inplace) {
            int maxnodes = 0; 
            for (int j=0; j<region[i].noderange.height; j++) {
                int x = j+region[i].leftbound;
                int submaxnodes = 0;
                for (int k=0; k<region[i].noderange.value[j]; k+=2) {
                    bool start = false;                    
                    for (int m=region[i].nodelist.value[j+k*region[i].nodelist.height]; m<=region[i].nodelist.value[j+(k+1)*region[i].nodelist.height]; m++) {
                        int y = m;
                        // Max nodes
                        if (!start && mask.value[y+x*mask.height]) {
                            start = true;
                        }
                        if (start && (!mask.value[y+x*mask.height] || y == region[i].nodelist.value[j+(k+1)*region[i].nodelist.height])) {
                            start = false;
                            // End of node pair, add two
                            submaxnodes += 2;
                        }              
                    }
                }
                if (submaxnodes > maxnodes) {
                    maxnodes = submaxnodes;
                }
            }
            
            // Form fields
            region_union[i].height_nodelist = region[i].nodelist.height;
            region_union[i].width_nodelist = std::max(maxnodes,2);
            region_union[i].nodelist.resize(region_union[i].height_nodelist*region_union[i].width_nodelist);
            region_union[i].noderange.resize(region_union[i].height_nodelist);
            
            // Copy Bounds
            region_union[i].leftbound = region[i].leftbound;
            region_union[i].rightbound = region[i].rightbound;
            region_union[i].upperbound = region[i].upperbound;
            region_union[i].lowerbound = region[i].lowerbound;
        } else {
            // Put totalpoints to zero for the in-place case
            region_union[i].totalpoints = 0; 
        }
        
        // Allocate buffers to store nodelist
        std::vector<int> nodelist_buffer;
        nodelist_buffer.resize(region_union[i].width_nodelist);
        int noderange_buffer = 0;

        // Get union
        for (int j=0; j<region[i].noderange.height; j++) {
            int x = j+region[i].leftbound;

            // Find updated nodelist and noderange for this column
            noderange_buffer = 0;
            bool start = false;
            for (int k=0; k<region[i].noderange.value[j]; k+=2) {                                           
                int node_top = 0;
                int node_bottom = 0;
                for (int m=region[i].nodelist.value[j+k*region[i].nodelist.height]; m<=region[i].nodelist.value[j+(k+1)*region[i].nodelist.height]; m++) {
                    int y = m;
                    if (!start && mask.value[y+x*mask.height]) {
                        start = true;
                        node_top = y;
                    }
                    if (start && (!mask.value[y+x*mask.height] || y == region[i].nodelist.value[j+(k+1)*region[i].nodelist.height])) {
                        start = false;
                        // End of node pair
                        if (y == region[i].nodelist.value[j+(k+1)*region[i].nodelist.height] && mask.value[y+x*mask.height]) {
                            node_bottom = y;
                        } else {
                            node_bottom = y-1;
                        }

                        // Store nodepair in buffer
                        nodelist_buffer[noderange_buffer] = node_top;
                        nodelist_buffer[noderange_buffer+1] = node_bottom;
                                                
                        // Update noderange in buffer
                        noderange_buffer += 2;
                        
                        // Update totalpoints
                        region_union[i].totalpoints += node_bottom - node_top + 1;
                    }              
                }
            }     

            // Store buffers
            for (int k=0; k<noderange_buffer; k++) {
                region_union[i].nodelist[j+k*region_union[i].height_nodelist] = nodelist_buffer[k];                                        
            }
            // Store noderange_buffer
            region_union[i].noderange[j] = noderange_buffer;
        }
    }
}

// THREAD SAFE 
void cholesky(std::vector<double> &mat,bool &positivedef,const int &size_mat) {
    // Stores cholesky decomposition in lower triangular portion of matrix in-place
    for (int i=0; i<size_mat; i++) {
        if (i>0) {
            for (int j=size_mat-1; j>=i; j--) {
                double mat_buffer = 0.0;
                for (int k=0; k<i; k++) {
                    mat_buffer += mat[j+k*size_mat]*mat[i+k*size_mat];
                }
                mat[j+i*size_mat] -= mat_buffer;
            }
        }
        if (mat[i+i*size_mat] > LAMBDA) { 
            double diag_sqrt = sqrt(mat[i+i*size_mat]);
            for (int j=i; j<size_mat; j++) {
                mat[j+i*size_mat] /= diag_sqrt;
            } 
        } else {
            positivedef = false;
            break;
        }
    }
}

// THREAD SAFE 
void forwardsub(std::vector<double> &vec,const std::vector<double> &mat,const int &size_mat) {
    // Assumes mat is lower triangular, stores solution in place
    vec[0] /= mat[0];
    for (int i=1; i<size_mat; i++) {
        double vec_buffer = 0.0;
        for (int j=0; j<i; j++) {
            vec_buffer += mat[i+j*size_mat]*vec[j];
        }
        vec[i] = (vec[i]-vec_buffer)/mat[i+i*size_mat];
    }
}

// THREAD SAFE 
void backwardsub(std::vector<double> &vec,const std::vector<double> &mat,const int &size_mat) {
    // Backward substitution is used for upper triangular matrices. 
    // However, this function is used with cholesky decomposition which 
    // stores the cholesky matrix in the lower triangular portion of the matrix.
    // So, this function assumes the upper-triangular matrix input is transposed
    // and stored in the lower triangular portion of mat. Lastly, the solution
    // is stored in-place.
    vec[size_mat-1] /= mat[(size_mat-1)+(size_mat-1)*size_mat];
    for (int i=size_mat-2; i>-1; i--) {
        double vec_buffer = 0.0;
        for (int j=i+1; j<size_mat; j++) {
            // switch indices because only lower triangular portion is G
            vec_buffer += mat[j+i*size_mat]*vec[j]; 
        }
        vec[i] = (vec[i]-vec_buffer)/mat[i+i*size_mat];
    }
}

// NOT THREAD SAFE 
OUT interp_qbs(double &interp,const double &x_tilda,const double &y_tilda,const class_double_array &plot_interp,const class_logical_array &mask,const int &offset_x,const int &offset_y,const int &border_bcoef) {
// This function performs biquintic b-spline interpolation. It makes sure that the
// rounded desired interpolation position (x_tilda and y_tilda) is located within 
// the mask before interpolating. If the position is successfully interpolated, 
// this function returns success;

// The b-spline coefficient plots are generally a subset of the entire image and 
// therefore require offsets as inputs. They usually also involve a border which 
// is used to expand the plots before obtaining the b-spline coefficients to mitigate 
// ringing.

    // Static memory buffers
    // Note that static memory will be shared between threads, so don't use this if multithreading is present
    static double x_vec_buffer[6];
    static double y_vec_buffer[6];
    static double QK_B_QKT_buffer[36];

    int x_tilda_floor = (int)floor(x_tilda);
    int y_tilda_floor = (int)floor(y_tilda);
    
    int x_tilda_round = (int)ncorr_round(x_tilda);
    int y_tilda_round = (int)ncorr_round(y_tilda);
        
    // Get bounds of the desired b-spline coefficients used for interpolation
    int top = y_tilda_floor-offset_y+border_bcoef-2;
    int left = x_tilda_floor-offset_x+border_bcoef-2;
    int bottom = y_tilda_floor-offset_y+border_bcoef+3;
    int right = x_tilda_floor-offset_x+border_bcoef+3;
        
    // Make sure top, left, bottom, and right are within the b-spline 
    // coefficient array. top, left, bottom and right are the bounding 
    // box of the b-spline coefficients used for interpolation of this
    // point;
    if (top >= 0 &&
        left >= 0 && 
        bottom < plot_interp.height &&
        right < plot_interp.width &&
        x_tilda_round >= 0 && x_tilda_round < mask.width &&
        y_tilda_round >= 0 && y_tilda_round < mask.height &&
        mask.value[y_tilda_round + x_tilda_round*mask.height]) {   
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

        // Compute QK*B_coef*QK^T 
        double b0 = plot_interp.value[top+left*plot_interp.height];
        double b1 = plot_interp.value[top+1+left*plot_interp.height];
        double b2 = plot_interp.value[top+2+left*plot_interp.height];
        double b3 = plot_interp.value[top+3+left*plot_interp.height];
        double b4 = plot_interp.value[top+4+left*plot_interp.height];
        double b5 = plot_interp.value[top+5+left*plot_interp.height];
        double b6 = plot_interp.value[top+(left+1)*plot_interp.height];
        double b7 = plot_interp.value[top+1+(left+1)*plot_interp.height];
        double b8 = plot_interp.value[top+2+(left+1)*plot_interp.height];
        double b9 = plot_interp.value[top+3+(left+1)*plot_interp.height];
        double b10 = plot_interp.value[top+4+(left+1)*plot_interp.height];
        double b11 = plot_interp.value[top+5+(left+1)*plot_interp.height];
        double b12 = plot_interp.value[top+(left+2)*plot_interp.height];
        double b13 = plot_interp.value[top+1+(left+2)*plot_interp.height];
        double b14 = plot_interp.value[top+2+(left+2)*plot_interp.height];
        double b15 = plot_interp.value[top+3+(left+2)*plot_interp.height];
        double b16 = plot_interp.value[top+4+(left+2)*plot_interp.height];
        double b17 = plot_interp.value[top+5+(left+2)*plot_interp.height];
        double b18 = plot_interp.value[top+(left+3)*plot_interp.height];
        double b19 = plot_interp.value[top+1+(left+3)*plot_interp.height];
        double b20 = plot_interp.value[top+2+(left+3)*plot_interp.height];
        double b21 = plot_interp.value[top+3+(left+3)*plot_interp.height];
        double b22 = plot_interp.value[top+4+(left+3)*plot_interp.height];
        double b23 = plot_interp.value[top+5+(left+3)*plot_interp.height];
        double b24 = plot_interp.value[top+(left+4)*plot_interp.height];
        double b25 = plot_interp.value[top+1+(left+4)*plot_interp.height];
        double b26 = plot_interp.value[top+2+(left+4)*plot_interp.height];
        double b27 = plot_interp.value[top+3+(left+4)*plot_interp.height];
        double b28 = plot_interp.value[top+4+(left+4)*plot_interp.height];
        double b29 = plot_interp.value[top+5+(left+4)*plot_interp.height];
        double b30 = plot_interp.value[top+(left+5)*plot_interp.height];
        double b31 = plot_interp.value[top+1+(left+5)*plot_interp.height];
        double b32 = plot_interp.value[top+2+(left+5)*plot_interp.height];
        double b33 = plot_interp.value[top+3+(left+5)*plot_interp.height];
        double b34 = plot_interp.value[top+4+(left+5)*plot_interp.height];
        double b35 = plot_interp.value[top+5+(left+5)*plot_interp.height];

        // Compute QK*B*QK^T
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

        // Interpolate value
        interp = (x_vec_buffer[0]*QK_B_QKT_buffer[0]+x_vec_buffer[1]*QK_B_QKT_buffer[6]+x_vec_buffer[2]*QK_B_QKT_buffer[12]+x_vec_buffer[3]*QK_B_QKT_buffer[18]+x_vec_buffer[4]*QK_B_QKT_buffer[24]+x_vec_buffer[5]*QK_B_QKT_buffer[30])*y_vec_buffer[0]+
                 (x_vec_buffer[0]*QK_B_QKT_buffer[1]+x_vec_buffer[1]*QK_B_QKT_buffer[7]+x_vec_buffer[2]*QK_B_QKT_buffer[13]+x_vec_buffer[3]*QK_B_QKT_buffer[19]+x_vec_buffer[4]*QK_B_QKT_buffer[25]+x_vec_buffer[5]*QK_B_QKT_buffer[31])*y_vec_buffer[1]+
                 (x_vec_buffer[0]*QK_B_QKT_buffer[2]+x_vec_buffer[1]*QK_B_QKT_buffer[8]+x_vec_buffer[2]*QK_B_QKT_buffer[14]+x_vec_buffer[3]*QK_B_QKT_buffer[20]+x_vec_buffer[4]*QK_B_QKT_buffer[26]+x_vec_buffer[5]*QK_B_QKT_buffer[32])*y_vec_buffer[2]+
                 (x_vec_buffer[0]*QK_B_QKT_buffer[3]+x_vec_buffer[1]*QK_B_QKT_buffer[9]+x_vec_buffer[2]*QK_B_QKT_buffer[15]+x_vec_buffer[3]*QK_B_QKT_buffer[21]+x_vec_buffer[4]*QK_B_QKT_buffer[27]+x_vec_buffer[5]*QK_B_QKT_buffer[33])*y_vec_buffer[3]+
                 (x_vec_buffer[0]*QK_B_QKT_buffer[4]+x_vec_buffer[1]*QK_B_QKT_buffer[10]+x_vec_buffer[2]*QK_B_QKT_buffer[16]+x_vec_buffer[3]*QK_B_QKT_buffer[22]+x_vec_buffer[4]*QK_B_QKT_buffer[28]+x_vec_buffer[5]*QK_B_QKT_buffer[34])*y_vec_buffer[4]+
                 (x_vec_buffer[0]*QK_B_QKT_buffer[5]+x_vec_buffer[1]*QK_B_QKT_buffer[11]+x_vec_buffer[2]*QK_B_QKT_buffer[17]+x_vec_buffer[3]*QK_B_QKT_buffer[23]+x_vec_buffer[4]*QK_B_QKT_buffer[29]+x_vec_buffer[5]*QK_B_QKT_buffer[35])*y_vec_buffer[5];
        return SUCCESS;
    }    
    return FAILED;
}

// THREADSAFE    
void expand_filt(class_double_array &plot_extrap,const ncorr_class_inverseregion &inverseregion) {
// NOTE that this requires input arrays to be 2x2 or greater in size since average filter uses
// an average of the four adjacent pixels without bounds checking.

    // Create mask and paint buffers
    std::vector<bool> mask_buffer_new;
    std::vector<bool> mask_buffer_old;
    std::vector<double> paint_buffer;

    mask_buffer_new.resize(plot_extrap.height*plot_extrap.width);
    mask_buffer_old.resize(plot_extrap.height*plot_extrap.width);
    paint_buffer.resize(plot_extrap.height*plot_extrap.width);     
    
    // Initialize mask buffers
    for (int i=0; i<inverseregion.noderange.height; i++) {
        int x = i;
        for (int j=0; j<inverseregion.noderange.value[i]; j+=2) {
            for (int k=inverseregion.nodelist.value[i+j*inverseregion.nodelist.height]; k<=inverseregion.nodelist.value[i+(j+1)*inverseregion.nodelist.height]; k++) {
                int y = k;
                mask_buffer_new[y+x*plot_extrap.height] = true;
                mask_buffer_old[y+x*plot_extrap.height] = true;
            }
        }
    }

    // Expand data ----------------------------------------------//
    int totalcounter = 1; // Keeps track of how many points were added by the border operation; if this value is zero then the image is full
    while (totalcounter > 0) {
        totalcounter = 0;
        // Expand by 1 pixel border
        for (int i=0; i<inverseregion.noderange.height; i++) {
            int x = i;
            for (int j=0; j<inverseregion.noderange.value[i]; j+=2) {
                for (int k=inverseregion.nodelist.value[i+j*inverseregion.nodelist.height]; k<=inverseregion.nodelist.value[i+(j+1)*inverseregion.nodelist.height]; k++) {
                    int y = k;
                    if (mask_buffer_old[y + x*plot_extrap.height]) {
                        // If value is inactive, search four neighbors and take their average value and replace this point with that value
                        double sum = 0.0;
                        int counter = 0;
                        
                        // Up
                        if (y > 0) {
                            if (!mask_buffer_old[y-1 + x*plot_extrap.height]) {
                                sum += plot_extrap.value[y-1 + x*plot_extrap.height];
                                counter++;
                                totalcounter++;
                            }
                        }
                        // Down
                        if (y < plot_extrap.height-1) {
                            if (!mask_buffer_old[y+1 + x*plot_extrap.height]) {
                                sum += plot_extrap.value[y+1 + x*plot_extrap.height];
                                counter++;
                                totalcounter++;
                            }
                        }
                        // Left
                        if (x > 0) {
                            if (!mask_buffer_old[y + (x-1)*plot_extrap.height]) {
                                sum += plot_extrap.value[y + (x-1)*plot_extrap.height];
                                counter++;
                                totalcounter++;
                            }
                        }
                        // Right
                        if (x < plot_extrap.width-1) {
                            if (!mask_buffer_old[y + (x+1)*plot_extrap.height]) {
                                sum += plot_extrap.value[y + (x+1)*plot_extrap.height];
                                counter++;
                                totalcounter++;
                            }
                        }
                        
                        // If counter is zero then no neighbors had any values, leave this point active and continue
                        if (counter > 0) {
                            plot_extrap.value[y + x*plot_extrap.height] = sum/(double)counter;
                            mask_buffer_new[y + x*plot_extrap.height] = false; // Inactivate point
                        }
                    }                
                }
            }
        }
        
        // Update mask_buffer_old with mask_buffer_new
        for (int i=0; i<inverseregion.noderange.height; i++) {
            int x = i;
            for (int j=0; j<inverseregion.noderange.value[i]; j+=2) {
                for (int k=inverseregion.nodelist.value[i+j*inverseregion.nodelist.height]; k<=inverseregion.nodelist.value[i+(j+1)*inverseregion.nodelist.height]; k++) {
                    int y = k;
                    mask_buffer_old[y + x*plot_extrap.height] = mask_buffer_new[y + x*plot_extrap.height];
                }
            }
        }
    }

    // Filter buffer_interdata_plot------------------------------//
    // Number of filter iterations will be based purely on a heuristic
    // Initialize paint buffer
    for (int i=0; i<inverseregion.noderange.height; i++) {
        int x = i;
        for (int j=0; j<inverseregion.noderange.value[i]; j+=2) {
            for (int k=inverseregion.nodelist.value[i+j*inverseregion.nodelist.height]; k<=inverseregion.nodelist.value[i+(j+1)*inverseregion.nodelist.height]; k++) {
                int y = k;
                paint_buffer[y+x*plot_extrap.height] = plot_extrap.value[y+x*plot_extrap.height];
            }
        }
    }
    
    // Start average filter - note that this requires at least a 2x2 array
    for (int i=0; i<FILTERITERATIONS; i++) {        
        for (int j=0; j<inverseregion.noderange.height; j++) {
            int x = j;
            if (x == 0) { // Left Column                
                for (int k=0; k<inverseregion.noderange.value[j]; k+=2) {
                    for (int l=inverseregion.nodelist.value[j+k*inverseregion.nodelist.height]; l<=inverseregion.nodelist.value[j+(k+1)*inverseregion.nodelist.height]; l++) {
                        int y = l;
                        if (y == 0) { //Top-left Corner
                            paint_buffer[y+x*plot_extrap.height] = (plot_extrap.value[y+x*plot_extrap.height]+plot_extrap.value[y+(x+1)*plot_extrap.height]+plot_extrap.value[y+1+x*plot_extrap.height]+plot_extrap.value[y+1+(x+1)*plot_extrap.height])/4;
                        } else if (y == plot_extrap.height-1) { //Bottom-left Corner
                            paint_buffer[y+x*plot_extrap.height] = (plot_extrap.value[y-1+x*plot_extrap.height]+plot_extrap.value[y-1+(x+1)*plot_extrap.height]+plot_extrap.value[y+x*plot_extrap.height]+plot_extrap.value[y+(x+1)*plot_extrap.height])/4;
                        } else { // Left side
                            paint_buffer[y+x*plot_extrap.height] = (plot_extrap.value[y-1+x*plot_extrap.height]+plot_extrap.value[y-1+(x+1)*plot_extrap.height]+plot_extrap.value[y+x*plot_extrap.height]+plot_extrap.value[y+(x+1)*plot_extrap.height]+plot_extrap.value[y+1+x*plot_extrap.height]+plot_extrap.value[y+1+(x+1)*plot_extrap.height])/6;
                        }
                    }
                }            
            } else if (x == plot_extrap.width-1) { // Right Column
                for (int k=0; k<inverseregion.noderange.value[k]; k+=2) {
                    for (int l=inverseregion.nodelist.value[k+k*inverseregion.nodelist.height]; l<=inverseregion.nodelist.value[k+(k+1)*inverseregion.nodelist.height]; l++) {
                        int y = l;
                        if (y == 0) { //Top-Right Corner
                            paint_buffer[y+x*plot_extrap.height] = (plot_extrap.value[y+(x-1)*plot_extrap.height]+plot_extrap.value[y+x*plot_extrap.height]+plot_extrap.value[y+1+(x-1)*plot_extrap.height]+plot_extrap.value[y+1+x*plot_extrap.height])/4;
                        } else if (y == plot_extrap.height-1) { //Bottom-Right Corner
                            paint_buffer[y+x*plot_extrap.height] = (plot_extrap.value[y-1+(x-1)*plot_extrap.height]+plot_extrap.value[y-1+x*plot_extrap.height]+plot_extrap.value[y+(x-1)*plot_extrap.height]+plot_extrap.value[y+x*plot_extrap.height])/4;
                        } else { // Right side
                            paint_buffer[y+x*plot_extrap.height] = (plot_extrap.value[y-1+(x-1)*plot_extrap.height]+plot_extrap.value[y-1+x*plot_extrap.height]+plot_extrap.value[y+(x-1)*plot_extrap.height]+plot_extrap.value[y+x*plot_extrap.height]+plot_extrap.value[y+1+(x-1)*plot_extrap.height]+plot_extrap.value[y+1+x*plot_extrap.height])/6;
                        }
                    }
                }    
            } else { // Middle Columns
                for (int k=0; k<inverseregion.noderange.value[k]; k+=2) {
                    for (int l=inverseregion.nodelist.value[k+k*inverseregion.nodelist.height]; l<=inverseregion.nodelist.value[k+(k+1)*inverseregion.nodelist.height]; l++) {
                        int y = l;
                        if (y == 0) { //Top-Side
                            paint_buffer[y+x*plot_extrap.height] = (plot_extrap.value[y+(x-1)*plot_extrap.height]+plot_extrap.value[y+x*plot_extrap.height]+plot_extrap.value[y+(x+1)*plot_extrap.height]+plot_extrap.value[y+1+(x-1)*plot_extrap.height]+plot_extrap.value[y+1+x*plot_extrap.height]+plot_extrap.value[y+1+(x+1)*plot_extrap.height])/6;
                        } else if (y == plot_extrap.height-1) { //Bottom-Side
                            paint_buffer[y+x*plot_extrap.height] = (plot_extrap.value[y-1+(x-1)*plot_extrap.height]+plot_extrap.value[y-1+x*plot_extrap.height]+plot_extrap.value[y-1+(x+1)*plot_extrap.height]+plot_extrap.value[y+(x-1)*plot_extrap.height]+plot_extrap.value[y+x*plot_extrap.height]+plot_extrap.value[y+(x+1)*plot_extrap.height])/6;
                        } else { // Middle Area
                            paint_buffer[y+x*plot_extrap.height] = (plot_extrap.value[y-1+(x-1)*plot_extrap.height]+plot_extrap.value[y-1+x*plot_extrap.height]+plot_extrap.value[y+(x-1)*plot_extrap.height]+plot_extrap.value[y-1+(x+1)*plot_extrap.height]+plot_extrap.value[y+x*plot_extrap.height]+plot_extrap.value[y+(x+1)*plot_extrap.height]+plot_extrap.value[y+1+(x-1)*plot_extrap.height]+plot_extrap.value[y+1+x*plot_extrap.height]+plot_extrap.value[y+1+(x+1)*plot_extrap.height])/9;
                        }
                    }
                }    
            }
        }

        // Repaint plot_extrap
        for (int j=0; j<inverseregion.noderange.height; j++) {
            int x = j;
            for (int k=0; k<inverseregion.noderange.value[j]; k+=2) {
                for (int l=inverseregion.nodelist.value[j+k*inverseregion.nodelist.height]; l<=inverseregion.nodelist.value[j+(k+1)*inverseregion.nodelist.height]; l++) {
                    int y = l;
                    plot_extrap.value[y+x*plot_extrap.height] = paint_buffer[y+x*plot_extrap.height];
                }
            }
        }
    }
}
    
// Debugging ------------------------------------------------------------//    

// NOT THREADSAFE
void imshow(double *img,int width,int height) {
    mxArray *rhs_imshow[2];
    
    // Copy array
    rhs_imshow[0] = mxCreateDoubleMatrix(height, width, mxREAL);   
    double *img_buf = mxGetPr(rhs_imshow[0]);
    for (int i=0; i<width; i++) {
        for (int j=0; j<height; j++) {
            img_buf[j+i*height] = img[j+i*height];
        }
    }    
    
    rhs_imshow[1] = mxCreateDoubleMatrix(0, 0, mxREAL);      
    
    // Call imshow - call figure first to make sure previous figure is not overwritten
    mexCallMATLAB(0,NULL,0,NULL,"figure");  
    mexCallMATLAB(0,NULL,2,rhs_imshow,"imshow");    
}

void imshow(int *img,int width,int height) {
    mxArray *rhs_imshow[2];
    
    // Copy array
    rhs_imshow[0] = mxCreateDoubleMatrix(height, width, mxREAL); 
    double *img_buf = mxGetPr(rhs_imshow[0]);
    for (int i=0; i<width; i++) {
        for (int j=0; j<height; j++) {
            img_buf[j+i*height] = (double)img[j+i*height];
        }
    }    
    
    rhs_imshow[1] = mxCreateDoubleMatrix(0, 0, mxREAL);      

    // Call imshow
    mexCallMATLAB(0,NULL,0,NULL,"figure");  
    mexCallMATLAB(0,NULL,2,rhs_imshow,"imshow"); 
}

void imshow(bool *img,int width,int height) {
    mxArray *rhs_imshow[2];
    
    // Copy array
    rhs_imshow[0] = mxCreateDoubleMatrix(height, width, mxREAL);  
    double *img_buf = mxGetPr(rhs_imshow[0]);
    for (int i=0; i<width; i++) {
        for (int j=0; j<height; j++) {
            img_buf[j+i*height] = (double)img[j+i*height];
        }
    }    
    
    rhs_imshow[1] = mxCreateDoubleMatrix(0, 0, mxREAL);      

    // Call imshow
    mexCallMATLAB(0,NULL,0,NULL,"figure");  
    mexCallMATLAB(0,NULL,2,rhs_imshow,"imshow"); 
}
