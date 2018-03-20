// This function creates a mask (logical array) basked on drawobject inputs, which can be a rectangle, ellipse, or polygon

#include <mex.h>
#include <math.h>
#include <vector>
#include <algorithm>
#include "standard_datatypes.h"
#include "ncorr_datatypes.h"
#include "ncorr_lib.h"

// ----------------------------------------------------//
// Drawobjects input ----------------------------------//
// ----------------------------------------------------//

struct local_struct_drawobjects { 
    // Constructor
    local_struct_drawobjects() {
    }

    // Properties
    class_double_array pos_imroi; 
    std::string type;
    std::string addorsub;
};

void get_drawobjects(std::vector<local_struct_drawobjects> &drawobjects,const mxArray *prhs) {
    // Check input
    if ((mxGetM(prhs) == 1 || (mxGetM(prhs) == 0 && mxGetN(prhs) == 0)) && mxIsClass(prhs,"struct")) {
        for (int i=0; i<(int)mxGetN(prhs); i++) {
            mxArray *mat_pos_imroi = mxGetField(prhs,i,"pos_imroi");
            mxArray *mat_type = mxGetField(prhs,i,"type");
            mxArray *mat_addorsub = mxGetField(prhs,i,"addorsub");
            if (mat_pos_imroi != 0 &&
                mat_type != 0 &&
                mat_addorsub != 0) {
                // Form drawobject
                local_struct_drawobjects drawobject_template;                    
                get_double_array(drawobject_template.pos_imroi,mat_pos_imroi);
                get_string(drawobject_template.type,mat_type);
                get_string(drawobject_template.addorsub,mat_addorsub);
                
                // Test inputs
                if (drawobject_template.type.compare(0,4,"poly") == 0) {
                    // Input is 'poly' - width must be 2 and height must not be zero. If height is zero it causes a crash
                    if (drawobject_template.pos_imroi.width != 2 || drawobject_template.pos_imroi.height == 0) {
                        mexErrMsgTxt("'poly' width must be 2 with a height of at least 1.\n");
                    } 
                } else if (drawobject_template.type.compare(0,4,"rect") == 0) {
                    // Input is 'rect' - width must be 4 and height must be 1
                    if ((drawobject_template.pos_imroi.width != 4) || (drawobject_template.pos_imroi.height != 1)) {
                        mexErrMsgTxt("'rect' width must be 4 and height must be 1.\n");
                    }   
                } else if (drawobject_template.type.compare(0,7,"ellipse") == 0) {
                    // Input is 'ellipse' - width must be 4 and height must be 1
                    if ((drawobject_template.pos_imroi.width != 4) || (drawobject_template.pos_imroi.height != 1)) {
                        mexErrMsgTxt("'ellipse' width must be 4 and height must be 1.\n");
                    }  
                } else {
                    mexErrMsgTxt("'type' field must be either 'poly', 'rect' or 'ellipse'.\n");
                }                                                      
                if (drawobject_template.addorsub.compare(0,3,"add") != 0 && drawobject_template.addorsub.compare(0,3,"sub") != 0) {
                    mexErrMsgTxt("'addorsub' field must be either 'add' or 'sub'.\n");
                }    
                
                // Store drawobject_template
                drawobjects.push_back(drawobject_template);    
            } else {
                mexErrMsgTxt("Some fields are missing for drawobjects.\n");
            }        
        }
    } else {
        mexErrMsgTxt("drawobjects must be a row vector or empty of class 'struct'.\n");
    } 
}

// ----------------------------------------------------//
// Other Structures -----------------------------------//
// ----------------------------------------------------//

struct local_struct_node {
    // Constructor
    local_struct_node(const int length_buffer) {
        length = 0;
        values.resize(length_buffer);
    }
    
    // Properties
    std::vector<int> values; 
    int length;
};

// ----------------------------------------------------//
// Main Class -----------------------------------------//
// ----------------------------------------------------//

class class_formmask {
public:
    // Constructor
    class_formmask(mxArray *plhs[ ],const mxArray *prhs [ ]);
    
    // Methods
    void analysis();

private:
    // Properties
    // Inputs:
    std::vector<local_struct_drawobjects> drawobjects;  // local datatype
    class_logical_array mask;                           // standard datatype
    
    // Outputs: None
};

class_formmask::class_formmask(mxArray *plhs[ ],const mxArray *prhs[ ]) {    
    // Get inputs ------------------------------------//
    // input 1: drawobjects
    get_drawobjects(drawobjects,prhs[0]);
    // input 2: mask    
    get_logical_array(mask,prhs[1]);   
    
    // Form/set Outputs -------------------------------//
    // outputs: None - inputs are modified in-place
    
    // Get Outputs ------------------------------------//
    // outputs: None - inputs are modified in-place
}

// ----------------------------------------------------//
// Main Class Methods ---------------------------------//
// ----------------------------------------------------//

void class_formmask::analysis() {    
    // Reset values
    mask.reset();  
    
    // Find maximum possible number of nodes per column for polygonal regions. 
    // Check each drawobject and take the max over all of them. This buffer
    // holds the nodes to paint for a column. I'm finding the max here so
    // I don't have to constantly resize it.
    int length_polynode_max = 0; 
    for (int i=0; i<(int)drawobjects.size(); i++) {
        if (drawobjects[i].type.compare(0,4,"poly") == 0) {
            // The number of points in the polygon is the absolute upperbound for the number of nodes
            if (drawobjects[i].pos_imroi.height > length_polynode_max) {
                length_polynode_max = drawobjects[i].pos_imroi.height;
            }
        }
    }
    // Form node buffer - this buffer will hold the nodes to paint in a column for polygonal regions
    local_struct_node polynode_buffer(length_polynode_max);

    // "add" drawobjects will paint the area with 1's while "sub" drawobjects will paint the area with 0's
    // Cycle over drawobjects and paint mask depending on fields
    for (int i=0; i<(int)drawobjects.size(); i++) {
        if (drawobjects[i].type.compare(0,4,"poly") == 0) {
            // Get bounds from max and min of x coordinates in pos_imroi. Make sure they do not go outside the mask
            // Note if polygon is empty then this part will freeze. Possibly fix later.
            int leftbound = (int)std::max(ceil(*std::min_element(drawobjects[i].pos_imroi.value,drawobjects[i].pos_imroi.value+drawobjects[i].pos_imroi.height)),0.0);
            int rightbound = (int)std::min(floor(*std::max_element(drawobjects[i].pos_imroi.value,drawobjects[i].pos_imroi.value+drawobjects[i].pos_imroi.height)),(double)mask.width-1.0);
            
            if (drawobjects[i].addorsub.compare(0,3,"add") == 0) {
                for (int j=leftbound; j<=rightbound; j++) {
                    int x = j;
                    polynode_buffer.length = 0;

                    // Find nodes for polygon - code from alienryderflex.com
                    int idx_node = drawobjects[i].pos_imroi.height-1; // Subtract one since the points given are closed
                    for (int k=0; k<drawobjects[i].pos_imroi.height; k++) {
                        if ((drawobjects[i].pos_imroi.value[k] < x && drawobjects[i].pos_imroi.value[idx_node] >= x) || drawobjects[i].pos_imroi.value[idx_node] < x && drawobjects[i].pos_imroi.value[k] >= x) {
                            polynode_buffer.values[polynode_buffer.length] = (int)ncorr_round(drawobjects[i].pos_imroi.value[k+drawobjects[i].pos_imroi.height]+(x-drawobjects[i].pos_imroi.value[k])/(drawobjects[i].pos_imroi.value[idx_node]-drawobjects[i].pos_imroi.value[k])*(drawobjects[i].pos_imroi.value[idx_node+drawobjects[i].pos_imroi.height]-drawobjects[i].pos_imroi.value[k+drawobjects[i].pos_imroi.height]));
                            ++polynode_buffer.length;
                        }
                        idx_node = k;
                    }

                    // Sort nodes
                    std::sort(polynode_buffer.values.begin(),polynode_buffer.values.begin()+polynode_buffer.length);

                    // Paint nodes - Must check to make sure y coords do not exceed mask
                    for (int k=0; k<polynode_buffer.length; k+=2) {
                        if (polynode_buffer.values[k] >= mask.height) { 
                            // top node is beyond the bottom of the mask
                            break;
                        }
                        if (polynode_buffer.values[k+1] > 0) { 
                            // bottom node comes after the top of the mask
                            if (polynode_buffer.values[k] < 0) { 
                                // top node comes before the top of the mask
                                polynode_buffer.values[k] = 0;
                            }
                            if (polynode_buffer.values[k+1] >= mask.height) { 
                                // bottom node comes after the bottom of the mask
                                polynode_buffer.values[k+1] = mask.height-1;
                            }
                            // At this point the nodes are within the mask bounds, so paint them
                            for (int l=polynode_buffer.values[k]; l<=polynode_buffer.values[k+1]; l++) { 
                                int y = l;
                                mask.value[y + x*mask.height] = true;
                            }
                        }
                    }
                }
            } else if (drawobjects[i].addorsub.compare(0,3,"sub") == 0) {
                for (int j=leftbound; j<=rightbound; j++) {
                    int x = j;
                    polynode_buffer.length = 0;

                    //Find nodes for polygon - code from alienryderflex.com
                    int idx_node = drawobjects[i].pos_imroi.height-1; // Subtract one since the points given are closed
                    for (int k=0; k<drawobjects[i].pos_imroi.height; k++) {
                        if ((drawobjects[i].pos_imroi.value[k] < x && drawobjects[i].pos_imroi.value[idx_node] >= x) || drawobjects[i].pos_imroi.value[idx_node] < x && drawobjects[i].pos_imroi.value[k] >= x) {
                            polynode_buffer.values[polynode_buffer.length] = (int)ncorr_round(drawobjects[i].pos_imroi.value[k+drawobjects[i].pos_imroi.height]+(x-drawobjects[i].pos_imroi.value[k])/(drawobjects[i].pos_imroi.value[idx_node]-drawobjects[i].pos_imroi.value[k])*(drawobjects[i].pos_imroi.value[idx_node+drawobjects[i].pos_imroi.height]-drawobjects[i].pos_imroi.value[k+drawobjects[i].pos_imroi.height]));
                            ++polynode_buffer.length;
                        }
                        idx_node = k;
                    }

                    // Sort nodes
                    std::sort(polynode_buffer.values.begin(),polynode_buffer.values.begin()+polynode_buffer.length);

                    // Paint nodes - Must check to make sure y coords do not exceed mask
                    for (int k=0; k<polynode_buffer.length; k+=2) {
                        if (polynode_buffer.values[k] >= mask.height) { 
                            // top node is beyond the bottom of the mask
                            break;
                        }
                        if (polynode_buffer.values[k+1] > 0) { 
                            // bottom node comes after the top of the mask
                            if (polynode_buffer.values[k] < 0) { 
                                // top node comes before the top of the mask
                                polynode_buffer.values[k] = 0;
                            }
                            if (polynode_buffer.values[k+1] >= mask.height) { 
                                // bottom node comes after the bottom of the mask
                                polynode_buffer.values[k+1] = mask.height-1;
                            }
                            // At this point the nodes are within the mask bounds, so paint them
                            for (int l=polynode_buffer.values[k]; l<=polynode_buffer.values[k+1]; l++) { 
                                int y = l;
                                mask.value[y + x*mask.height] = false;
                            }
                        }
                    }
                }
            }
        } else if (drawobjects[i].type.compare(0,4,"rect") == 0) {
            // WIDTH IS NOT TRADITIONAL PIXEL WIDTH - it is "already subtracted by 1"
            int upperbound = (int)std::max(ceil(drawobjects[i].pos_imroi.value[1]),0.0);
            int lowerbound = (int)std::min(floor(drawobjects[i].pos_imroi.value[1]+drawobjects[i].pos_imroi.value[3]),(double)mask.height-1.0);
            int leftbound = (int)std::max(ceil(drawobjects[i].pos_imroi.value[0]),0.0);
            int rightbound = (int)std::min(floor(drawobjects[i].pos_imroi.value[0]+drawobjects[i].pos_imroi.value[2]),(double)mask.width-1.0);
            
            if (drawobjects[i].addorsub.compare(0,3,"add") == 0) {
                for (int j=leftbound; j<=rightbound; j++) {
                    int x = j;
                    // Paint nodes
                    for (int k=upperbound; k<=lowerbound; k++) {
                        int y = k;
                        mask.value[y + x*mask.height] = true;
                    }
                }
            } else if (drawobjects[i].addorsub.compare(0,3,"sub") == 0) {
                for (int j=leftbound; j<=rightbound; j++) {
                    int x = j;
                    // Paint nodes
                    for (int k=upperbound; k<=lowerbound; k++) {
                        int y = k;
                        mask.value[y + x*mask.height] = false;
                    }
                } 
            }
        } else if (drawobjects[i].type.compare(0,7,"ellipse") == 0) {
            // WIDTH IS NOT TRADITIONAL PIXEL WIDTH - it is already subtracted by 1
            int leftbound = (int)std::max(ceil(drawobjects[i].pos_imroi.value[0]),0.0);
            int rightbound = (int)std::min(floor(drawobjects[i].pos_imroi.value[0]+drawobjects[i].pos_imroi.value[2]),(double)mask.width-1.0);
            double a = drawobjects[i].pos_imroi.value[2]/2.0; // half-width
            double b = drawobjects[i].pos_imroi.value[3]/2.0; // half-height

            if (drawobjects[i].addorsub.compare(0,3,"add") == 0) {
                for (int j=leftbound; j<=rightbound; j++) {
                    int x = j;
                    int upperbound = (int)std::max(ceil(drawobjects[i].pos_imroi.value[1]+b-b*sqrt(1.0-pow((((double)j-drawobjects[i].pos_imroi.value[0]-a)/a),2))),0.0);
                    int lowerbound = (int)std::min(floor(drawobjects[i].pos_imroi.value[1]+b+b*sqrt(1.0-pow((((double)j-drawobjects[i].pos_imroi.value[0]-a)/a),2))),double(mask.height)-1.0);
                    
                    //Paint nodes
                    for (int k=upperbound; k<=lowerbound; k++) {
                        int y = k;
                        mask.value[y + x*mask.height] = true;
                    }
                }
            } else if (drawobjects[i].addorsub.compare(0,3,"sub") == 0) {
                for (int j=leftbound; j<=rightbound; j++) {
                    int x = j;
                    int upperbound = (int)std::max(ceil(drawobjects[i].pos_imroi.value[1]+b-b*sqrt(1.0-pow((((double)j-drawobjects[i].pos_imroi.value[0]-a)/a),2))),0.0);
                    int lowerbound = (int)std::min(floor(drawobjects[i].pos_imroi.value[1]+b+b*sqrt(1.0-pow((((double)j-drawobjects[i].pos_imroi.value[0]-a)/a),2))),double(mask.height)-1.0);

                    //Paint nodes
                    for (int k=upperbound; k<=lowerbound; k++) {
                        int y = k;
                        mask.value[y + x*mask.height] = false;
                    }
                }
            }
        }
    }
}

void mexFunction(int nlhs,mxArray *plhs[ ],int nrhs,const mxArray *prhs[ ]) {
    // Make sure number of input and output arguments are correct
    if (nrhs == 2 && nlhs == 0) {
        // Create formmask
        class_formmask formmask(plhs,prhs);
        
        // Run analysis and assign outputs
        formmask.analysis();
    } else {
        mexErrMsgTxt("Only two inputs and no output arguments allowed.\n");
    }
}
