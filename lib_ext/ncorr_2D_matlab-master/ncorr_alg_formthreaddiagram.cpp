// This function forms a thread diagram and a preview of the thread diagram overlaid
// on the image. This function guarantees each region is 4-way contiguous

#include <mex.h>
#include <math.h>
#include <vector>
#include <list>
#include "standard_datatypes.h" 
#include "ncorr_datatypes.h"

// ----------------------------------------------------//
// Main Class -----------------------------------------//
// ----------------------------------------------------//

class class_formthreaddiagram {
public:
    // Constructor
    class_formthreaddiagram(mxArray *plhs [ ],const mxArray *prhs [ ]);
    
    // Methods
    void analysis();
    void analyzepoint(const int &x_new,const int &y_new,const int &num_queue);

private:
    // Properties
    // Inputs:
    class_double_array threaddiagram;           // standard datatype
    class_double_array preview_threaddiagram;   // standard datatype
    class_integer_array generators;             // standard datatype
    class_logical_array regionmask;             // standard datatype
    std::vector<ncorr_class_img> img;           // ncorr datatype   
    
    // Other:
    std::vector<std::list<std::vector<int> > > queue;
    std::vector<int> coords_buffer;
};

class_formthreaddiagram::class_formthreaddiagram(mxArray *plhs [ ],const mxArray *prhs[ ]) {
    // Form inputs ------------------------------------//
    // input 1: threaddiagram
    get_double_array(threaddiagram,prhs[0]); 
    // input 2: preview_threaddiagram
    get_double_array(preview_threaddiagram,prhs[1]); 
    // input 3: generators
    get_integer_array(generators,prhs[2]); 
    // input 4: regionmask
    get_logical_array(regionmask,prhs[3]);  
    // input 5: image (reduced)
    get_imgs(img,prhs[4]); 
    
    // Check inputs - these are rudimentary 
    if (img[0].gs.height == threaddiagram.height &&
        img[0].gs.width == threaddiagram.width &&
        img[0].gs.height == preview_threaddiagram.height &&
        img[0].gs.width == preview_threaddiagram.width &&
        img[0].gs.height == regionmask.height &&
        img[0].gs.width == regionmask.width &&
        generators.width == 2 || (generators.width == 0 && generators.height == 0)) {
        // Resize queue - one per generator
        queue.resize(generators.height);
        
        // Resize coords_buffer - this is used to push coordinates into the queue
        coords_buffer.resize(2);        
        
        // Form/set Outputs -------------------------------//
        // outputs: None - inputs are modified in-place
        
        // Get Outputs ------------------------------------//
        // outputs: None - inputs are modified in-place
    } else {
        mexErrMsgTxt("Inputs incorrect.\n");
    }    
}

// ----------------------------------------------------//
// Main Class Methods ---------------------------------//
// ----------------------------------------------------//

void class_formthreaddiagram::analysis() {
    // Initialize threaddiagram to -1; Note: the entire overlay is modified so it does not need to be reset
    for (int i=0; i<threaddiagram.width; i++) { 
        for (int j=0; j<threaddiagram.height; j++) { 
            threaddiagram.value[j+i*threaddiagram.height] = -1.0;
        }
    }
        
    // Initialize queues by inserting generators and painting them    
    for (int i=0; i<generators.height; i++) {
        // Paint
        threaddiagram.value[generators.value[i+generators.height] + generators.value[i]*threaddiagram.height] = i; 
        
        // Store coords
        coords_buffer[0] = generators.value[i];                      // X-coord
        coords_buffer[1] = generators.value[i+generators.height];    // Y-coord
        
        // Push back
        queue[i].push_back(coords_buffer);    
    }      
    
    // Initialize total_queue, this keeps track of the total number of 
    // elements in all the queues. When this is zero, exit.
    int total_queue = generators.height; 
    while (total_queue > 0) {
        // Cycle over each queue and add/paint four neighbors (up, right, down, left)
        for (int i=0; i<(int)queue.size(); i++) {
            // Make sure queue is not empty
            if (queue[i].size() > 0) {              
                // analyzepoint adds points using push_back, so this preserves the front
                analyzepoint(queue[i].front()[0]-1,queue[i].front()[1],i);
                analyzepoint(queue[i].front()[0],queue[i].front()[1]-1,i);
                analyzepoint(queue[i].front()[0]+1,queue[i].front()[1],i);
                analyzepoint(queue[i].front()[0],queue[i].front()[1]+1,i);
                               
                // Pop front element since it's finished
                queue[i].pop_front();
            }
        }           
        
        // Cycle over queue and recalculate total_queue
        total_queue = 0;
        for (int i=0; i<(int)queue.size(); i++) { 
            total_queue += (int)queue[i].size();
        }        
    }   
        
    // Go over entire mask to form preview_threaddiagram. Add a double highlight for points on the boundary between regions,
    // and then add just a normal highlight for the interior regions
    for (int i=0; i<regionmask.width; i++) {        
        int x = i;
        for (int j=0; j<regionmask.height; j++) {
            int y = j;
            if (regionmask.value[y+x*regionmask.height]) {
                // This is an interior point 
                if (x == 0 || x ==  regionmask.width-1 || y == 0 || y == regionmask.height-1) {
                    // This is a point on the boundary of regionmask, since its true, its automatically a boundary point
                    // Use a double highlight
                    preview_threaddiagram.value[y+x*preview_threaddiagram.height] = 2*img[0].max_gs;
                } else {
                    // This is an interior point of the region mask, see if a 9 neighborhood region contains two different values; if so it's a boundary point
                    double val_center = threaddiagram.value[y+x*threaddiagram.height];    
                    bool same = true;
                    for (int k=-1; k<=1; k++) {
                        for (int l=-1; l<=1; l++) {
                            double val_buffer = threaddiagram.value[(y+l)+(x+k)*threaddiagram.height];
                            if (val_center != val_buffer) {
                                same = false;
                                break;
                            }
                        }
                        if (!same) {
                            // This is a "between" point; double highlight
                            preview_threaddiagram.value[y+x*preview_threaddiagram.height] = 2*img[0].max_gs;
                        } else {
                            // This is an interior point; single highlight
                            preview_threaddiagram.value[y+x*preview_threaddiagram.height] = img[0].gs.value[y+x*img[0].gs.height]+img[0].max_gs;
                        }
                    }
                }    
            } else {
                // Point is outside region - Set this as just a normal greyscale value
                preview_threaddiagram.value[y+x*preview_threaddiagram.height] = img[0].gs.value[y+x*img[0].gs.height];
            }                  
        }    
    }    
}

void class_formthreaddiagram::analyzepoint(const int &x_new,const int &y_new,const int &num_queue) {
// This function checks to make sure point is valid and then paints/adds it
    if (y_new >= 0 && y_new < threaddiagram.height && 
        x_new >= 0 && x_new < threaddiagram.width) {
        // Check to make sure value is still -1 and within the regionmask
        if (threaddiagram.value[y_new + x_new*threaddiagram.height] == -1.0 &&
            regionmask.value[y_new + x_new*regionmask.height]) {    
            // Paint
            threaddiagram.value[y_new + x_new*threaddiagram.height] = (double)num_queue;   

            // Store coords
            coords_buffer[0] = x_new;    // X-coord
            coords_buffer[1] = y_new;    // Y-coord

            // Push back 
            queue[num_queue].push_back(coords_buffer);   
        }
    }
}

void mexFunction(int nlhs,mxArray *plhs[ ],int nrhs,const mxArray *prhs[ ]) {
    if (nrhs == 5 && nlhs == 0) {
        // Create formthreaddiagram
        class_formthreaddiagram formthreaddiagram(plhs,prhs);
        
        // Run analysis and assign outputs
        formthreaddiagram.analysis();
    } else {
        mexErrMsgTxt("Only five inputs and no output arguments.\n");
    }
}
