// These are function declarations used to obtain ncorr-type inputs to mex files
// Some functions are not thread safe.

#ifndef NCORR_DATATYPES_H
#define NCORR_DATATYPES_H

#include <vector>
#include "standard_datatypes.h"

// ----------------------------------------------------------------------//
// NCORR_CLASS_IMG STUFF ------------------------------------------------//
// ----------------------------------------------------------------------//

// This is for an ncorr_class_img input ---------------------------------// 
class ncorr_class_img {
public:
    // Constructor    
    ncorr_class_img();    // THREADSAFE
    
    // Properties
    std::string type;
    class_double_array gs; 
    double max_gs;       
    class_double_array bcoef; 
    int border_bcoef;
};

// NOT THREADSAFE
void get_imgs(std::vector<ncorr_class_img> &images,const mxArray *prhs);

// ----------------------------------------------------------------------//
// NCORR_CLASS_ROI STUFF ------------------------------------------------//
// ----------------------------------------------------------------------//

class ncorr_class_region {
public:
    // Constructor
    ncorr_class_region();    // THREADSAFE

    // Properties
    class_integer_array nodelist;
    class_integer_array noderange;
    int upperbound;
    int lowerbound;
    int leftbound;
    int rightbound;
    int totalpoints;

    // Methods
    void alloc(const int &h,const int &w);    // NOT THREADSAFE
    void free();                              // NOT THREADSAFE    
};

// NOT THREADSAFE
void get_region(std::vector<ncorr_class_region> &region,const mxArray *prhs);

// Cirroi Stuff ---------------------------------------------------------//
struct struct_cirroi {
    // Constructor
    struct_cirroi();         // THREADSAFE
    
    // Properties
    ncorr_class_region region;
    class_logical_array mask; 
    int radius;
    int x;
    int y;
};

struct struct_info_cirroi {
    // Constructor
    struct_info_cirroi();    // THREADSAFE
    
    // Main properties
    ncorr_class_region region;
    class_logical_array mask; 
    int radius;
    int x;
    int y;
    
    // Additional Properties
    class_logical_array mask_buffer; // Used as additional storage so you dont have to overwrite original mask
    class_integer_array circletemplate;
    class_integer_array queue_buffer;
    std::vector<int> queue_nodelist;
    std::vector<int> queue_nodeindex; 
    std::vector<int> activelines; 
};

// This is for an ncorr_class_roi input ---------------------------------//
class ncorr_class_roi {
public:
    // Constructor
    ncorr_class_roi();    // THREADSAFE
        
    // Properties
    class_logical_array mask;
    std::vector<ncorr_class_region> region;
    std::vector<struct_cirroi> cirroi;
    
    // Methods
    void set_cirroi(const int &radius_i,const int &thread_total);              // NOT THREADSAFE
    void update_cirroi(const int &num_region,const int &thread_num);           // THREADSAFE
    void get_cirroi(const int &x_i,const int &y_i,const int &num_region,const bool &subsettrunc,const int &thread_num);    // THREADSAFE
    bool withinregion(const int &x_i,const int &y_i,const int &num_region);    // THREADSAFE

private:
    // Properties
    std::vector<struct_info_cirroi> info_cirroi;
};

// NOT THREADSAFE
void get_rois(std::vector<ncorr_class_roi> &rois,const mxArray *prhs);

// Inherited Types ------------------------------------------------------//

class ncorr_class_inverseregion : public ncorr_class_region {
public:
    // Constructor
    ncorr_class_inverseregion(ncorr_class_region &region,const int &border_extrap); // NOT THREADSAFE
};

#endif /* NCORR_DATATYPES_H */
