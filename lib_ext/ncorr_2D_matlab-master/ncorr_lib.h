// These are functions declarations used in ncorr. 
  
#ifndef NCORR_LIB_H
#define NCORR_LIB_H 

#include <vector>
#include <string> 
#include "standard_datatypes.h"
#include "ncorr_datatypes.h"

#define WAITBAR_UPDATE 1000  // Number of waitbar increments until an update is shown
#define LAMBDA 0.0000000001  // Cutoff for values approximately zero
#define FILTERITERATIONS 50  // Number of average filters applied when extrapolating data

enum OUT {CANCELLED = -1, FAILED = 0, SUCCESS = 1}; // These let ncorr know whether analysis was successful, failed, or was cancelled so it can be handled appropriately.

// ----------------------------------------------------------------------//
// WAITBAR STUFF --------------------------------------------------------//
// ----------------------------------------------------------------------//

class class_waitbar {
public:
    // Constructor
    class_waitbar();    // THREADSAFE
    // Destructor
    ~class_waitbar();   // NOT THREADSAFE    
    
    // Methods
    void start(const int &num_img,const int &total_imgs,const int &computepoints);     // NOT THREADSAFE
    void increment();                                                                  // THREADSAFE
    bool updateandcheck();                                                             // NOT THREADSAFE
    
private:
    // Properties
    int computepoints;
    int counter_total;
    int counter_sub;
    mxArray *lhs_waitbar_create[1]; 
    mxArray *lhs_getappdata[1];       
    mxArray *rhs_getappdata[2]; 
    mxArray *rhs_waitbar_update[2];
    mxArray *mat_waitbar_update_1;
    mxArray *mat_getappdata_1;
    double *fraction_waitbar_update;
    
    // Methods
    bool getappdata();    // NOT THREADSAFE                                                
};

// ----------------------------------------------------------------------//
// FUNCTIONS ------------------------------------------------------------//
// ----------------------------------------------------------------------//

double ncorr_round(const double &r);                // THREADSAFE                
int sign(const double &r);                          // THREADSAFE    
int mod_pos(const int &i,const int &n);             // THREADSAFE    
void form_boundary(std::vector<std::vector<int> > &vec_boundary,const std::vector<int> &point_topleft,const class_logical_array &mask,int &direc);     // THREADSAFE    

// Use this as a threadsafe version of a region, since it uses vectors.
struct vec_struct_region {
public:
    // Constructor
    vec_struct_region(); // THREADSAFE    

    // Properties
    std::vector<int> nodelist;
    std::vector<int> noderange;
    int height_nodelist;
    int width_nodelist;
    int upperbound;
    int lowerbound;
    int leftbound;
    int rightbound;
    int totalpoints;
};
    
void form_regions(std::vector<vec_struct_region> &region,bool &removed,const class_logical_array &mask,const int &cutoff,const bool &preservelength);               // THREADSAFE
void form_union(std::vector<vec_struct_region> &region_union,const std::vector<ncorr_class_region> &region,const class_logical_array &mask,const bool &inplace);    // THREADSAFE    

void cholesky(std::vector<double> &mat,bool &positivedef,const int &size_mat);    // THREADSAFE    
void forwardsub(std::vector<double> &vec,const std::vector<double> &mat,const int &size_mat);        // THREADSAFE    
void backwardsub(std::vector<double> &vec,const std::vector<double> &mat,const int &size_mat);       // THREADSAFE    
OUT interp_qbs(double &interp,const double &x_tilda,const double &y_tilda,const class_double_array &plot_interp,const class_logical_array &mask,const int &offset_x,const int &offset_y,const int &border_bcoef); // NOT THREADSAFE    
void expand_filt(class_double_array &plot_extrap,const ncorr_class_inverseregion &inverseregion);    // THREADSAFE    

// Debugging ------------------------------------------------------------//

void imshow(double *img,int width,int height);    // NOT THREADSAFE
void imshow(int *img,int width,int height);       // NOT THREADSAFE
void imshow(bool *img,int width,int height);      // NOT THREADSAFE

#endif /* NCORR_LIB_H */
