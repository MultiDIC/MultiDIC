// These are function definitions used to obtain standard inputs to mex files
// Note that all mx* functions are not thread safe. Also note that there are 
// some "iffy" things in place here, namely how ints and logicals are passed, 
// but the implementation used is simple and works on the majority of systems.

#include <mex.h>
#include "standard_datatypes.h"

// ----------------------------------------------------------------------//
// This is for string input ---------------------------------------------//
// ----------------------------------------------------------------------//

// NOT THREAD SAFE 
void get_string(std::string &string,const mxArray *mat_buf) { 
    // NOTE: This is potentially dangerous because user could change string length
    string = std::string(mxArrayToString(mat_buf));
}

// ----------------------------------------------------------------------//
// This is for a double scalar input ------------------------------------//
// ----------------------------------------------------------------------//

// NOT THREAD SAFE 
void get_double_scalar(double &scalar,const mxArray *mat_buf) {
    // Check input
    if (mxIsClass(mat_buf,"double")) {
        if (mxGetN(mat_buf) == 1 &&  mxGetM(mat_buf) == 1) {
            // At this point, input is correct
            scalar =  *mxGetPr(mat_buf);
        } else {
            mexErrMsgTxt("Double scalar is not of size == [1 1].\n");
        }
    } else {
        mexErrMsgTxt("Double scalar is not double.\n");
    }
}

// ----------------------------------------------------------------------//
// This is for an integer scalar input ----------------------------------//
// ----------------------------------------------------------------------//

// NOT THREAD SAFE 
void get_integer_scalar(int &scalar,const mxArray *mat_buf) {
    // Check input - NOTE: that most systems use 32 bit ints.
    // To my knowledge both 32 bit and 64 bit unix and windows
    // use ILP32 (32 bit systems), LP64 (64 bit unix systems), 
    // and LLP64 (64 bit windows) which specify 32 bit ints. 
    // This is from stackoverflow. However, be aware that 
    // int doesnt necessarily have to be 32 bit.
    if (mxIsClass(mat_buf,"int32")) {
        if (mxGetN(mat_buf) == 1 &&  mxGetM(mat_buf) == 1) {
            // At this point, input is correct
            scalar = *((int *)mxGetData(mat_buf));
        } else {
            mexErrMsgTxt("Integer scalar is not of size == [1 1].\n");
        }
    } else {
        mexErrMsgTxt("Integer scalar is not int32.\n");
    }
}

// ----------------------------------------------------------------------//
// This is for a logical scalar input -----------------------------------//
// ----------------------------------------------------------------------//

// NOT THREAD SAFE 
void get_logical_scalar(bool &scalar,const mxArray *mat_buf) {
    // Check input - NOTE: that proper implementation would
    // technically use mxLogical, but in both 32/64 bit linux
    // and windows systems I've tested, bool has worked. 
    // Maybe update this in the future.
    if (mxIsClass(mat_buf,"logical")) {
        if (mxGetN(mat_buf) == 1 &&  mxGetM(mat_buf) == 1) {
            // At this point, input is correct
            scalar = *mxGetLogicals(mat_buf);
        } else {
            mexErrMsgTxt("Logical scalar is not of size == [1 1].\n");
        }
    } else {
        mexErrMsgTxt("Logical scalar is not boolean.\n");
    }
}

// ----------------------------------------------------------------------//
// This is for a double array input -------------------------------------//
// ----------------------------------------------------------------------//

// THREAD SAFE 
class_double_array::class_double_array() {
    width = 0;
    height = 0;
    value = NULL;
}
   
// THREAD SAFE  
void class_double_array::reset() {
    for (int i=0; i<width; i++) { 
        for (int j=0; j<height; j++) {
            value[j+i*height] = 0.0;
        }
    }
}

// NOT THREAD SAFE 
void class_double_array::alloc(const int &h,const int &w) {
    if (value == NULL) {
        width = w;
        height = h;
        value = (double *)mxCalloc(h*w,sizeof(double));
    } else {
        mexErrMsgTxt("Memory has already been allocated when attempting to alloc.\n");
    }
}

// NOT THREAD SAFE 
void class_double_array::free() {
    if (value != NULL) {
        width = 0;
        height = 0;
        mxFree(value);
        value = NULL;
    } else {
        mexErrMsgTxt("Memory has not been allocated yet when attempting to free.\n");    
    }
}

// NOT THREAD SAFE 
void get_double_array(class_double_array &array,const mxArray *mat_buf) {
    // Check input
    if (mxIsClass(mat_buf,"double")) {
        // At this point input is correct
        array.width = (int)mxGetN(mat_buf);
        array.height = (int)mxGetM(mat_buf);
        array.value = mxGetPr(mat_buf);
    } else {
        mexErrMsgTxt("Double array is not double.\n");
    }
}

// ----------------------------------------------------------------------//
// This is for an integer array input -----------------------------------//
// ----------------------------------------------------------------------//

// THREAD SAFE 
class_integer_array::class_integer_array() {
    width = 0;
    height = 0;
    value = NULL;
}

// THREAD SAFE 
void class_integer_array::reset() {
    for (int i=0; i<width; i++) { 
        for (int j=0; j<height; j++) {
            value[j+i*height] = 0;
        }
    }
}

// NOT THREAD SAFE
void class_integer_array::alloc(const int &h,const int &w) {
    if (value == NULL) {
        width = w;
        height = h;
        value = (int *)mxCalloc(h*w,sizeof(int));
    } else {
        mexErrMsgTxt("Memory has already been allocated when attempting to alloc.\n");
    }
}

// NOT THREAD SAFE
void class_integer_array::free() {
    if (value != NULL) {
        width = 0;
        height = 0;
        mxFree(value);
        value = NULL;
    } else {
        mexErrMsgTxt("Memory has not been allocated yet when attempting to free.\n");    
    }
}

// NOT THREAD SAFE
void get_integer_array(class_integer_array &array,const mxArray *mat_buf) {
    // Check input - check notes for the int scalar about using the
    // native int type.
    if (mxIsClass(mat_buf,"int32")) {
        // At this point input is correct
        array.width = (int)mxGetN(mat_buf);
        array.height = (int)mxGetM(mat_buf);
        array.value = (int *)mxGetData(mat_buf);
    } else {
        mexErrMsgTxt("Integer array is not int32.\n");
    }
}

// ----------------------------------------------------------------------//
// This is for a logical array input ------------------------------------//
// ----------------------------------------------------------------------//

// THREAD SAFE
class_logical_array::class_logical_array() {
    width = 0;
    height = 0;
    value = NULL;
}

// THREAD SAFE
void class_logical_array::reset() {
    for (int i=0; i<width; i++) { 
        for (int j=0; j<height; j++) {
            value[j+i*height] = false;
        }
    }
}

// NOT THREAD SAFE 
void class_logical_array::alloc(const int &h,const int &w) {
    if (value == NULL) {
        width = w;
        height = h;
        value = (bool *)mxCalloc(h*w,sizeof(bool));
    } else {
        mexErrMsgTxt("Memory has already been allocated when attempting to alloc.\n");
    }
}

// NOT THREAD SAFE 
void class_logical_array::free() {
    if (value != NULL) {
        width = 0;
        height = 0;
        mxFree(value);
        value = NULL;
    } else {
        mexErrMsgTxt("Memory has not been allocated yet when attempting to free.\n");    
    }
}

// NOT THREAD SAFE 
void get_logical_array(class_logical_array &array,const mxArray *mat_buf) {
    // Check input - check notes for the logical scalar about using the
    // native bool type.
    if (mxIsClass(mat_buf,"logical")) {
        // At this point input is correct
        array.width = (int)mxGetN(mat_buf);
        array.height = (int)mxGetM(mat_buf);
        array.value = mxGetLogicals(mat_buf);
    } else {
        mexErrMsgTxt("Logical array is not bool.\n");
    }
}
