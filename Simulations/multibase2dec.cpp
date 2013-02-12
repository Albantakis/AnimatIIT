#include "mex.h"
#include "math.h"
#include "matrix.h"

void MultiBaseToDec(int decimalValue, int* state_size, mwSize state_size_vec_size, double* multibase_digits);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    
    int decimalValue;
    double *state_size_vec_double, double *state_values;
    mwSize num_nodes;
    
    decimalValue =  mxGetScalar(prhs[0]);
    state_size_vec_double = mxGetPr(prhs[1]);
    state_size_vec_size = mxGetN(prhs[1]);

    int state_size_vec[state_size_vec_size];
    for(int i = 0; i < state_size_vec_size; i++)
        state_size_vec[i] = state_size_vec_double[i];
    
    plhs[0] = mxCreateDoubleMatrix(state_size_vec_size, 1, mxREAL);
    double* multibase_array = mxGetPr(plhs[0]);
    
    decToMultiBase(decimalValue, state_size_vec, state_size_vec_size, multibase_array);

}

void decToMultiBase(int decimalValue, int* state_size, mwSize state_size_vec_size, double* multibase_digits)
{
    
    for(mwSize i = 0; i < state_size_vec_size; i++) {
        
        *multibase_digits = decimalValue % *state_size;
        decimalValue = floor(decimalValue/ *state_size);
        multibase_digits++;
        state_size++;
    }
    
}