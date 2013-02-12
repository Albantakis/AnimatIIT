#include "mex.h"
#include "math.h"
#include "matrix.h"

void decToBinary(double decimalValue, mwSize binaryValueSize, double *binaryValue);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    
    double decimalValue;      
    mwSize binaryValueSize;
    
    decimalValue = mxGetScalar(prhs[0]);
    binaryValueSize = mxGetScalar(prhs[1]);
    
    plhs[0] = mxCreateDoubleMatrix(binaryValueSize, 1, mxREAL);
    double *binaryArray = mxGetPr(plhs[0]);
    
    //mxArray *binaryValueArray = mxCreateDoubleMatrix(binaryValueSize, 1, mxREAL);

    decToBinary(decimalValue, binaryValueSize, binaryArray);

}

void decToBinary(double decimalValue, mwSize binaryValueSize, double *binaryValue)
{
    
    for(mwSize i = 0; i < binaryValueSize; i++) {
        
        *binaryValue = decimalValue - floor(decimalValue/2)*2;
        decimalValue = floor(decimalValue/2);
        binaryValue++;
    }
    
}