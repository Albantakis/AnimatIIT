#include "mex.h"
#include "math.h"
#include "matrix.h"
#include "string.h"

void KLD(double* prob1, double *prob2, mwSize dist_size, double* kld_value);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    
    double* prob1, *prob2;
    mwSize dist_size;
    
    prob1 = mxGetPr(prhs[0]);
    prob2 = mxGetPr(prhs[1]);
    dist_size = mxGetM(prhs[0]);
    
    plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
    double* kld_value = mxGetPr(plhs[0]);
    
    
    //mxArray *binaryValueArray = mxCreateDoubleMatrix(binaryValueSize, 1, mxREAL);

    KLD(prob1, prob2, dist_size, kld_value);
}

void KLD(double* prob1, double *prob2, mwSize dist_size, double* kld_value)
{
    
    double h1 = 0;
    for (int i = 0; i < dist_size; i++) {
        
        bool changed = false;
        if (prob2[i] == 0) {
            prob2[i] = 1;
            changed = true;
        }
        
        h1 += log2(prob2[i]) * prob1[i];
        
        if (changed)
            prob2[i] = 0;
    }
    
    double h2 = 0;
    for (int i = 0; i < dist_size; i++) {
        
        bool changed = false;
        if (prob1[i] == 0) {
            prob1[i] = 1;
            changed = true;
        }
        
        h2 += log2(prob1[i]) * prob1[i];

        if (changed)
            prob1[i] = 0;
    }
    
    *kld_value = h2 - h1;
    
}


// matlab function

// prob2(prob2==0) = 1; % avoid log0 when computing entropy
// H1 = - sum(prob.*log2(prob2)) ;
// 
// prob(prob==0) = 1;
// H2 = - sum(prob.*log2(prob));
// H = H1 - H2;