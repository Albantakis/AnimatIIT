#include "mex.h"
#include "math.h"
#include "matrix.h"
#include "string.h"
#include "stdio.h"
#include "conversion.h"

/* 
 * function prob_prod = prob_prod_comp(prob1,prob2,M,x0_p1,op_fb)
 * prob1 - probability distribution as a column vector
 * prob2 - probability distribution as a column vector
 * M - a subset of the elements as row vector of ints
 * x0_p1 - the int which codes for the state of...
 */
void prob_prod_comp(double *prob1, int prob1Size, double *prob2, int prob2Size, 
                    int *elements, int elementsSize, int *elementsSubset, int elementsSubsetSize, double *probOut, int opFB);

void decToBinary(int decimalValue, int binaryValueSize, int *binaryValue);

int binaryToDec(int *binaryValue, int binaryValueSize);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    
    double *prob1, *prob2, *probOut, *temp;
    int prob1Size, prob2Size, elementsSize, elementsSubsetSize, opFB;
        
    prob1 = mxGetPr(prhs[0]);
    prob1Size = mxGetM(prhs[0]);
    
    prob2 = mxGetPr(prhs[1]);
    prob2Size = mxGetM(prhs[1]);
    
    temp = mxGetPr(prhs[2]);
    elementsSize = mxGetN(prhs[2]);
    int elements[elementsSize];
    for(int i = 0; i < elementsSize; i++)
        elements[i] = temp[i];
    
    temp = mxGetPr(prhs[3]);
    elementsSubsetSize = mxGetN(prhs[3]);
    int elementsSubset[elementsSubsetSize];
    for(int i = 0; i < elementsSubsetSize; i++)
        elementsSubset[i] = temp[i];
    
    opFB = mxGetScalar(prhs[4]);
    
    plhs[0] = mxCreateDoubleMatrix(pow(2,elementsSize), 1, mxREAL);
    probOut = mxGetPr(plhs[0]);
    
    prob_prod_comp(prob1, prob1Size, prob2, prob2Size, elements, elementsSize, 
            elementsSubset, elementsSubsetSize, probOut, opFB);

}

void prob_prod_comp(double* prob1, int prob1Size, double* prob2, int prob2Size, 
                    int* elements, int elementsSize, int* elementsSubset, int elementsSubsetSize, double* probOut, int opFB) {

    if (!prob2Size)
        memcpy(probOut, prob1, sizeof(double) * prob1Size);
    
    else if (!prob1Size)
        memcpy(probOut, prob2, sizeof(double) * prob2Size);
    
    else {
        
        //find size of complement set
        int elementsSubsetComplementSize = elementsSize - elementsSubsetSize;
    
        
        //(re)build elementsSubset and elementsSubsetComplement in terms of the size M
        // for example if M = {1,3,4}, then the subsets can only contain {1,2,3}
        
        int elementsSubsetComplement[elementsSubsetComplementSize];
        
        int elementsSubIndex = 0;
        int elementsSubCompIndex = 0;
        for (int i = 0; i < elementsSize; i++) {
            if (elementsSubIndex < elementsSubsetSize && elementsSubset[elementsSubIndex] == elements[i]) {
                elementsSubset[elementsSubIndex] = i + 1;
                elementsSubIndex++;
            }
            else if (elementsSubCompIndex < elementsSubsetComplementSize)  {
                elementsSubsetComplement[elementsSubCompIndex] = i + 1;
                elementsSubCompIndex++;
            }
        }
        
//         for (int i = 0; i < elementsSubsetSize; i++)
//             mexPrintf("SUBSET ELEMENT: %d\n", elementsSubset[i]);
//         
//         for (int i = 0; i < elementsSubsetComplementSize; i++)
//             mexPrintf("SUBSET COMP ELEMENT: %d\n", elementsSubsetComplement[i]);

        if (opFB == 3) {
            
//             double prob1_s[2^elementsSubsetSize][2^elementsSubsetSize];
//             prob1_s = reshape(prob1,[2^N1 2^N1]);
//             prob2_s = reshape(prob2,[2^N2 2^N2]);
//             prob_prod = zeros(2^N,2^N);
//             for i = 1:2^N
//                 xp_bs = trans2(i-1,N);
//                 xp_bs1 = xp_bs(x0_p1);
//                 xp_bs2 = xp_bs(x0_p2);
//                 xp_i1 = trans10(xp_bs1);
//                 xp_i2 = trans10(xp_bs2);
//                 for j=1: 2^N
//                     xf_bs = trans2(j-1,N);
//                     xf_bs1 = xf_bs(x0_p1);
//                     xf_bs2 = xf_bs(x0_p2);
//                     xf_i1 = trans10(xf_bs1);
//                     xf_i2 = trans10(xf_bs2);
//                     prob_prod(i,j) = prob1_s(xp_i1,xf_i1)*prob2_s(xp_i2,xf_i2);
//                 end
//             end
//             prob_prod = prob_prod(:);
        }// end if(opFB==3)
        else {
        
//     else
//         prob_prod = zeros(2^N,1);
//         x0_bs = zeros(N,1);
//         for i=1: 2^N1
//             x0_bs1 = trans2(i-1,N1);
//             x0_bs(x0_p1) = x0_bs1;
//             for j=1: 2^N2
//                 x0_bs2 = trans2(j-1,N2);
//                 x0_bs(x0_p2) = x0_bs2;
//                 x0_i = trans10(x0_bs);
//                 prob_prod(x0_i) = prob1(i)*prob2(j);
//             end
//         end
//         
//     end
//         
                    
            int x0_bs[elementsSize];
            
            for(int i = 0; i < pow(2,elementsSubsetSize); i++) {
                
                int x0_bs1[elementsSubsetSize];
                decToBinary(i, elementsSubsetSize, x0_bs1);
                
                for (int k = 0; k < elementsSubsetSize; k++) {
                    
                	x0_bs[elementsSubset[k] - 1] = x0_bs1[k];
                }
                
                for (int j = 0; j < pow(2,elementsSubsetComplementSize); j++) {
                    
                    int x0_bs2[elementsSubsetComplementSize];
                    decToBinary(j, elementsSubsetComplementSize, x0_bs2);
                    
                    
                    for (int k = 0; k < elementsSubsetComplementSize; k++) {
                    
                        x0_bs[elementsSubsetComplement[k] - 1] = x0_bs2[k];
                    }
                        
                    int x0_i = binaryToDec(x0_bs, elementsSize);
                    //mexPrintf("%d", x0_i);
                    probOut[x0_i] = prob1[i]*prob2[j];
                    
                } //end for  
            } //end for
        }//end else
    } //end else
}// end function

void decToBinary(int decimalValue, int binaryValueSize, int *binaryValue)
{
    
    for(mwSize i = 0; i < binaryValueSize; i++) {
        
        *binaryValue = decimalValue - floor(decimalValue/2)*2;
        decimalValue = floor(decimalValue/2);
        binaryValue++;
    }  
}

int binaryToDec(int *binaryValue, int binaryValueSize) {

    int decValue = 0;
    for(int i = 0; i < binaryValueSize; i++) {
        if (binaryValue[i]) {
            decValue = decValue + pow(2,i);
        }
    }
    return decValue;
}