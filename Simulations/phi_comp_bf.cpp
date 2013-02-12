#include "mex.h"
#include "math.h"
#include "matrix.h"
#include "string.h"
#include "stdio.h"
//#include "engine.h"

using namespace std;


/*
 * Function declarations
 */
int convElementListToState(int *elements, int size);

int factorial(int n);

int nchoosek(int n, int k);

int normalization(mxArray *pastPartition1, mxArray *pastPartition2, mxArray *currentPartition1, mxArray *currentPartition2);

int bipartition(int *set, int setSize, int nMax, int option, mxArray *B1, mxArray *B2);

int* convertMxArrayToIntArray(const mxArray *inputMxArray);

void decToBinary(int decimalValue, int binaryValueSize, int *binaryValue);

int binaryToDec(int *binaryValue, int binaryValueSize);

int min(int int1, int int2);

double klDiv(double *probDist1, double *probDist2, int supportSize);

void prob_prod_comp(double *prob1, int prob1Size, double *prob2, int prob2Size, 
                    int *elements, int elementsSize, int *elementsSubset, int elementsSubsetSize, double *probOut, int opFB);

void phi_comp_bf(int *options, int *system, int systemSize, int *thisPartition, int thisPartitionSize, 
                    int *pastPartition, int pastPartitionSize,
                    int *futurePartition, int futurePartitionSize,
                    int *currentState, int currentStateSize,
                    double *transProbMatrix, const mxArray *pastRepertoires, const mxArray *futureRepertoires,
                    mxArray *MIP, mxArray *probProdMIP, mxArray *prob, double *phiMIP);

/*
 * Gateway function
 */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    
    //unpack options struct
    const mxArray *mxOptions = prhs[0];
    int options[3];
//     options[0] = mxGetScalar(mxGetField(mxOptions, 0, "op_context"));
//     options[1] = mxGetScalar(mxGetField(mxOptions, 0, "op_whole"));
//     options[2] = mxGetScalar(mxGetField(mxOptions, 0, "op_min"));
    
    options[0] = 0;
    options[1] = 0;
    options[2] = 1;
    
    int *system = convertMxArrayToIntArray(prhs[1]);
    int systemSize = mxGetN(prhs[1]);
    
    int *thisPartition = convertMxArrayToIntArray(prhs[2]);
    int thisPartitionSize = mxGetN(prhs[2]);
    
    int *pastPartition = convertMxArrayToIntArray(prhs[3]);
    int pastPartitionSize = mxGetN(prhs[3]);
    
    int *futurePartition = convertMxArrayToIntArray(prhs[4]);
    int futurePartitionSize = mxGetN(prhs[4]);
    
    int *currentState = convertMxArrayToIntArray(prhs[5]);
    int currentStateSize = mxGetN(prhs[5]);
    
    double *transProbMatrix = mxGetPr(prhs[6]);
    //we're going to need to get the size to implement the op_min = 0 option
    
    const mxArray *pastRepertoires = prhs[8];
    const mxArray *futureRepertoires = prhs[9];
    
    plhs[0] = mxCreateDoubleMatrix(2, 1, mxREAL);
    double *phiMIP = mxGetPr(plhs[0]);

    mxArray *prob = plhs[1];
    mxArray *probProdMIP = plhs[2];
    mxArray *MIP = plhs[3];
    
    phi_comp_bf(options, system, systemSize, thisPartition, thisPartitionSize, 
                    pastPartition, pastPartitionSize,
                    futurePartition, futurePartitionSize,
                    currentState, currentStateSize,
                    transProbMatrix, pastRepertoires, futureRepertoires,
                    MIP, probProdMIP, prob, phiMIP);
    
    
}

/* 
 * function [phi_MIP prob prob_prod_MIP MIP] = phi_comp_bf(options,M,x0,xp,xf,x0_s,p,b_table,BRs,FRs)
 */
void phi_comp_bf(int *options, int *system, int systemSize, int *thisPartition, int thisPartitionSize, 
                    int *pastPartition, int pastPartitionSize,
                    int *futurePartition, int futurePartitionSize,
                    int *currentState, int currentStateSize,
                    double *transProbMatrix, const mxArray *pastRepertoires, const mxArray *futureRepertoires,
                    mxArray *MIP, mxArray *probProdMIP, mxArray *prob, double *phiMIP) {
    
    
    int numPastPartitions, numFuturePartitions, numCurrentPartitions;
    mxArray *pastPartitions1, *pastPartitions2, *futurePartitions1, *futurePartitions2, 
            *currentPartitions1, *currentPartitions2, *probProducts;
    double *phiCandidates, *probPastFullSystem, *probFutureFullSystem;
    
    int op_context = options[0];
    int op_whole = options[1];
    int op_min = options[2];


    // unpartitioned transition repertoire
    if (op_context == 0) {
        
        int probDims[] = {2,1};
        prob = mxCreateCellArray(2, probDims);
        int nIndices = 2; int indicesPast[] = {convElementListToState(thisPartition, thisPartitionSize), 
                                             convElementListToState(pastPartition, pastPartitionSize)};
        int indexPast = mxCalcSingleSubscript(pastRepertoires, nIndices, indicesPast);
        probPastFullSystem = mxGetPr(mxGetCell(pastRepertoires,indexPast));
        
        int indicesFuture[] = {convElementListToState(thisPartition, thisPartitionSize), 
                           convElementListToState(futurePartition, futurePartitionSize)};
        int indexFuture = mxCalcSingleSubscript(futureRepertoires, nIndices, indicesFuture);
        probFutureFullSystem = mxGetPr(mxGetCell(futureRepertoires, indexFuture));
        
//         mxArray *probTemp1 = mxCreateDoubleMatrix(pow(2,systemSize),1,mxREAL);
//         mxSetData(probTemp1, (void *)probPastFullSystem);
//         int indices1[] = {0,0};
//         int index = mxCalcSingleSubscript(prob, 2, indices1);
//         mxSetCell(prob, index, probTemp1);
//         
//         mxArray *probTemp2 = mxCreateDoubleMatrix(pow(2,systemSize),1,mxREAL);
//         mxSetData(probTemp2, (void *)probFutureFullSystem);
//         int indices2[] = {1,0};
//         index = mxCalcSingleSubscript(prob, 2, indices2);
//         mxSetCell(prob, index, probTemp2);
//         
//         
    } else {}//TODO: put in progressive option
// 
// 
//     
//     if (pastPartitionSize !=0)
//         numPastPartitions = bipartition(pastPartition, pastPartitionSize, pastPartitionSize, 0, pastPartitions1, pastPartitions2);
//     else {
//         int dimSizes[] = {1,1};
//         pastPartitions1 = mxCreateCellArray(2,dimSizes);
//         int numSubs = 2; int subs[] = {0,0};
//         int index = mxCalcSingleSubscript(pastPartitions1, numSubs, subs);
//         //mxArray *cellPtr = mxGetCell(pastPartitions1, index);
//         //*cellPtr = NULL;
//         //cellPtr = mxGetCell(pastPartitions2, index);
//         //*cellPtr = NULL;
//         numPastPartitions = 1;
//     }
//     
//     if (futurePartitionSize !=0)
//         numFuturePartitions = bipartition(futurePartition, futurePartitionSize, futurePartitionSize, 0, futurePartitions1, futurePartitions2);
//     else {
//         int dimSizes[] = {1,1};
//         futurePartitions1 = mxCreateCellArray(2,dimSizes);
//         int numSubs = 2; int subs[] = {0,0};
//         int index = mxCalcSingleSubscript(futurePartitions1, numSubs, subs);
//         //mxArray *cellPtr = mxGetCell(futurePartitions1, index);
//         //*cellPtr = NULL;
//         //cellPtr = mxGetCell(futurePartitions2, index);
//         //*cellPtr = NULL;
//         numFuturePartitions = 1;
//     }
//     
//     numCurrentPartitions = bipartition(thisPartition, thisPartitionSize, thisPartitionSize, 0, currentPartitions1, currentPartitions2);
//     
//     if(op_min == 0) {}//TODO: this optoin
//     else {
//         
//         phiCandidates[numPastPartitions, numCurrentPartitions, numFuturePartitions, 2];
//         int numSubs  = 3;
//         int dimSizes[] = {numPastPartitions, numCurrentPartitions, 2};
//         probProducts = mxCreateCellArray(numSubs, dimSizes);
//         
//         for(int i = 0; i < numPastPartitions; i++) {
//             
//             mxArray *selectedPastPartition1 = mxGetCell(pastPartitions1, i);
//             mxArray *selectedPastPartition2 = mxGetCell(pastPartitions2, i);
//             
//             for (int j = 0; j < numCurrentPartitions; j++) {
//              
//                 mxArray *selectedCurrentPartition1 = mxGetCell(currentPartitions1, j);
//                 mxArray *selectedCurrentPartition2 = mxGetCell(currentPartitions2, j);
//                 
//                 int norm = normalization(selectedPastPartition1, selectedPastPartition2, selectedCurrentPartition1, selectedCurrentPartition2);
//                 
//                 double *probPartition1, *probPartition2;
//                 int probPartition1Size, probPartition2Size;
//                 double phi;
//                 
//                 for (int bf = 0; bf < 2; bf++) {
//                     
//                     if (norm != 0) {
//                         
//                         int *selectedCurrentPartition1Int, *selectedPastPartition1Int;
//                         
//                         if(bf == 0) {
//                             
//                             selectedCurrentPartition1Int = convertMxArrayToIntArray(selectedCurrentPartition1);
//                             selectedPastPartition1Int = convertMxArrayToIntArray(selectedPastPartition1);
//                             int numSubs = 2; int indices[numSubs];
//                             indices[0] = convElementListToState(selectedCurrentPartition1Int, mxGetN(selectedCurrentPartition1));
//                             indices[1] = convElementListToState(selectedPastPartition1Int, mxGetN(selectedPastPartition1));
//                             int index = mxCalcSingleSubscript(pastRepertoires, numSubs, indices);
//                             probPartition1 = mxGetPr(mxGetCell(pastRepertoires, index));
//                             probPartition1Size = mxGetM(mxGetCell(pastRepertoires, index));
//                             
//                             int *selectedCurrentPartition2Int = convertMxArrayToIntArray(selectedCurrentPartition2);
//                             int *selectedPastPartition2Int = convertMxArrayToIntArray(selectedPastPartition2);
//                             indices[0] = convElementListToState(selectedCurrentPartition2Int, mxGetN(selectedCurrentPartition2));
//                             indices[1] = convElementListToState(selectedPastPartition2Int, mxGetN(selectedPastPartition2));
//                             index = mxCalcSingleSubscript(pastRepertoires, numSubs, indices);
//                             probPartition2 = mxGetPr(mxGetCell(pastRepertoires, index));
//                             probPartition2Size = mxGetM(mxGetCell(pastRepertoires, index));
//                         }
//                         
//                         double *probProd;
//                         prob_prod_comp(probPartition1, probPartition1Size, probPartition2, probPartition2Size, 
//                             pastPartition, pastPartitionSize, selectedCurrentPartition1Int, mxGetN(selectedCurrentPartition1), probProd, 0);
//                         phi = klDiv(probPastFullSystem, probProd, pow(2,systemSize));
//                         int nSubs = 3; int subs[] = {i,j,bf};
//                         int index = mxCalcSingleSubscript(probProducts, nSubs, subs);
//                         mxArray *tempArray = mxCreateDoubleMatrix(pow(2,pastPartitionSize),1,mxREAL);
//                         mxSetData(tempArray,(void *)probProd);
//                         //mxArray *probProdsPtr = mxGetCell(probProducts, index);
//                         mxSetCell(probProducts, index, tempArray);
//                         
//                    }
//                    else {phi = DBL_MAX;}
//                    
//                    phiCandidates[i,j,bf,1] = phi;
//                    phiCandidates[i,j,bf,2] = phi/norm;
//                    
//                 }
//             }
//         }
//         
//         int nDims = 3; int dims1[] = {2,2,2};
//         MIP = mxCreateCellArray(nDims,dims1);
//         double phiMIPTemp[] = {0,0};
//         phiMIP = phiMIPTemp;
//         nDims = 2; int dims2[] = {2,1};
//         probProdMIP = mxCreateCellArray(nDims,dims2);
//         
//         for(int bf = 0; bf < 2; bf++) {
//             
//             double phiNormMin = DBL_MAX;
//             double phiMin = DBL_MAX;
//             int iMin, jMin;
//             
//             for(int i = 0; i < numPastPartitions; i++) {
//                 
//                 for(int j = 0; j < numCurrentPartitions; j++) {
//                     
//                     if (phiCandidates[i,j,bf,2] <= phiNormMin && phiCandidates[i,j,bf,1] <= phiMin) {
//                         phiMin = phiCandidates[i,j,bf,1];
//                         phiNormMin = phiCandidates[i,j,bf,2];
//                         iMin = i;
//                         jMin = j;
//                     }
//                 }   
//             
//             
//             phiMIP[bf] = phiMin;
//             int nDims = 3; int indices[] = {iMin, jMin, bf};
//             int index = mxCalcSingleSubscript(probProducts, nDims, indices);
//             mxArray *probProbEntry = mxGetCell(probProducts, index);
//             int indices2[] = {bf,0};
//             index = mxCalcSingleSubscript(probProdMIP,2,indices2);
//             mxSetCell(probProdMIP, index, probProbEntry);
//             
//             int indices3[] = {0,iMin};
//             index = mxCalcSingleSubscript(pastPartitions1, 2, indices3);
//             mxArray *temp = mxGetCell(pastPartitions1, index);
//             int indices4[] = {0,0,bf};
//             index = mxCalcSingleSubscript(MIP, 3, indices4);
//             mxSetCell(MIP, index, temp);
//             
//             int indices5[] = {0,iMin};
//             index = mxCalcSingleSubscript(pastPartitions2, 2, indices5);
//             temp = mxGetCell(pastPartitions1, index);
//             int indices6[] = {1,0,bf};
//             index = mxCalcSingleSubscript(MIP, 3, indices6);
//             mxSetCell(MIP, index, temp);
//             
//             int indices7[] = {0,jMin};
//             index = mxCalcSingleSubscript(currentPartitions1, 2, indices7);
//             temp = mxGetCell(pastPartitions1, index);
//             int indices8[] = {0,1,bf};
//             index = mxCalcSingleSubscript(MIP, 3, indices8);
//             mxSetCell(MIP, index, temp);
//             
//             int indices9[] = {0,jMin};
//             index = mxCalcSingleSubscript(currentPartitions2, 2, indices9);
//             temp = mxGetCell(pastPartitions1, index);
//             int indices10[] = {1,1,bf};
//             index = mxCalcSingleSubscript(MIP, 3, indices10);
//             mxSetCell(MIP, index, temp);
//             
//             }
//             
//         }
//         
//     }


}

int convElementListToState(int *elements, int size) {
    
//     int stateValue = 0;
//     
//     for(int i=0; i < size; i++) {
//         
//         stateValue = stateValue + pow(2,elements[i]-1);
//     }
    
    return 0;
}

int bipartition(int *set, int setSize, int nMax, int option, mxArray *B1, mxArray *B2) {
    
    
    if(option == 1)
        nMax = floor(nMax/2);
    
    int nB = 0;
    for(int i = 0; i <= nMax; i++)
        nB = nB + nchoosek(setSize,i);
    
    int numDims = 2; int cellDims[2] = {nB,1};
    B1 = mxCreateCellArray(numDims,cellDims);
    B2 = mxCreateCellArray(numDims,cellDims);
    
    int iB = 0;
    
    //open up MATLAB engine to use the nchoosek function
    // TODO: !!! CAN WE USE mexCallMATLAB INSTEAD!?
//     Engine *ep;
//     
//     if (!(ep = engOpen("\0"))) {
// 		fprintf(stderr, "\nCan't start MATLAB engine\n");
// 		//return EXIT_FAILURE;
// 	}
    
    for (int i = 0; i <= nMax; i++) {
        
        mxArray *B1mxArray;
        mxArray *B2mxArray;
        //double *B1, *B2;
        
        if(i == 0) {
            
            //*B1cellPtr = NULL;
            //memcpy(B2mxArray, set, setSize * sizeof(*set));
            mxSetData(B2mxArray, set);
            //B2cellPtr = set;... used memcpy above
            iB++;
        } else {
            
            //B1mxArray = mxCreateDoubleMatrix(1,i,mxREAL);
            //B2mxArray = mxCreateDoubleMatrix(1,setSize - i,mxREAL);
            
            int B1intArray[i];
            int B2intArray[setSize - i];
            
//             mxArray *iMATLAB = mxCreateDoubleMatrix(1, 1, mxREAL);
//             memcpy((void *)mxGetPr(iMATLAB), (void *)i, sizeof(i));
//             mxArray *setSizeMATLAB = mxCreateDoubleMatrix(1, 1, mxREAL);
//             memcpy((void *)mxGetPr(setSizeMATLAB), (void *)setSize, sizeof(setSize));
            
//             engPutVariable(ep, "i", iMATLAB);
//             engPutVariable(ep, "setSize", setSizeMATLAB);
            
            mxArray *iMATLAB;
            int iVec[] = {i};
            mxSetData(iMATLAB, iVec);
            
            int nlhs = 1;
            int nrhs = 2;
            
            mxArray *output_array[nlhs], *input_array[nrhs];
            int vec[setSize];
            for(int j = 0; j < setSize; j++)
                vec[j] = j+1;
            
            mxArray *vecMATLAB;
            mxSetData(vecMATLAB, vec);
            input_array[0] = vecMATLAB;
            input_array[1] = iMATLAB;
            
            mexCallMATLAB(nlhs, output_array, nrhs, input_array, "nchoosek");
            //engEvalString(ep, "C_b = nchoosek(1:setSize,i);");
            
            //int *chooseSets;
            mxArray *C_bmxArray = output_array[0];
//             double *temp = mxGetPr(C_bmxArray);
            int chooseSetsSizeCols = mxGetN(C_bmxArray);
            int chooseSetsSizeRows = mxGetM(C_bmxArray);
            int chooseSetSize = chooseSetsSizeCols*chooseSetsSizeRows;
//             int chooseSets[chooseSetSize];
//             for(int j = 0; j < chooseSetSize; j++)
//                 chooseSets[j] = temp[j];
            int *chooseSet = convertMxArrayToIntArray(C_bmxArray);
            
            for(int j = 0; j < chooseSetsSizeRows; j++) {

                int chooseSetIndex = 0;
                int B1index = 0;
                int B2index = 0;
                for(int k =0; k < setSize; k++) {
                    if(j + chooseSetIndex*chooseSetsSizeRows < chooseSetSize 
                            && chooseSet[j + chooseSetIndex*chooseSetsSizeRows] == set[k]) {
                        
                        B1intArray[B1index] = set[k];
                        chooseSetIndex++;
                        B1index++;
                    } else {
                        B2intArray[B2index] = set[k];
                        B2index++;
                    }      
                    iB++;
                }
             }
// 
//             mxDestroyArray(iMATLAB);
//             mxDestroyArray(setSizeMATLAB);
            
            mxSetData(B1mxArray,B1intArray);
            mxSetData(B2mxArray,B2intArray);
            
            mxSetCell(B1, iB, B1mxArray);
            mxSetCell(B2, iB, B2mxArray);
        }
        
    }
    
    //engEvalString(ep, "close;");
    
    return nB;
}

int nchoosek(int n,int k) {
    
    return factorial(n)/(factorial(k) * factorial(n-k));
}

int factorial(int n) {
    
    if (n == 1 || n == 0)
        return 1;
    return n * factorial(n - 1);
}

int normalization(mxArray *pastPartition1, mxArray *pastPartition2, mxArray *currentPartition1, mxArray *currentPartition2) {
    
    int pastPartition1Size = mxGetN(pastPartition1);
    int pastPartition2Size = mxGetN(pastPartition2);
    int currentPartition1Size = mxGetN(currentPartition1);
    int currentPartition2Size = mxGetN(currentPartition2);
    
    return min(currentPartition1Size, pastPartition2Size) + min(currentPartition2Size, pastPartition1Size);
}

int* convertMxArrayToIntArray(const mxArray *inputMxArray) {

    double *temp = mxGetPr(inputMxArray);
    int outputSizeCols = mxGetN(inputMxArray);
    int outputSetsSizeRows = mxGetM(inputMxArray);
    int outputSize = outputSizeCols*outputSetsSizeRows;
    int outputArray[outputSize];
    for(int i = 0; i < outputSize; i++)
        outputArray[i] = temp[i];
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

double klDiv(double *probDist1, double *probDist2, int supportSize) {
    
    double crossEntropy =0; double entropyDist1 = 0;
    
    for(int i = 0; i < supportSize; i++) {
        
        if (probDist2[i] != 0)
            crossEntropy -= probDist1[i] * log(probDist2[i]) / log(2);
        
        if (probDist2[i] != 0)
            entropyDist1 -= probDist1[i] * log(probDist1[i]) / log(2);
    }
    
    return crossEntropy - entropyDist1;
}

int min(int int1, int int2) {
    
    if (int1 > int2)
        return int2;
    else
        return int1;
}