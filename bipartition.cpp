#include "engine.h"

int bipartition(int *set, int setSize, int nMax, int option, mxArray *B1, mxArray *B2) {
    
    
    if(option == 1)
        nMax = floor(nMax/2);
    
    nB = 0;
    for(int i = 0; i <= nMax; i++)
        nB = nB + nchoosek(setSize,i);
    
    int cellDims[2] = (nB,1);
    B1 = mxCreateCellArray (2,cellDims);
    B2 = mxCreateCellArray (2,cellDims);
    
    int iB = 0;
    
    //open up MATLAB engine to use the nchoosek function
    // TODO: !!! CAN WE USE mexCallMATLAB INSTEAD!?
    Engine *ep;
    
    if (!(ep = engOpen("\0"))) {
		fprintf(stderr, "\nCan't start MATLAB engine\n");
		//return EXIT_FAILURE;
	}
    
    for (int i = 0; i <= nMax; i++) {
        
        B1cellPtr = mxGetCell(B1, iB);
        B2cellPtr = mxGetCell(B2, iB);
        
        if(i == 0) {
            
            B1cellPtr = NULL;
            B2cellPtr = set; //do we need memcpy here?
            iB++;
        } else {
            
            mxArray *iMATLAB = mxCreateDoubleMatrix(1, 1, mxREAL);
            memcpy((void *)mxGetPr(iMATLAB), (void *)i, sizeof(i));
            mxArray *setSizeMATLAB = mxCreateDoubleMatrix(1, 1, mxREAL);
            memcpy((void *)mxGetPr(setSizeMATLAB), (void *)setSize, sizeof(setSize));
            
            engPutVariable(ep, "i", iMATLAB);
            engPutVariable(ep, "setSize", setSizeMATLAB);
            
            engEvalString(ep, "C_b = nchoosek(1:setSize,i);");
            
            int *chooseSets;
            mxArray *C_bmxArray = engGetVariable(ep,"C_b");
            double *temp = mxGetPr(C_bmxArray);
            int chooseSetsSizeCols = mxGetN(C_bmxArray);
            int chooseSetsSizeRows = mxGetM(C_bmxArray);
            int chooseSetSize = chooseSetsSizeCols*chooseSetsSizeRows;
            int chooseSets[chooseSetSize];
            for(int j = 0; j < chooseSetSize; j++)
                chooseSets[j] = temp[j];
            
            for(int j = 0; j < C_bSizeRows; j++) {

                int chooseSetIndex = 0;
                int B1index = 0;
                int B2index = 0;
                for(int k =0; k < setSize; k++) {
                    if((j + chooseSetIndex*C_bSizeRows) < chooseSetSize) 
                            && chooseSet[j + chooseSetIndex*C_bSizeRows] == set[k]) {
                        
                        B1cellPtr[B1index] = set[k];
                        chooseSetIndex++;
                        B1index++;
                    } else {
                        B2cellPtr[B2index] = set[k];
                        B2index++;
                    }      
                    iB++;
                }
            }
            
            

            mxDestroyArray(iMATLAB);
            mxDestroyArray(setSize);
            
        }
    }
    
    engEvalString(ep, "close;");
    
    return nB;
}
    
int nchoosek(n,k) {
    
    return factorial(n)/(factorial(k) * factorial(n-k));
}

int factorial(int n) {
    
    if (n == 1 || n == 0)
        return 1;
    return n * factorial(n - 1);
}