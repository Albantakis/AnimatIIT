#include "conversion.h"
#include "math.h"


int convElementListToState(int *elements, int size) {
    
    int stateValue = 0;
    
    for(int i=0; i < size; i++) {
        
        stateValue = stateValue + pow(2,elements[i]-1);
    }
}

void decToBinary(int decimalValue, int binaryValueSize, int *binaryValue)
{
    
    for(int i = 0; i < binaryValueSize; i++) {
        
        *binaryValue = decimalValue - floor(decimalValue/2)*2;
        decimalValue = floor(decimalValue/2);
        binaryValue++;
    }  
}

int binaryToDec(int *binaryValue, int binaryValueSize) {

    int decValue = 0;
    for(int i = 0; i < binaryValueSize; i++) {
        
        if (binaryValue[i])
            decValue = decValue + 2^i;
    }
    
    return decValue;
}