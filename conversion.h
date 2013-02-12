#ifndef CONVERSION
#define CONVERSION


/* Convert an array of elements, identified by their number, to a decimal
 * decimal number describing that state. For example, passing in the array 
 * {1, 3, 5} returns the number 2^0 + 2^2 + 2^4 = 22
 */

int convElementListToState(int *elements, int size);

void decToBinary(int decimalValue, int binaryValueSize, int *binaryValue);

int binaryToDec(int *binaryValue, int binaryValueSize);

#endif