// Minh Vu - 979167331 - minhvu@pdx.edu

#include <iostream>
#include <fstream>
#include <string>
using namespace std;

int ConvertBinaryCharArray2Decimal(char* BinCharArray) 		// convert an array of character (8 characters) to its decimal base value
{
	int a, number;
	int len = 7;
	int sum = 0;
	for (int i = 0; i <= len; i++)
	{
		a = 1;
		number = (BinCharArray[i] - '0'); 			// convert char to num
		a = a << (len - i);
		sum += number * a; 
	}
	return sum;
}
void Write2File(char *buffer, char* filename, int &flagwrite)			// Write values to binary file
{
	int decimalValue = 0;
	ofstream outputfile;
	if (flagwrite == 0) 			// check the flag to see if this is the 1st time to write in this file
	{
		outputfile.open(filename); 					// just open the file to write in from beginning
	}
	else
	{
		outputfile.open(filename, ios::app); 		// continue to append to the file
	}
	decimalValue = ConvertBinaryCharArray2Decimal(buffer); // convert from binary to decimal base
	unsigned char output[] = { decimalValue }; 					// set it as unsigned char (binary)
	outputfile.write( (const char *) output, 1);   // 1 byte
	flagwrite = 1; 							// already open and write in the file
	outputfile.close();
	
}
int main(int argc, char * argv[])
{
	ifstream huffmantable(argv[1]);							//receive the file name as the 1st argument
	string var1, var2;
	int numOfLines = 0;
	string line;
	while(getline(huffmantable,line))						// count the number of lines in the 1st huffman table file
		++numOfLines;
	huffmantable.clear();									// clear and reset the ios to beginning of the file
	huffmantable.seekg(0, ios::beg);						//
	string *key;											// key and codeword are similar to the idea of hashmap (or dictionary)
	key = new string[numOfLines];							// I am not sure why I could not use vector in C++ anymore
	string *codeword;										
	codeword = new string[numOfLines];
	int i = 0;
	while(huffmantable >> var1 >> var2)						// read the huffman table
	{
		key[i] = var1;										// store the file in 2 places, 1 is the key, 1 is the codeword
		codeword[i] = var2;
		++i;
	}

	for(int i = 0; i < numOfLines; i++)						// display the huffman table on screen, for debugging purpose only
		cout << key[i] << " " << codeword[i] << endl;

	ifstream input(argv[2]);								// receive the 2nd argument as input file
	char newchar;
	string InputFileString;
	while((newchar = input.get()) != EOF)					// read the whole file characters into memory at InputFileString variable
	{
		 InputFileString += newchar;
	}

	ofstream myfile;
	string result;

	for(int i = 0; i < InputFileString.length(); i++)			// compare the input with the huffman table to map the whole result	
	{															// into a variable (type string)
		for(int j = 0; j < numOfLines; j++)		
		{

			if((int)InputFileString[i]  == stoi(key[j]))					// compare character and the keyword
			{
				result += codeword[j];										// write the codeword down if it is a match
			}
		}
	}

	cout << result;												//
	cout << endl;												//		Debugging purposes
	cout << result.length() << endl;							//
	
	int remain = result.length() % 8;							// check to see if the whole variable is enough to form whole bytes
	if(remain != 0)												// if not, padding zeros to form whole bytes
	{
		int zeros = 8 - remain;
		for(int i = 0; i < zeros; ++i)
		{
			result += '0';
		}
	}


	char * out;													// convert the 'result' variable from type: string to char*
	out = new char[result.length()];
	for(int i = 0; i < result.length(); ++i )
		out[i] = result[i];
	cout << out << endl;

	char* buff;													// Attempt to write result to binary file
	buff = new char[8];											// write 8 characters each time 
	int count = 0;					
	int flag = 0;								
	for(int i = 0, j = 0; i < result.length(); ++i)				// loop through the whole result (in char*)
	{	
		buff[j] = out[i];												
		++j;
		++count;
		if(count % 8 == 0)										// get enough 8 characters
		{
			Write2File(buff,argv[3],flag);					// write them to file
			j = 0;													// reset values
			count = 0;												//
		}
	}
	
	delete[] key;
	delete[] codeword;
	delete[] out;
	delete[] buff;
	return 0;
}