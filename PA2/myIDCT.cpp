// Student: Minh Vu
// ID: 979167331
// File: myIDCT
//----------------------------
#include <cstdio>
#include <cstdlib>
#include <math.h>
#include <cstdarg>
#include <iostream>
#define MACROBLOCK 16
#define BLOCK 8
#define DCTCOSINE(a, b) (cos(((2.0 * b + 1.0) * a * M_PI) / 16.0))
using namespace std;
// Structure of pgm file
struct pgmFile							// hold data read from file
{
	char *type;
	int rows;
	int cols;
	float qscale;
	int qt[BLOCK][BLOCK];
	double block[BLOCK][BLOCK];			
};
// Some essential variables (global variables)
double block[BLOCK][BLOCK];				// 'main' block, use to store data of a particular block
pgmFile pgmData;
unsigned char **MemMacroBlock; 			// Macroblock in memory to hold values
int decoder[BLOCK*BLOCK];				// use for inverse zigzag order

//generate the decoder: generate 2 8x8 arrays. One is [1,2,3,4,..63]. The order is the 1st one but in zigzag order
// decoder holds the different values (in position) between these 2 arrays. 
void generateDecoder()
{
	int sample[BLOCK*BLOCK];
	int origin[BLOCK][BLOCK];						//[1,2,3,4,..,63]
	int count = 0;
	for(int m = 0; m < BLOCK;m++)					// generate the 1st array
		for(int n = 0; n < BLOCK;n++)
		{
			origin[m][n] = count;
			count++;
		}
	int k = 0;
	for(int i = 0; i < BLOCK*BLOCK; i++,k++)
		sample[i] = k;
	int res[BLOCK*BLOCK];						// 2nd array in zigzag order
	int j = 0;
	for( int i = 0; i < BLOCK*BLOCK; i++)		// run zigzag reorder to generate the 2nd array
	{
		if(i % 2 == 1)
		{
			for(int y = min(i,BLOCK - 1); y >= max(0,i - (BLOCK -1)); y--)
				res[j++] = origin[i-y][y];					
		}
		else {
			for(int x = min(i,BLOCK - 1); x >= max(0,i- (BLOCK -1));x--)
				res[j++] = origin[x][i-x];					
		}
	}
	for(int i = 0; i < BLOCK*BLOCK; i++)
	{
		int key = sample[i];
		for(int j = 0; j < BLOCK*BLOCK; j++)
		{
			if(key == res[j])
			{
				decoder[i] = j - i;			// determine the right position for each value
			}

		}
	}
}
// Inverse zigzag order for each block (after read from file to memory) 
void InverseZigZagOrder()
{
	int temp[BLOCK*BLOCK];
	int final[BLOCK*BLOCK];
	int n = 0;
	for(int i = 0; i < BLOCK; i++)
		for(int j = 0; j < BLOCK; j++)
		{
			temp[n] =pgmData.block[i][j];			// copy block to 1D array
			n++;
		}
	for(int i = 0; i < BLOCK*BLOCK; i++)			// inverse zigzag order of the 1D array 
	{
		int key = decoder[i];
		if(key == 0)
			final[i] = temp[i];
		else
			final[i] = temp[i+key];
	}
	n = 0;
	for(int i = 0; i < BLOCK; i++)
		for(int j = 0; j < BLOCK; j++)
			pgmData.block[i][j] = final[n++];		// copy 1D array back to block

}
// read the quantize file to memory
void handleQuantizeFile(FILE *image,FILE *qf)
{
	char* type = new char[sizeof(char)];
	rewind(image);
	fscanf(image, "%s", type);
	fscanf(image, "%d", &pgmData.cols);
	fscanf(image, "%d", &pgmData.rows);
	fscanf(image, "%f", &pgmData.qscale);
}
// multiply the block with quantize matrix and qvalue
void multiplyQuantize(double matrix[BLOCK][BLOCK], int qt[BLOCK][BLOCK], float qscale)
{
	for (int i = 0; i < BLOCK; i++)
	{
		for (int j = 0; j < BLOCK; j++)
		{
			matrix[i][j] = matrix[i][j]*qt[i][j] * qscale;
		}
	}
}
//copy blocks to macroblock 
void copyBlockValueToMemMacro(int count)
{

	int BlockInitX = ((count / 4) * 16) + (((count % 4 == 1) || (count % 4 == 3)) ? 8 : 0);
	int BlockInitY = ((count % 4 == 2) || (count % 4 == 3)) ? 8 : 0;
	int BlockEndX = BlockInitY + 8;
	int BlockEndY = BlockInitX + 8;
	for (int m = BlockInitY,i = 0; m < BlockEndX; m++, i++)
	{
		for (int n = BlockInitX,j = 0; n < BlockEndY; n++, j++)
		{
			MemMacroBlock[m][n] = (unsigned char)pgmData.block[i][j];
		}
	}
}
// write the whole macro block to file
void writeMacroblockToFile(FILE *output)
{
	int i, j;
	for (i = 0; i < MACROBLOCK; i++)
	{
		for (j = 0; j < pgmData.cols; j++)
		{
			fprintf(output, "%c", MemMacroBlock[i][j]);
		}
	}
}
// do IDCT by using the IDCT equation
void applyIDCTequation(double block[BLOCK][BLOCK])
{
	int u, v, x, y;
	double CuCv;
	double SumRow, SumCol;
	float tempBlock[BLOCK][BLOCK];			// temp block

	for (u = 0; u < BLOCK; u++)
	{
		for (v = 0; v < BLOCK; v++)
		{
			SumRow = 0;
			for (x = 0; x < BLOCK; x++)
			{
				SumCol = 0;
				for (y = 0; y < BLOCK; y++)
				{
					if ((x == 0) && (y == 0))			//Cu,Cv == 0
					{
						CuCv = 0.5;
					}
					else if ((x > 0) && (y > 0))		//Cu,Cv > 0
					{
						CuCv = 1.0;
					}
					else
					{
						CuCv = 1 / sqrt(2);
					}
					SumCol = SumCol + (CuCv * block[x][y] * DCTCOSINE(x, u) * DCTCOSINE(y, v)); // DCT transform formula
				}
				SumRow = SumRow + SumCol;
			}
			tempBlock[u][v] = SumRow / 4.0;				// store all calculation to temp Block
		}
	}
	for (u = 0; u < BLOCK; u++)
	{
		for (v = 0; v < BLOCK; v++)
		{
			block[u][v] = tempBlock[u][v];				// copy temp Block to 'main' Block
		}
	}
}

// read the quantize file and save to memory
void readQuantMatrix(FILE *qf)
{
	for (int i = 0; i < BLOCK; i++)
	{
		for (int j = 0; j < BLOCK; j++)
		{
			fscanf(qf, "%d", &pgmData.qt[i][j]);
		}
	}
}
// reset the offset of the pixel to [-127,128]
void resetOffset()
{
	int i, j;
	for (i = 0; i < BLOCK; i++)
	{
		for (j = 0; j < BLOCK; j++)
		{
			pgmData.block[i][j] -= 127; 
		}
	}
}
// round and crop coefficients
void roundcrop()
{
	int i, j;
	for (i = 0; i < BLOCK; i++)
	{
		for (j = 0; j < BLOCK; j++)
		{
			if (pgmData.block[i][j] < 0)
			{
				pgmData.block[i][j] = 0;
			}
			if (pgmData.block[i][j] > 255)
			{
				pgmData.block[i][j] = 255;
			}
			pgmData.block[i][j] = round(pgmData.block[i][j]);
		}
	}
}
// Decompress process in general
void Decompress()
{
	resetOffset();

	multiplyQuantize(pgmData.block, pgmData.qt, pgmData.qscale); 

	applyIDCTequation(pgmData.block);

	roundcrop();
}

// general IDCT function: read block from file then do IDCT on each block
void handleIDCTprocess(FILE *image, FILE *qf, FILE *output)
{
	int countBlock = 0;
	int BlockOffset;
	float readvalue;
	// jump to next line
	getc(image);
	
	MemMacroBlock = new unsigned char*[MACROBLOCK * sizeof(char *)];				// Allocate memory
	unsigned char* temp = new unsigned char[pgmData.cols * sizeof(char)];			// temp value

	for (int i = 0; i < MACROBLOCK; i++)
	{
		MemMacroBlock[i] = temp + (i * pgmData.cols);
	}

	int numBlock = (pgmData.rows / MACROBLOCK) * (pgmData.cols / MACROBLOCK) * 4;	// number of blocks
	int numTempBlock = 4 * (pgmData.cols / MACROBLOCK);								// number of temp blocks

	for (countBlock = 0; countBlock < numBlock; countBlock++)
	{
		fscanf(image, "%f", &readvalue);		// read the start position of the block
		fscanf(image, "%f", &readvalue);		// read the end position of the block
		for (int x = 0; x < BLOCK; x++)
		{
			for (int y = 0; y < BLOCK; y++)
			{
				fscanf(image, "%f", &readvalue);		// read the value
				pgmData.block[x][y] = readvalue;		// copy the value into memory
			}
		}
		InverseZigZagOrder();							// inverse Zigzag order on each block
		Decompress();									// do Decompression process on each block
		BlockOffset = countBlock % numTempBlock;
		copyBlockValueToMemMacro(BlockOffset);			// save each block to memory
		if (BlockOffset == (numTempBlock - 1))			// only write when got the whole macroblock
		{
			writeMacroblockToFile(output);				// write down the whole macroblock from memory to file
		}
	}
}
// Main function
int main(int argc, char *argv[])
{
	FILE *image = fopen(argv[1], "r"); 					// read the image
	FILE *qf = fopen(argv[2], "r"); 					// read quantize file
	FILE *output = fopen(argv[3], "wb");				// open output file
	generateDecoder();									
	handleQuantizeFile(image,qf);
	readQuantMatrix(qf);
	// write the first 3 lines to the output file
	fprintf(output, "%s\n", "P5");
	fprintf(output, "%d %d\n", pgmData.cols, pgmData.rows);
	fprintf(output, "%s\n", "255"); 
	// write the rest
	handleIDCTprocess(image, qf, output); 
	fclose(image);
	fclose(qf);
	fclose(output);
	return 0;
}