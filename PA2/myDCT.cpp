// Student: Minh Vu
// ID: 979167331
// File: myDCT
//----------------------------
#include <cstdio>
#include <cstdlib>
#include <math.h>
#include <cstdarg>
#include <iostream>
#define MACROBLOCK 16
#define BLOCK 8
#define PI 3.14159265359
#define DCTCOSINE(a, b) (cos(((2.0 * b + 1.0) * a * PI) / 16.0))
using namespace std;
// Structure of PGM 
struct pgmFile
{
	int rows;
	int cols;
	int format;
	unsigned char **value;
};
// Some essential variables (global variables)
int qMatrix[BLOCK][BLOCK]; 							// the matrix contains all quantize values
float qscale;										// the quantize scale
unsigned char macroBlock[MACROBLOCK][MACROBLOCK];   // Macro Block
double block[BLOCK][BLOCK];   						// block
pgmFile pgmData;									// the pgm data

// read the PGM file
void handlePGMFile(FILE *image)				
{
	char* type = new char[sizeof(char)]; 
	rewind(image); 							// start to read from the beginning
	fscanf(image, "%s", type); 				// read the type of the file
	fscanf(image, "%d", &pgmData.rows); 	// copy number of row
	fscanf(image, "%d", &pgmData.cols);		// copy number of colummn
	fscanf(image, "%d", &pgmData.format);	// copy the maxval
	pgmData.value = new unsigned char*[pgmData.cols * sizeof(char *)];			// allocate memory
	unsigned char *temp = new unsigned char[pgmData.rows * pgmData.cols * sizeof(char)];		// temp value

	for (int i = 0; i < pgmData.cols; i++)
	{
		pgmData.value[i] = temp + (i * pgmData.rows);
	}
	getc(image);
	// copy all values into struct pgmData (memory) 
	for (int x = 0; x < pgmData.cols; x++)
	{
		for (int y = 0; y < pgmData.rows; y++)
		{
			pgmData.value[x][y] = getc(image);
		}
	}
}
// copy all the values into qMatrix
void handleQuantFile(FILE *qf)					
{
	int x,y;
	for (x = 0; x < BLOCK; x++)
	{
		for (y = 0; y < BLOCK; y++)
		{
			fscanf(qf, "%d", &qMatrix[x][y]);
		}
	}
}
// convert quantized values to [0,255]
int ConvertTo255(int quantvalue)			
{	
	int i = 0;
	if (quantvalue < -127)
	{
		i = -127;
	}
	else if (quantvalue > 128)
	{
		i = 128;
	}
	else {
		i = quantvalue;
	}
	i = i + 127;
	return i;
}
// quantize the dctMatrix
void quantize(double dctMatrix[BLOCK][BLOCK], int qMatrix[BLOCK][BLOCK], float qscale)
{
	int quantvalue;
	for (int i = 0; i < BLOCK; i++)
	{
		for (int j = 0; j < BLOCK; j++)
		{
			quantvalue = (int)round(dctMatrix[i][j] / (qMatrix[i][j] * qscale)); 
			dctMatrix[i][j] = ConvertTo255(quantvalue);
		}
	}
}
// calculate the DCT transform equation
void calculateDCT(double block[BLOCK][BLOCK])
{
	int u, v, x, y;
	double CuCv;				// CuCv = Cu * Cv in the formala
	double SumRow, SumCol;
	double tempBlock[BLOCK][BLOCK];

	for (u = 0; u < BLOCK; u++)
	{
		for (v = 0; v < BLOCK; v++)
		{
			if ((u == 0) && (v == 0))		//Cu,Cv == 0
			{
				CuCv = 0.5;
			}
			else if ((u > 0) && (v > 0))	// Cu,Cv > 0
			{
				CuCv = 1.0;
			}
			else		
			{
				CuCv = 1 / sqrt(2);
			}
			SumRow = 0;
			for (x = 0; x < BLOCK; x++)
			{
				SumCol = 0;
				for (y = 0; y < BLOCK; y++)
				{
					SumCol = SumCol + (block[x][y] * DCTCOSINE(u, x) * DCTCOSINE(v, y));		// DCT transform formula
				}
				SumRow = SumRow + SumCol;
			}
			tempBlock[u][v] = (CuCv / 4.0) * SumRow;
		}
	}
	//copy the transformed values back to the block
	for (u = 0; u < BLOCK; u++)
	{
		for (v = 0; v < BLOCK; v++)
		{
			block[u][v] = tempBlock[u][v];
		}
	}
}
// zigzag reodering: convert the block to zigzag order
void dozigzag()
{
	int res[BLOCK*BLOCK];			// temp 1D array 8x8
	int j = 0;
	//begin zigzag order
	for( int i = 0; i < BLOCK*BLOCK; i++)
	{
		if(i % 2 == 1)
		{
			for(int y = min(i,BLOCK - 1); y >= max(0,i - (BLOCK -1)); y--)
				res[j++] = block[i-y][y];					// copy from block to 1D temp
		}
		else {
			for(int x = min(i,BLOCK - 1); x >= max(0,i- (BLOCK -1));x--)
				res[j++] = block[x][i-x];					// copy from block to 1D temp
		}
	}
	int k = 0;
	// copy result back to block
	for(int i = 0; i < BLOCK; ++i)
		for(int j = 0; j < BLOCK; ++j)
			block[i][j] = res[k++];							// copy from 1D temp to the block
}
// Write the block to file
void writetofile(int BlockBegin, int BlockEnd, FILE *output)
{
	fprintf(output, "%d %d\n", BlockBegin, BlockEnd);
	for (int s = 0; s < BLOCK; s++)
	{
		for (int t = 0; t < BLOCK; t++)
		{
			fprintf(output, "%5d", (unsigned char)block[s][t]);
		}
		fprintf(output, "\n");
	}
}
// handle data at Block level: calculate number of blocks, perform DCT formulat, quantize values and write each block to output file
void handleBlock(int YmBpos, int XmBpos, FILE *output)
{
	int XmB = 0, YmB = 0;
	int numMBlockX = MACROBLOCK / BLOCK;					// number of block per macroblock in col
	int numMBlockY = MACROBLOCK / BLOCK;					// number of block per one macroblock in row
	int BlockInitX, BlockInitY, BlockEndX, BlockEndY;
	for (int i = 0; i < numMBlockX; i++)
	{
		for (int j = 0; j < numMBlockY; j++)
		{
			BlockInitX = BLOCK * i;							// start position of block in col
			BlockEndX = BlockInitX + BLOCK;					// end position of block in col
			for (int m = BlockInitX; m < BlockEndX; m++)
			{
				BlockInitY = BLOCK * j;						// start position of block in row
				BlockEndY = BlockInitY + BLOCK;				// end position of block in row
				for (int n = BlockInitY; n < BlockEndY; n++)
				{
					block[XmB][YmB] = macroBlock[m][n];		// copy data to the block
					YmB++;
				}
				XmB++;
				YmB = 0;
			}
			calculateDCT(block);							// calculate the DCT transform 
			quantize(block, qMatrix, qscale);				// quantize the block
			dozigzag();										// transform the block to zigzag order
			writetofile(BlockInitY + (16 * YmBpos), BlockInitX + (16 * XmBpos), output);	// write each block to output file
			XmB = 0;
			YmB = 0;
		}
		XmB = 0;
		YmB = 0;
	}
}
// handle data at MacroBlock level: calculate possible number of macbroblock, split each into blocks 
void handleMacroBlock(unsigned char **pgm, int rowSize, int colSize, FILE *output)
{
	int XmB = 0, YmB = 0;

	int numMBlockX = colSize / MACROBLOCK;				// number of macroblock in col
	int numMBlockY = rowSize / MACROBLOCK;				// number of macroblock in row
	for (int i = 0; i < numMBlockX; i++)
	{
		for (int j = 0; j < numMBlockY; j++)
		{
			int MBlockInitX = (MACROBLOCK * i);				// start position of the macroblock in col
			int MBlockCapX = MBlockInitX + MACROBLOCK;		// end postion of the macroblock in col
			for (int m = MBlockInitX; m < MBlockCapX; m++)
			{
				int MBlockInitY = (MACROBLOCK * j);			// start position of the macroblock in row
				int MBlockCapY = MBlockInitY + MACROBLOCK;  // end postion of the macroblock in row
				for (int n = MBlockInitY; n < MBlockCapY; n++)
				{
					macroBlock[XmB][YmB] = pgm[m][n];		// copy data
					YmB++;
				}
				XmB++;
				YmB = 0;
			}
			handleBlock(j, i, output);						// handle problem at block level
			XmB = 0;
			YmB = 0;
		}
		XmB = 0;
		YmB = 0;
	}
}
// Main function
int main(int argc, char *argv[])
{
	int temp;
	temp = sscanf(argv[3], "%f", &qscale); 	// read quantization scale
	FILE *image = fopen(argv[1], "r"); 		// open pgm file
	FILE *qf = fopen(argv[2], "r"); 		// open quantize file
	FILE *output = fopen(argv[4], "w"); 	// open output file
	handleQuantFile(qf);
	handlePGMFile(image);
	// write the first 3 lines to output file
	fprintf(output, "%s\n", "MYDCT");
	fprintf(output, "%d %d\n", pgmData.rows, pgmData.cols);
	fprintf(output, "%f\n", qscale);
	// write the rest
	handleMacroBlock(pgmData.value, pgmData.rows, pgmData.cols, output);
	fclose(image);
	fclose(qf);
	fclose(output);
}
