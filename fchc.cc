#include <iostream>
#include <stdlib.h>
#include <time.h>
#include <omp.h>
#include <math.h>
#include <fstream>
#include <string>
#include <stdio.h>

#define C1 1
#define C2 (C1<<1)
#define C3 (C1<<2)
#define C4 (C1<<3)
#define C5 (C1<<4)
#define C6 (C1<<5)
#define C7 (C1<<6)
#define C8 (C1<<7)
#define C9 (C1<<8)
#define C10 (C1<<9)
#define C11 (C1<<10)
#define C12 (C1<<11)
#define C13 (C1<<12)
#define C14 (C1<<13)
#define C15 (C1<<14)
#define C16 (C1<<15)
#define C17 (C1<<16)
#define C18 (C1<<17)
#define C19 (C1<<18)
#define C20 (C1<<19)
#define C21 (C1<<20)
#define C22 (C1<<21)
#define C23 (C1<<22)
#define C24 (C1<<23)
#define OBS (C1<<24)

#define PI 3.14159265358979323846

using namespace std;

int C[24] = {
	C1,  C2,  C3,  C4,
	C5,  C6,  C7,  C8,
	C9,  C10, C11, C12,
	C13, C14, C15, C16,
	C17, C18, C19, C20,
	C21, C22, C23, C24
};

int Reverse[24] = {
	C4, C3, C2, C1,
	C8, C7, C6, C5,
	C12,C11,C10,C9,
	C16,C15,C14,C13,
	C20,C19,C18,C17,
	C24,C23,C22,C21
};

int CCELL_DIR[6][6] = {
	{ C1, C2, C5, C6, C9, C10 },
	{ C1, C3, C13, C14, C17, C18 },
	{ C5, C7, C13, C15, C21, C22 },
	{ C3, C4, C7, C8, C11, C12 },
	{ C2, C4, C15, C16, C19, C20 },
	{ C6, C8, C14, C16, C23, C24 }
};
const int c[24][4] = {
	{ 1,1,0,0 },
	{ 1,-1,0,0 },
	{ -1,1,0,0 },
	{ -1,-1,0,0 },
	
	{ 1,0,1,0 },
	{ 1,0,-1,0 },
	{ -1,0,1,0 },
	{ -1,0,-1,0 },

	{ 1,0,0,1 },
	{ 1,0,0,-1 },
	{ -1,0,0,1 },
	{ -1,0,0,-1 },

	{ 0,1,1,0 },
	{ 0,1,-1,0 },
	{ 0,-1,1,0 },
	{ 0,-1,-1,0 },

	{ 0,1,0,1 },
	{ 0,1,0,-1 },
	{ 0,-1,0,1 },
	{ 0,-1,0,-1 },

	{ 0,0,1,1 },
	{ 0,0,1,-1 },
	{ 0,0,-1,1 },
	{ 0,0,-1,-1 }
};

int switchBits(int n, int i, int j)
{
	int a, b;
	a = n & C[i];
	b = n & C[j];
	if (a > 0 && b == 0)
	{
		n ^= C[i];
		n |= C[j];
	}
	else if (b > 0 && a == 0)
	{
		n |= C[i];
		n ^= C[j];
	}
	return n;
}

int sigma(int n, int* q, int i, int j)
{
	int a = q[0];
	int b = q[1];
	int c = q[2];
	int d = q[3];

	if (i == 1)
	{
		q[0] = ( a + b + c - d ) / 2;
		q[1] = ( a + b - c + d ) / 2;
		q[2] = ( a - b + c + d ) / 2;
		q[3] = (-a + b + c + d ) / 2;

		n = switchBits(n, 1, 21);
		n = switchBits(n, 2, 22);
		n = switchBits(n, 5, 17);
		n = switchBits(n, 6, 18);
		n = switchBits(n, 8, 12);
		n = switchBits(n, 11, 15);
	}
	else
	{
		q[0] = ( a + b + c + d ) / 2;
		q[1] = ( a + b - c - d ) / 2;
		q[2] = ( a - b + c - d ) / 2;
		q[3] = ( a - b - c + d ) / 2;

		n = switchBits(n, 1, 20);
		n = switchBits(n, 2, 23);
		n = switchBits(n, 5, 16);
		n = switchBits(n, 6, 19);
		n = switchBits(n, 9, 12);
		n = switchBits(n, 10, 15);
	}
	return n;
}

int S(int n, int* q, int i, int j)
{
	switch (i)
	{
	case 1:
		q[0] = -q[0];
		n = switchBits(n, 0, 2);
		n = switchBits(n, 1, 3);
		n = switchBits(n, 4, 6);
		n = switchBits(n, 5, 7);
		n = switchBits(n, 8, 10);
		n = switchBits(n, 9, 11);
		break;
	case 2:
		q[1] = -q[1];
		n = switchBits(n, 0, 1);
		n = switchBits(n, 2, 3);
		n = switchBits(n, 12, 14);
		n = switchBits(n, 13, 15);
		n = switchBits(n, 16, 18);
		n = switchBits(n, 17, 19);
		break;
	case 3:
		q[2] = -q[2];
		n = switchBits(n, 4, 5);
		n = switchBits(n, 6, 7);
		n = switchBits(n, 12, 13);
		n = switchBits(n, 14, 15);
		n = switchBits(n, 20, 22);
		n = switchBits(n, 21, 23);
		break;
	case 4:
		q[3] = -q[3];
		n = switchBits(n, 8, 9);
		n = switchBits(n, 10, 11);
		n = switchBits(n, 16, 17);
		n = switchBits(n, 18, 19);
		n = switchBits(n, 20, 21);
		n = switchBits(n, 22, 23);
		break;
	default:
		break;
	}
	return n;
}

int P(int n, int*q, int i, int j)
{
	int a = i - 1;
	int b = j - 1;
	int qa = q[a];
	q[a] = q[b];
	q[b] = qa;

	if (i == 1)
	{
		if (j == 2)
		{
			n = switchBits(n, 1, 2);
			n = switchBits(n, 4, 12);
			n = switchBits(n, 5, 13);
			n = switchBits(n, 6, 14);
			n = switchBits(n, 7, 15);
			n = switchBits(n, 8, 16);
			n = switchBits(n, 9, 17);
			n = switchBits(n, 10, 18);
			n = switchBits(n, 11, 19);
		}
		else if (j == 3)
		{
			n = switchBits(n, 0, 12);
			n = switchBits(n, 1, 14);
			n = switchBits(n, 2, 13);
			n = switchBits(n, 3, 15);
			n = switchBits(n, 5, 6);
			n = switchBits(n, 8, 20);
			n = switchBits(n, 9, 21);
			n = switchBits(n, 10, 22);
			n = switchBits(n, 11, 23);
		}
		else if (j == 4)
		{
			n = switchBits(n, 0, 16);
			n = switchBits(n, 1, 18);
			n = switchBits(n, 2, 17);
			n = switchBits(n, 3, 19);
			n = switchBits(n, 4, 20);
			n = switchBits(n, 5, 22);
			n = switchBits(n, 6, 21);
			n = switchBits(n, 7, 23);
			n = switchBits(n, 9, 10);
		}
	}
	else if (i == 2)
	{
		if (j == 3)
		{
			n = switchBits(n, 0, 4);
			n = switchBits(n, 1, 5);
			n = switchBits(n, 2, 6);
			n = switchBits(n, 3, 7);
			n = switchBits(n, 13, 14);
			n = switchBits(n, 16, 20);
			n = switchBits(n, 17, 21);
			n = switchBits(n, 18, 22);
			n = switchBits(n, 19, 23);
		}
		else if (j == 4)
		{
			n = switchBits(n, 0, 8);
			n = switchBits(n, 1, 9);
			n = switchBits(n, 2, 10);
			n = switchBits(n, 3, 11);
			n = switchBits(n, 12, 20);
			n = switchBits(n, 13, 22);
			n = switchBits(n, 14, 21);
			n = switchBits(n, 15, 23);
			n = switchBits(n, 17, 18);
		}
	}
	else if (i == 3 && j == 4)
	{
		n = switchBits(n, 4, 8);
		n = switchBits(n, 5, 9);
		n = switchBits(n, 6, 10);
		n = switchBits(n, 7, 11);
		n = switchBits(n, 12, 16);
		n = switchBits(n, 13, 17);
		n = switchBits(n, 14, 18);
		n = switchBits(n, 15, 19);
		n = switchBits(n, 21, 22);
	}
	return n;
}

int* momenta(int node)
{
	int bit = 1;
	int*q = new int[4]{ 0,0,0,0 };

	for (int i = 0; i < 24; ++i)
	{
		if (bit & node)
		{
			for (int j = 0; j < 4; ++j)
			{
				q[j] += c[i][j];
			}
		}
		bit <<= 1;
	}
	return q;
}

int(*iso[])(int n, int*q, int i, int j) = { S,P,sigma };

void goBack(int**steps, int* nodes, int*q, int step, int length)
{
	for (int i = 1; i < length; i++)
	{
		for (int j = step - 1; j >= 0; --j)
		{
			nodes[i] = iso[steps[j][0]](nodes[i], q, steps[j][1], steps[j][2]);
		}
	}
}

int* newNode(int n, int**steps)
{
	int*q = momenta(n);

	int step = 0;
	
	// invert negative q[i]
	for (int i = 0; i < 4; i++)
	{
		if (q[i] < 0)
		{
			//S changes sign of q[i];
			n = S(n, q, i+1, 0);
			//we record each step to the steps field, since we will go back
			steps[step][0] = 0;
			steps[step][1] = i+1;
			steps[step][2] = 0;
			++step;
		}
	}
	// sort q[i] from high to low
	int highIndex;
	int highValue;
	for (int i = 0; i < 4; ++i)
	{
		highIndex = i;
		highValue = q[i];
		for (int j = i+1; j < 4; ++j)
		{
			if (q[j] > highValue)
			{
				highValue = q[j];
				highIndex = j;
			}
		}
		if (highIndex > i)
		{
			n = P(n, q, i + 1, highIndex + 1);
			steps[step][0] = 1;
			steps[step][1] = i + 1;
			steps[step][2] = highIndex + 1;
			++step;
		}
	}
	// to fulfill second condition:
	if (q[3] > 0)
	{
		if (q[0] + q[3] == q[1] + q[2])
		{
			n = sigma(n, q, 2, 0);
			steps[step][0] = 2;
			steps[step][1] = 2;
			steps[step][2] = 0;
			++step;
		}
		else if (q[0] + q[3] > q[1] + q[2])
		{
			n = sigma(n, q, 1, 0);
			steps[step][0] = 2;
			steps[step][1] = 1;
			steps[step][2] = 0;
			++step;
		}
	}
	if (q[3] < 0)
	{
		n = S(n, q, 4, 0);
		//we record each step to the steps field, since we will go back
		steps[step][0] = 0;
		steps[step][1] = 4;
		steps[step][2] = 0;
		++step;
	}


	int* nodes;
	//class 12
	if (q[0] == 0)
	{
		int length = 13;
		nodes = new int[length];
		nodes[0] = length - 1;
		//S3 S1 P34 P12,  S4 S1 P34 P12,  S3 S2 P34 P12, S4 S2 P34 P12,  
		//S2 S1 P24 P13,  S4 S1 P24 P13,  S3 S2 P24 P13, S4 S3 P24 P13,
		//S2 S1 P23 P14,  S3 S1 P23 P14,  S4 S2 P23 P14, S4 S3 P23 P14

		nodes[1] = P(n, q, 1, 2);
		nodes[1] = P(nodes[1], q, 3, 4);
		nodes[1] = S(nodes[1], q, 1, 0);
		nodes[1] = S(nodes[1], q, 3, 0);

		nodes[2] = P(n, q, 1, 2);
		nodes[2] = P(nodes[2], q, 3, 4);
		nodes[2] = S(nodes[2], q, 1, 0);
		nodes[2] = S(nodes[2], q, 4, 0);

		nodes[3] = P(n, q, 1, 2);
		nodes[3] = P(nodes[3], q, 3, 4);
		nodes[3] = S(nodes[3], q, 2, 0);
		nodes[3] = S(nodes[3], q, 3, 0);

		nodes[4] = P(n, q, 1, 2);
		nodes[4] = P(nodes[4], q, 3, 4);
		nodes[4] = S(nodes[4], q, 2, 0);
		nodes[4] = S(nodes[4], q, 4, 0);

		nodes[5] = P(n, q, 1, 3);
		nodes[5] = P(nodes[5], q, 2, 4);
		nodes[5] = S(nodes[5], q, 1, 0);
		nodes[5] = S(nodes[5], q, 2, 0);

		nodes[6] = P(n, q, 1, 3);
		nodes[6] = P(nodes[6], q, 2, 4);
		nodes[6] = S(nodes[6], q, 1, 0);
		nodes[6] = S(nodes[6], q, 4, 0);

		nodes[7] = P(n, q, 1, 3);
		nodes[7] = P(nodes[7], q, 2, 4);
		nodes[7] = S(nodes[7], q, 2, 0);
		nodes[7] = S(nodes[7], q, 3, 0);

		nodes[8] = P(n, q, 1, 3);
		nodes[8] = P(nodes[8], q, 2, 4);
		nodes[8] = S(nodes[8], q, 3, 0);
		nodes[8] = S(nodes[8], q, 4, 0);

		nodes[9] = P(n, q, 1, 4);
		nodes[9] = P(nodes[9], q, 2, 3);
		nodes[9] = S(nodes[9], q, 1, 0);
		nodes[9] = S(nodes[9], q, 2, 0);
		
		nodes[10] = P(n, q, 1, 4);
		nodes[10] = P(nodes[10], q, 2, 3);
		nodes[10] = S(nodes[10], q, 1, 0);
		nodes[10] = S(nodes[10], q, 3, 0);

		nodes[11] = P(n, q, 1, 4);
		nodes[11] = P(nodes[11], q, 2, 3);
		nodes[11] = S(nodes[11], q, 2, 0);
		nodes[11] = S(nodes[11], q, 4, 0);

		nodes[12] = P(n, q, 1, 4);
		nodes[12] = P(nodes[12], q, 2, 3);
		nodes[12] = S(nodes[12], q, 3, 0);
		nodes[12] = S(nodes[12], q, 4, 0);

		goBack(steps, nodes, q, step, length);
		
	}
	//class 11
	else if (q[1] == 0)
	{
		int length = 7;
		nodes = new int[length];
		nodes[0] = length - 1;
		//S4 S2 P23,   S4 S3 P23,   S3 S2 P24,   
		//S4 S3 P24,   S3 S2 P34,   S4 S2 P34

		nodes[1] = P(n, q, 2, 3);
		nodes[1] = S(nodes[1], q, 2, 0);
		nodes[1] = S(nodes[1], q, 4, 0);

		nodes[2] = P(n, q, 2, 3);
		nodes[2] = S(nodes[2], q, 3, 0);
		nodes[2] = S(nodes[2], q, 4, 0);

		nodes[3] = P(n, q, 2, 4);
		nodes[3] = S(nodes[3], q, 2, 0);
		nodes[3] = S(nodes[3], q, 3, 0);

		nodes[4] = P(n, q, 2, 4);
		nodes[4] = S(nodes[4], q, 3, 0);
		nodes[4] = S(nodes[4], q, 4, 0);

		nodes[5] = P(n, q, 3, 4);
		nodes[5] = S(nodes[5], q, 2, 0);
		nodes[5] = S(nodes[5], q, 3, 0);

		nodes[6] = P(n, q, 3, 4);
		nodes[6] = S(nodes[6], q, 2, 0);
		nodes[6] = S(nodes[6], q, 4, 0);
		

		goBack(steps, nodes, q, step, length);

	}
	//class 10
	else if (q[2] == 0 && q[0]==q[1])
	{
		int length = 7;
		nodes = new int[length];
		nodes[0] = length - 1;
		//S3,P34,P12,   S4,P34,P12,   S4,S3,sigma1,   S4,S3,P34,P12,sigma1
		//S4,S3,sigma2,   P34,P12,sigma2
		nodes[1] = P(n, q, 1, 2);
		nodes[1] = P(nodes[1], q, 3, 4);
		nodes[1] = S(nodes[1], q, 3, 0);

		nodes[2] = P(n, q, 1, 2);
		nodes[2] = P(nodes[2], q, 3, 4);
		nodes[2] = S(nodes[2], q, 4, 0);

		nodes[3] = sigma(n, q, 1, 0);
		nodes[3] = S(nodes[3], q, 3, 0);
		nodes[3] = S(nodes[3], q, 4, 0);

		nodes[4] = sigma(n, q, 1, 0);
		nodes[4] = P(nodes[4], q, 1, 2);
		nodes[4] = P(nodes[4], q, 3, 4);
		nodes[4] = S(nodes[4], q, 3, 0);
		nodes[4] = S(nodes[4], q, 4, 0);

		nodes[5] = sigma(n, q, 2, 0);
		nodes[5] = S(nodes[5], q, 3, 0);
		nodes[5] = S(nodes[5], q, 4, 0);

		nodes[6] = sigma(n, q, 2, 0);
		nodes[6] = P(nodes[6], q, 1, 2);
		nodes[6] = P(nodes[6], q, 3, 4);

		goBack(steps, nodes, q, step, length);
	}
	//class 9
	else if (q[2] == 0 && q[0] > q[1])
	{
		int length = 4;
		nodes = new int[length];
		nodes[0] = length - 1;
		
		nodes[1] = S(n, q, 3, 0);
		nodes[1] = S(nodes[1], q, 4, 0);
		
		nodes[2] = P(n, q, 3, 4);
		nodes[2] = S(nodes[2], q, 3, 0);
		
		nodes[3] = P(n, q, 3, 4);
		nodes[3] = S(nodes[3], q, 4, 0);

		goBack(steps, nodes, q, step, length);
	}
	//class 8
	else if (q[0]==q[1] && q[1] == q[2] && q[2] > q[3] && q[3]==0)
	{
		int length = 5;
		nodes = new int[length];
		nodes[0] = length - 1;
		//1
		nodes[1] = P(n, q, 1, 2);
		nodes[1] = P(nodes[1], q, 2, 3);
		//2
		nodes[2] = P(n, q, 1, 3);
		nodes[2] = P(nodes[2], q, 2, 3);
		//3
		nodes[3] = P(n, q, 1, 2);
		nodes[3] = P(nodes[3], q, 2, 3);
		nodes[3] = S(nodes[3], q, 4, 0);

		//4
		nodes[4] = P(n, q, 1, 3);
		nodes[4] = P(nodes[4], q, 2, 3);
		nodes[4] = S(nodes[4], q, 4, 0);

		goBack(steps, nodes, q, step, length);
	}
	//class 6 or 7
	else if (q[0]>q[1] && q[1]==q[2] && q[2] > q[3] && q[3] == 0)
	{
		//class 6
		if (q[0] == 2*q[1])
		{
			int length = 5;
			nodes = new int[length];

			nodes[0] = length - 1;

			nodes[1] = sigma(n, q, 1, 0);
			nodes[1] = S(nodes[1], q, 4, 0);

			nodes[2] = sigma(n, q, 2, 0);
			nodes[2] = S(nodes[2], q, 4, 0);

			nodes[3] = sigma(n, q, 1, 0);
			nodes[3] = P(nodes[3], q, 2, 3);
			nodes[3] = S(nodes[3], q, 4, 0);

			nodes[4] = sigma(n, q, 2, 0);
			nodes[4] = P(nodes[4], q, 2, 3);
			nodes[4] = S(nodes[4], q, 4, 0);

			goBack(steps, nodes, q, step, length);

		}
		//class 7
		else
		{
			int length = 2;
			nodes = new int[length];

			nodes[0] = length - 1;

			nodes[1] = P(n, q, 2, 3);
			nodes[1] = S(nodes[1], q, 4, 0);

			goBack(steps, nodes, q, step, length);
		}
		
	}
	//class 5
	else if (q[0]==q[1] && q[1]>q[2] && q[2] > q[3] && q[4] == 0)
	{
		int length = 2;
		nodes = new int[length];

		nodes[0] = length - 1;
		//1
		nodes[1] = P(n, q, 1, 2);
		nodes[1] = S(nodes[1], q, 4, 0);

		goBack(steps, nodes, q, step, length);
	}
	//class 3,4
	else if (q[0]>q[1] && q[1]>q[2] && q[2] > q[3] && q[3] == 0)
	{
		//class 3
		if (q[0]==q[1]+q[2])
		{
			int length = 3;
			nodes = new int[length];
			
			nodes[0] = length - 1;

			nodes[1] = sigma(n, q, 1, 0);
			nodes[1] = S(nodes[1], q, 4, 0);

			nodes[2] = sigma(n, q, 2, 0);
			nodes[2] = S(nodes[2], q, 4, 0);

			goBack(steps, nodes, q, step, length);

		}
		//class 4
		else
		{
			int length = 2;
			nodes = new int[length];

			nodes[0] = length - 1;

			nodes[1] = S(n, q, 4, 0);

			goBack(steps, nodes, q, step, length);
		}
	}
	//class 2
	else if (q[0]==q[1] && q[1]==q[2] && q[2]>q[3] && q[3]>0)
	{
		int length = 3;
		nodes = new int[length];

		nodes[0] = length - 1;

		nodes[1] = P(n, q, 1, 2);
		nodes[1] = P(nodes[1], q, 2, 3);

		nodes[2] = P(n, q, 1, 3);
		nodes[2] = P(nodes[2], q, 2, 3);

		goBack(steps, nodes, q, step, length);
	}
	//class 1
	else if (q[0]==q[1] && q[1]>q[2] && q[2] > q[3] && q[3] > 0)
	{
		int length = 2;
		nodes = new int[length];

		nodes[0] = length - 1;

		nodes[1] = P(n, q, 1, 2);

		goBack(steps, nodes, q, step, length);
	}
	else
		nodes = nullptr;
	return nodes;
}


void fillTable(int** table)
{
	int n;
	
	int ** steps = new int*[20];
	for(int i = 0; i<20; ++i)
	   steps[i] = new int[3];	

//#pragma omp parallel for private (n, steps)
	for ( n = 0; n < OBS; ++n)
	{
		table[n] = newNode(n, steps);
		if (!(n % 100000))
			cout << " n = " << n << endl;
	}
}

int PeriodicBC(int n, int N)
{
	if (n < 0)
		return N - 1;
	else if (n < N)
		return n;
	else
		return 0;
}

void Propagation(int***from, int***to, int**table, int X, int Y, int Z, int R = 1)
{
	int x,y,z;
	int x1, y1, z1;
	int x2, y2, z2;

	R = X / 2;
	int R2 = R*R;
	int i;
	int n;
	int new_x, new_y, new_z;

	bool stat = (R > 0);

#pragma omp parallel for private (n, new_x, new_y, new_z, x, y, z, x1, y1, z1, x2, y2, z2, i)
	for ( x = 1; x < X-1; ++x)
	{
		if (stat)
		{
			x1 = x - R;
			x2 = x1 * x1;
		}
		for ( y = 1; y < Y-1; ++y)
		{
			if (stat)
			{
				y1 = y - R;
				y2 = y1 * y1;
			}
			for ( z = 1; z < Z-1; ++z)
			{
				if (stat)
				{
					z1 = z - R;
					z2 = z1 * z1;
				}
				if (!stat || x2 + y2 + z2 < R2)
				{
				n = from[x][y][z];
				for ( i = 0; i < 24; ++i)
				{
					if (n & C[i])
					{
						new_x = PeriodicBC(x + c[i][0], X);
						new_y = PeriodicBC(y + c[i][1], Y);
						new_z = PeriodicBC(z + c[i][2], Z);
						if (to[new_x][new_y][new_z] & OBS)
							to[new_x][new_y][new_z] |= Reverse[i];
						else
							to[new_x][new_y][new_z] |= C[i];
					}
				}
				from[x][y][z] = 0;
				}
			}
		}
	}
}

void compute_velocity(int***grid, double****v, double****mean, double*****gamma, int X, int Y, int Z, int side, int I, int J, int K)
{
	double N = side*side*side;

	int i, j, k, l, m, n;

	int x, y, z;

#pragma omp parallel for private (i,j,k,l,m,n,x,y,z)

	for (i = 0; i < I; ++i)
	{
		for (j = 0; j < J; ++j)
		{
			for (k = 0; k < K; ++k)
			{
				for (x = i*side; x < (i + 1)*side; ++x)
				{
					for (y = j*side; y < (j + 1)*side; ++y)
					{
						for (z = k*side; z < (k + 1)*side; ++z)
						{
							n = grid[x][y][z];
							for (m = 0; m < 24; ++m)
							{
								if (n & C[m])
								{
									for (l = 0; l < 3; ++l)
									{
										v[i][j][k][l] += c[m][l];
									}
								}
							}
						}
					}
				}
				for (l = 0; l < 3; l++)
				{
					v[i][j][k][l] /= N;
					mean[i][j][k][l] += v[i][j][k][l];
				}
				for (l = 0; l < 3; l++)
				{
					for ( m = 0; m < 3; m++)
					{
						gamma[i][j][k][l][m] += (v[i][j][k][l] * v[i][j][k][m]);
					}
				}
			}
		}
	}
}

string Print(double****v, int t, int I, int J, int K, int side, string file_name)
{
	int half = side/2;
	int i,j,k,l;

	ofstream out;
	file_name += to_string(t) + ".dat";
	out.open(file_name);
	
	for (i = 0; i < I; ++i)
	{
		for (j = 0; j < J; ++j)
		{
			for (k = 0; k < K; ++k)
			{
				out << half + i*side << "\t" << half + j*side << "\t" << half + k*side;
				for (l = 0; l < 3; l++)
				{
					out << "\t" << v[i][j][k][l];
				}
				out << endl;
			}
		}
	}
	out.close();
	return file_name;
}

string Print(double*****g, int t, int I, int J, int K, int side, string file_name)
{
	int half = side / 2;
	int i, j, k, l, m;

	ofstream out;
	file_name += to_string(t) +".dat";
	out.open(file_name);

	for (i = 0; i < I; ++i)
	{
		for (j = 0; j < J; ++j)
		{
			for (k = 0; k < K; ++k)
			{
				out << half + i*side << "\t" << half + j*side << "\t" << half + k*side;
				for (l = 0; l < 3; l++)
				{
					for (m = 0; m < 3; m++)
					{
						out << "\t" << g[i][j][k][l][m];
					}
				}
				out << endl;
			}
		}
	}
	out.close();
	return file_name;
}

void Collision(int***grid, int**table, int t, int X, int Y, int Z, int R = 0)
{
	int x,y,z;
	int x1, y1, z1;
	int x2, y2, z2;

	int n;
	int k;
	
	R = X / 2;
	int R2 = R*R;

#pragma omp parallel for private (x,y,z,x1,y1,z1,x2,y2,z2,n,k)

	for (x = 0; x < X; ++x)
	{
		if (R > 0)
		{
			x1 = x - R;
			x2 = x1 * x1;
		}
		for (y = 0; y < Y; ++y)
		{
			if (R > 0)
			{
				y1 = y - R;
				y2 = y1 * y1;
			}
			for (z = 0; z < Z; ++z)
			{
				if (R > 0)
				{
					z1 = z - R;
					z2 = z1 * z1;
				}
				if (R == 0 || x2 + y2 + z2 < R2)
				{
					n = grid[x][y][z];
					if (n & OBS)
						break;
					k = (t % table[n][0]) + 1;
					grid[x][y][z] = table[n][k];
					++t;
				}
			}
		}
	}
}

/* ARRAY OF NODES - THE MICROSCOPIC GRID */
int***newGrid(int X, int Y, int Z)
{
	int***grid = new int**[X];
	for (int x = 0; x < X; x++)
	{
		grid[x] = new int*[Y];
		for (int y = 0; y < Y; y++)
		{
			grid[x][y] = new int[Z]();

		}
	}
	return grid;
}

#define Xflow (C1+C2+C5+C6+C9+C10)
/* PARTICLES WITH VELOCITY IN X-DIRECTION */
void set_initial(int***grid, int X, int Y, int Z, int R = 0)
{
	int x,y,z;
	

	x = 1;	
	if (R == 0)
	{
		for (z = 1; z < Z - 1; ++z)
			for (y = 1; y < Y - 1; ++y)
				grid[x][y][z] = Xflow;
	}
	else
	{
		int y1, z1;
		int Yh = Y / 2;
		int Zh = Z / 2;
		int R2 = R*R-1;

		for (y = 0; y < Y; y++)
		{
			y1 = y - Yh;
			for (z = 0; z < Z; z++)
			{
				z1 = z - Zh;
				if (y1*y1 + z1*z1 < R2)
				{
					grid[x][y][z] = Xflow;
				}
			}
		}
	}
}

/* allocate array for velocity */
double****newVelocity(int I, int J, int K)
{
	int i,j,k;
	double****v = new double***[I];
	for(i=0; i<I; ++i)
	{
		v[i] = new double**[J];
		for(j=0; j<J; ++j)
		{
			v[i][j] = new double*[K];
			for(k=0; k<K; ++k)
			{
				v[i][j][k] = new double[3]();
			}
		}
	}
	return v;
}

/* allocate array for covariance tensor */
double*****newGamma(int I, int J, int K)
{
	int i, j, k, l;
	double*****g = new double****[I];
	for (i = 0; i < I; ++i)
	{
		g[i] = new double***[J];
		for (j = 0; j < J; j++)
		{
			g[i][j] = new double**[K];
			for ( k = 0; k < K; k++)
			{
				g[i][j][k] = new double*[3];
				for (l = 0; l < 3; ++l)
				{
					g[i][j][k][l] = new double[3]();
				}
			}
		}
	}
	return g;
}

/**/
void plot(string file_name, int X, int Y, int Z)
{
	const char* f = file_name.c_str();
	FILE * pipe = popen("gnuplot -persistent", "w");

	fprintf(pipe, "reset\n");
	fprintf(pipe, "set terminal pngcairo\n");
	fprintf(pipe, "set output '%s.png'\n", f);
	fprintf(pipe, "set view equal xyz\n");
	fprintf(pipe, "set xlabel \"X\"\n");
	fprintf(pipe, "set ylabel \"Y\"\n");
	fprintf(pipe, "set zlabel \"Z\"\n");
	fprintf(pipe, "set xrange [-1:%i]\n",X);
	fprintf(pipe, "set yrange [-1:%i]\n",Y);
	fprintf(pipe, "set zrange [-1:%i]\n",Z);
	fprintf(pipe, "set ticslevel 0\n");
	fprintf(pipe, "splot \"%s\" with vectors\n", f);

	pclose(pipe);	
}



void set_round_plate_obstacle(int*** a, int R, int x, int X, int Y, int Z)
{
	int Sz = Z/2;
	int Sy = Y/2;

	int R2 = R*R;
	int z,y;
	int z_s, y_s;

	for(z=-Sz; z < Sz; ++z)
		for(y=-Sy; y < Sy; ++y)
		{
			if(z*z + y*y < R2)
				a[x][y+Sy][z+Sz] |= OBS;
		}
}

void set_sphere_obstacle(int*** a, int R, int Sx, int X, int Y, int Z)
{
	int Sy = Y / 2;
	int Sz = Z / 2;

	int R2 = R*R;

	int x, y, z;
	int x1, y1, z1;
	int x2, y2, z2;

	for (x = Sx - R; x < Sx + R; x++)
	{
		x1 = x - Sx;
		x2 = x1*x1;
		for (y = Sy - R; y < Sy + R; y++)
		{
			y1 = y - Sy;
			y2 = y1*y1;
			for (z = Sz - R; z < Sz + R; z++)
			{
				z1 = z - Sz;
				z2 = z1*z1;
				if (x2 + y2 + z2 < R2)
				{
					a[x][y][z] |= OBS;
				}
			}
		}
	}
}

//Xp,Yp,Xm,Ym,Zp,Zm
int CELL_DIR[6][6] = {
	{ C1, C2, C5, C6, C9, C10 },
	{ C1, C3, C13, C14, C17, C18 },
	{ C5, C7, C13, C15, C21, C22 },
	{ C3, C4, C7, C8, C11, C12 },
	{ C2, C4, C15, C16, C19, C20 },
	{ C6, C8, C14, C16, C23, C24 }
};

int flow_dir(int a, int b, int c = -1)
{
	int dir;
	int not_dir;
	int i, j, k;

	for (i = 0; i < 6; i++)
	{
		dir |= CELL_DIR[a][i];
		dir |= CELL_DIR[b][i];

		not_dir |= CELL_DIR[(a + 3) % 6][i];
		not_dir |= CELL_DIR[(b + 3) % 6][i];

		if (c > 0)
		{
			dir |= CELL_DIR[c][i];
			not_dir |= CELL_DIR[(c + 3) % 6][i];
		}

	}
	dir &= (~not_dir);
	return dir;
}

void set_initial_sphere(int***a, int X, int Y, int Z, int R = 0)
{
	R = X / 2;
	int R2out = R*R;
	int R2in = (R - 2)*(R - 2);
	int up = R*cos(PI / 4);
	int down = R*sin(PI / 4);
	
	int x, y, z;
	int x1, y1, z1;
	int x2, y2, z2;

	int i;
	int c;
#pragma omp parallel for private (x, y, z, x1, y1, z1, x2, y2, z2, i, c)
	for (x = 0; x < X; ++x)
	{
		x1 = x - R;
		x2 = x1*x1;
		for (y = 0; y < Y; ++y)
		{
			y1 = y - R;
			y2 = y1*y1;
			for (z = 0; z < Z; ++z)
			{
				z1 = z - R;
				z2 = z1*z1;
				if (x2 + y2 + z2 > R2in && x2 + y2 + z2 < R2out)
				{
					c = 0;
					if (z1 > up)
					{
						for (i = 0; i < 6; ++i)
						{
							c |= CELL_DIR[5][i];
						}
					}
					else if (z1 < -up)
					{
						for (i = 0; i < 6; i++)
						{
							c |= CELL_DIR[2][i];
						}
					}
					else if (y1 > up)
					{
						for (i = 0; i < 6; i++)
						{
							c |= CELL_DIR[4][i];
						}
					}
					else if (y1 < -up)
					{
						for (i = 0; i < 6; i++)
						{
							c |= CELL_DIR[1][i];
						}
					}
					else if (x1 > up)
					{
						for (i = 0; i < 6; i++)
						{
							c |= CELL_DIR[3][i];
						}
					}
					else if (x1 < -up)
					{
						for (i = 0; i < 6; i++)
						{
							c |= CELL_DIR[0][i];
						}
					}
					else if (x1 < down && x1 > -down)
					{
						c = flow_dir(y1 > 0 ? 4 : 1, z1 > 0 ? 5 : 2);
					}
					else if (y1 < down && y1 > -down)
					{
						c = flow_dir(x1 > 0 ? 3 : 0, z1 > 0 ? 5 : 2);
					}
					else if (z1 < down && z1 > -down)
					{
						c = flow_dir(x1 > 0 ? 3 : 0, y1 > 0 ? 4 : 1);
					}
					else
					{
						c = flow_dir(x1 > 0 ? 3 : 0, y1 > 0 ? 4 : 1, z1 > 0 ? 5 : 2);
					}
					a[x][y][z] = c;
				}
			}
		}
	}
}

void set_round_tunnel(int*** a, int R, int X, int Y, int Z)
{
	int x, y, z;

	int R2 = R*R - 1;

	int Yh = Y / 2;
	int Zh = Z / 2;
	int y2, z2;

	for (x = 0; x < X; ++x)
	{
		for (y = 0; y < Y; y++)
		{
			y2 = y - Yh;
			for (z = 0; z < Z; z++)
			{
				z2 = z - Zh;
				if (y2*y2 + z2*z2 >= R2)
				{
					a[x][y][z] |= OBS;
				}
			}
		}
	}
}

void set_square_tunnel(int***a, int X, int Y, int Z)
{
	int x,y,z;

	for(x=0; x < X; ++x)
	{
		for(z = 0; z < Z; ++z)
		{
			a[x][0][z] |= OBS;
			a[x][Y-1][z] |= OBS;
		}
		for(y = 0; y < Y; ++y)
		{
			a[x][y][0] |= OBS;
			a[x][y][Z-1] |= OBS;
		}
	}
}

void finalizeMean(double****mean, int side, int I, int J, int K, int div)
{
	int N = side*side*side;
	int i, j, k, l;
}

void compute_mean(double****mean, int I, int J, int K, int div)
{
	int i, j, k, l;

#pragma omp parallel for private (i, j, k, l)
	for ( i = 0; i < I; i++)
	{
		for ( j = 0; j < J; j++)
		{
			for ( k = 0; k < K; k++)
			{
				for ( l = 0; l < 3; l++)
				{
					mean[i][j][k][l] /= div;

				}
			}
		}
	}
}

void compute_covariance_tensor(double*****gamma, double****mean, int I, int J, int K, int div)
{
	int i, j, k, l, m;
	
#pragma omp parallel for private (i, j, k, l, m)
	for (i = 0; i < I; i++)
	{
		for (j = 0; j < J; j++)
		{
			for (k = 0; k < K; k++)
			{
				for (l = 0; l < 3; l++)
				{
					for ( m = 0; m < 3; m++)
					{
						gamma[i][j][k][l][m] = ( gamma[i][j][k][l][m] / div ) - ( mean[i][j][k][l] * mean [i][j][k][m] );
					}
				}
			}
		}
	}
}

void NullArray(double****v, int I, int J, int K)
{
	int i,j,k,l;
#pragma omp parallel for private (i, j, k, l)
	for(i = 0; i < I; ++i)
		for(j = 0; j < J; ++j)
			for(k = 0; k < K; ++k)
				for(l = 0; l < 3; ++l)
					v[i][j][k][l] = 0;
}

void NullArray(double*****g, int I, int J, int K)
{
	int i,j,k,l,m;
#pragma omp parallel for private (i, j, k, l, m)
	for(i = 0; i < I; ++i)
		for(j = 0; j < J; ++j)
			for(k = 0; k < K; ++k)
				for(l = 0; l < 3; ++l)
					for(m = 0; m < 3; ++m)
						g[i][j][k][l][m] = 0;
}

int main()
{
	
	int** table = new int*[OBS];	

	fillTable(table);
	
	time_t start = time(NULL);
	int T = 100;

	int X = 120;
	int Y = 120;
	int Z = 120;
	int side = 10;

	int I = X/side;
	int J = Y/side;
	int K = Z/side;

	// radius of tunnel
	int Rt = Y / 2;

	// radius of barrier
	int Rb = Y / 4;
	
	// X-coordinate of the barrier
	int Sx = X / 4;

	double****velocity = newVelocity(I, J, K);
	double****mean = newVelocity(I, J, K);

	double*****gamma = newGamma(I, J, K);


	int*** even = newGrid(X, Y, Z);
	int*** odd  = newGrid(X, Y, Z);


	/* INITIAL CONDITIONS */
	//exploding_cube(even,X,Y,Z);


	/* OBSTACLES */	
	//set_round_plate_obstacle(even,30,90, X,Y,Z);
	//set_round_plate_obstacle(odd,30,90, X,Y,Z);
	//set_sphere_obstacle(even, Rb, Sx, X, Y, Z);
	//set_sphere_obstacle(odd, Rb, Sx, X, Y, Z);


	/* TUNNEL AT THE SIDES MADE OF OBSTACLES */
	//set_round_tunnel(even, Rt, X, Y, Z);
	//set_round_tunnel(odd, Rt, X, Y, Z);
	//set_square_tunnel(even,X,Y,Z);
	//set_square_tunnel(odd,X,Y,Z);
	


	int nepar;
	string file_name;
	for(int t = 1; t <= T; ++t)
	{
		cout << "som v kroku " << t << endl;

		nepar = t&1;
		if(nepar)
		{
			set_initial_sphere(odd, X, Y, Z);
			Collision(odd, table, t, X, Y, Z);		
			Propagation(odd, even, table, X, Y, Z);
		}
		else
		{
			set_initial_sphere(even, X, Y, Z);
			Collision(even, table, t, X, Y, Z);
			Propagation(even, odd, table, X, Y, Z);
		}
		if(!(t%10))
		{
			compute_velocity(odd,velocity,mean,gamma,X,Y,Z,side,I,J,K);
	//		file_name = Print(velocity,t,I,J,K,side,"velocity");
	//		plot(file_name,X,Y,Z);
		}
		if(!(t%100))
			cout << time(NULL) - start << " seconds have past so far..." << endl;
		if (!(t % 1000) && t <= 5000)
		{
			int div = 100;
			compute_mean(mean, I, J, K, div);
			compute_covariance_tensor(gamma, mean, I, J, K, div);
			file_name = Print(mean, t, I, J, K, side, "mean_velocity");
			plot(file_name, X, Y, Z);
			file_name = Print(gamma, t, I, J, K, side, "covariance_tensor");
	
			NullArray(mean, I, J, K);
			NullArray(gamma, I, J, K);
		}
		if (t == 15000)
		{
			int div = 1000;
			compute_mean(mean, I, J, K, div);
			compute_covariance_tensor(gamma, mean, I, J, K, div);
			file_name = Print(mean, t, I, J, K, side, "mean_velocity");
			plot(file_name, X, Y, Z);
			file_name = Print(gamma, t, I, J, K, side, "covariance_tensor");
		}
	}
	time_t stop = time(NULL);

	cout << "trvalo to " << stop - start << " sekund." << endl;
	return 0;
}
