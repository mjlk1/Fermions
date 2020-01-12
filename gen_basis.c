/* Implementation of algorithm described in https://arxiv.org/abs/1908.02084 for two component Fermi system.

	Compilation: gcc gen_basis.c -o gen_basis -lm

	Input: H/R (harmonic trap/ring) nA nB energyMax 
	For single particle states in ring trap state k is written as 2k, whereas state -k is written as 2k-1.*/

#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>

#define MAX_BASIS_COUNT 5000000
#define SPECTRUM_COUNT 501

void writeState(const int nA, const int nB, int *p, int *stateA, int *stateB)
{
	for (int i=0;i<nA;i++)
		*(p+i) = *(stateA+i+1);

	for (int i=0;i<nB;i++)
		*(p+nA+i) = *(stateB+i+1);
}

void generateBasis(const int nA, const int nB, const double energies[], const int energyMax, int *p, int *basisCount)
{
	int stateA[nA+1];
	int stateB[nB+1];
	int indexA = nA, indexB = nB;
	double energyA, energyB;

	*basisCount = 0;

	for (int i=1;i<=nA;i++)
		stateA[i] = i-1;

	while (indexA>0)
	{
		energyA = 0.0;

		for (int i=1;i<=nA;i++)
			energyA += energies[stateA[i]];

		if (energyA<=energyMax)
		{	
			for (int i=1;i<=nB;i++)
				stateB[i] = i-1;

			indexB = nB;
			while (indexB>0)
			{
				energyB = 0.0;

				for (int i=1;i<=nB;i++)
					energyB += energies[stateB[i]];

				if (energyA+energyB<=energyMax)
				{
					writeState(nA,nB,p,stateA,stateB);
					p += nA+nB;
					*basisCount += 1;
					indexB = nB;
				}
				else
					indexB--;

				if (indexB>0)
				{
					stateB[indexB]++;
	
					for (int i=indexB+1;i<=nB;i++)
						stateB[i] = stateB[i-1]+1;
				}
			}
			indexA = nA;
		}
		else
			indexA--;

		if (indexA>0)
		{
			stateA[indexA]++;

			for (int i=indexA+1;i<=nA;i++)
				stateA[i] = stateA[i-1]+1;
		}
	}
}

char* concat(const char *s1, const char *s2)
{
    const size_t len1 = strlen(s1);
    const size_t len2 = strlen(s2);
    char *result = malloc(len1 + len2 + 1);
    memcpy(result, s1, len1);
    memcpy(result + len1, s2, len2 + 1);
    return result;
}

int main(int argc, const char *argv[])
{
	int nA, nB, energyMax;
	char problemType;

	scanf("%c %d %d %d", &problemType, &nA, &nB, &energyMax);

	int stride = nA+nB;
	int basisCount, basisCountTemp;

	double energies[SPECTRUM_COUNT];
	if (problemType=='H')
		for (int i=0;i<SPECTRUM_COUNT;i++)
			energies[i] = 0.5+i;
	else if (problemType=='R')
	{
		energies[0] = 0.0;
		int i = 1, j = 1;
		while (j<SPECTRUM_COUNT)
		{
			energies[j]=i*i;
			energies[j+1]=i*i;
			j += 2;
			i++;
		}
	}
	else
	{
		printf("%s", "Unknown problem type");
		return 1;
	}

	double groundStateEnergy = 0.0;

	for (int i=0;i<nA;i++)
		groundStateEnergy += energies[i];

	for (int i=0;i<nB;i++)
		groundStateEnergy += energies[i];

	int energyMin = (int)floor(groundStateEnergy);

	int rangeBasis[energyMax-energyMin];
	int *pBasis = (int*)malloc(MAX_BASIS_COUNT*sizeof(int)*stride);

	generateBasis(nA,nB,energies,energyMin,pBasis,&basisCount);
	int indexLast = stride*basisCount-1; // pBasis+indexLast points at the last generated element
	rangeBasis[0] = basisCount;

	bool x;
	int q;
	for (int e=energyMin+1;e<=energyMax;e++)
	{
		int *pTempBasis = (int*)malloc(MAX_BASIS_COUNT*sizeof(int)*stride);
		generateBasis(nA,nB,energies,e,pTempBasis,&basisCountTemp);

		for (int i=0;i<basisCountTemp;i++)
		{
			x = true;
			for (int j=0;j<basisCount;j++)
			{
				q = 0;
				for (int k=0;k<stride;k++)
					if (*(pTempBasis+i*stride+k)==*(pBasis+j*stride+k))
						q++;

				if (q == stride)
					x = false;
			}

			if (x)
			{
				for (int j=0;j<stride;j++)
					*(pBasis+indexLast+1+j) = *(pTempBasis+i*stride+j);
				indexLast += stride;
			}
		}
		free(pTempBasis);
		basisCount = basisCountTemp;
		rangeBasis[e-energyMin] = basisCount; 
	}

	char nAs[4];
	sprintf(nAs,"%d",nA);
	char nBs[4];
	sprintf(nBs,"%d",nB);
	char energyMaxs[10];
	sprintf(energyMaxs,"%d",energyMax);
	char problemTypes[4];
	sprintf(problemTypes,"%c",problemType);

	char filename1[1000] = "./BasisData/";

	strcpy(filename1,concat(filename1, problemTypes));
	strcpy(filename1,concat(filename1, "_nA="));
	strcpy(filename1,concat(filename1, nAs));
	strcpy(filename1,concat(filename1, "_nB="));
	strcpy(filename1,concat(filename1, nBs));
	strcpy(filename1,concat(filename1, "_Emax="));
	strcpy(filename1,concat(filename1, energyMaxs));
	strcpy(filename1,concat(filename1, ".basis"));

	FILE *fptr1;
	fptr1 = fopen(filename1,"w");

   	for (int i=0; i<basisCount;i++)
	{
		for (int j=0;j<stride;j++)
		{
			fprintf(fptr1,"%d", *(pBasis+i*stride+j));
			if (j<stride-1)
				fprintf(fptr1," ");
		}
		fprintf(fptr1,"\n");
	}

	fclose(fptr1);

	char filename2[1000] = "./BasisData/";
	
	strcpy(filename2,concat(filename2, problemTypes));
	strcpy(filename2,concat(filename2, "_nA="));
	strcpy(filename2,concat(filename2, nAs));
	strcpy(filename2,concat(filename2, "_nB="));
	strcpy(filename2,concat(filename2, nBs));
	strcpy(filename2,concat(filename2, "_Emax="));
	strcpy(filename2,concat(filename2, energyMaxs));
	strcpy(filename2,concat(filename2, ".basisrange"));

	FILE *fptr2;
	fptr2 = fopen(filename2,"w");

   	for (int i=0; i<=energyMax-energyMin;i++)
		fprintf(fptr2,"%d \n", rangeBasis[i]);

	fclose(fptr2);

	free(pBasis);

	return 0;
}