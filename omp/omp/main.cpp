#include "stdafx.h"
#include <omp.h>

void ReadFile(const char *fileName) {
	float **matrixA;
	int n, m;

	FILE *file = fopen(fileName, "r");
	if (file == NULL)
		printf("File not found\nPress any key to end");
	else {
		fscanf(file, "%d %d", &n, &m);
		matrixA = new float*[n];
		for (int i = 0; i < n; i++){
			matrixA[i] = new float[m];
		}
		for (int i = 0; i < n; i++)
			for (int j = 0; j < m; j++){
				fscanf(file, "%f", &matrixA[i][j]);
			}
		fclose(file);

		for (int i = 0; i < n; i++){
			for (int j = 0; j < m; j++)
				printf("%f ", matrixA[i][j]);
			printf("\n");
		}
	}
}

long long fibonachi(long int n){
	if (n == 0 || n == 1)
		return 0;
	return fibonachi(n - 1) + fibonachi(n - 2);
}

int _tmain(int argc, _TCHAR* argv[]) {
	time_t begin, end;

	time(&begin);

	const char *fileName = "C:\\Sources\\itmo-c++labs\\omp\\files\\matrixA.in";
	ReadFile(fileName);

	time(&end);
	printf("\n%f\n", difftime(end, begin));

	system("pause");

	return 0;
}