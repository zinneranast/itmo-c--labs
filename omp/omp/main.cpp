#include "stdafx.h"
#include <omp.h>

void GenerateMatrix(const char *fileName, int n, int m) {
	ofstream file;

	try {
		file.open(fileName);
		file << n << " " << m << endl;
		for (int i = 0; i < n; i++){
			for (int j = 0; j < m; j++)
				file << rand() % 200 - 100 << " ";
			file << endl;
		}

	} catch (exception &e) {
		printf("Error: Invalid file path\n", fileName);
	}
}

vector<vector<int>> ReadFile(const char *fileName) {
	int n, m;
	ifstream file;

	try {
		file.open(fileName);
		file >> n >> m;
		if (n < 0 || m < 0)
			throw new exception("Error: Invalid matrix dimensions\n");

		vector<vector<int>> matrix(n, vector<int>(m));

		for (int i = 0; i < n; i++)
			for (int j = 0; j < m; j++)
				file >> matrix[i][j];

		file.close();
		return matrix;

	}
	catch (exception &e) {
		printf("Error: File %s does not exist\n", fileName);
	}
}

void WriteFile(const char *fileName, vector<vector<int>> matrix) {
	int n, m;
	ofstream file;

	try {
		file.open(fileName);
		for (int i = 0; i < matrix.size(); i++){
			for (int j = 0; j < matrix[0].size(); j++)
				file << matrix[i][j] << " ";
			file << endl;
		}

	}
	catch (exception &e) {
		printf("Error: Invalid file path\n", fileName);
	}
}

int _tmain(int argc, _TCHAR* argv[]) {

	const char *fileName;

	fileName = "C:\\Sources\\itmo-c++labs\\omp\\files\\matrixA.in";
	GenerateMatrix(fileName, 1000, 1);
	vector<vector<int>> matrixA = ReadFile(fileName);

	fileName = "C:\\Sources\\itmo-c++labs\\omp\\files\\matrixB.in";
	GenerateMatrix(fileName, 1, 1000);
	vector<vector<int>> matrixB = ReadFile(fileName);
	
	int N = matrixA.size();
	int M = matrixB[0].size();
	int K = matrixA[0].size();

	if (K == matrixB.size()){
		vector<vector<int>> matrixC(N, vector<int>(M));
		int i, j, k, ijk;

		#ifdef _OPENMP

		double start_time, end_time;
		start_time = omp_get_wtime();

		#pragma omp parallel
		{
			#pragma omp single for private(i, j, k, ijk) schedule(static)
			//#pragma omp for private(i, j, k, ijk) schedule(dynamic)
			//#pragma omp for private(i, j, k, ijk) schedule(guided, 4)
			for (int ijk = 0; ijk < N*M*K; ijk++) {
				i = ijk / (M*K);
				j = ijk % M;
				k = ijk % K;
				matrixC[i][j] += matrixA[i][k] * matrixB[k][j];
				cout << matrixC[i][j] << endl;
			}
		}

		end_time = omp_get_wtime();
		printf("%.2f", end_time - start_time);

		#else

		time_t start_time, end_time;
		time(&start_time);

		for (int ijk = 0; ijk < N*M*K; ijk++) {
			i = ijk / (M*K);
			j = ijk % M;
			k = ijk % K;
			matrixC[i][j] += matrixA[i][k] * matrixB[k][j];
			cout << matrixC[i][j] << endl;
		}

		time(&end_time);
		printf("%.2f\n", difftime(end_time, start_time));
		
		#endif

		//WriteFile("C:\\Sources\\itmo-c++labs\\omp\\files\\matrixC.out", matrixC);
		system("pause");
		
		return 0;
	}
}