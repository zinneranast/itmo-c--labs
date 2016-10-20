#include "stdafx.h"
#include <omp.h>

vector<vector<int>> ReadFile(const char *fileName) {
	int n, m;
	ifstream file;

	try {
		file.open(fileName);
		file >> n >> m;
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

vector<vector<int>> MultiplyMatrixes(vector<vector<int>> matrixA, vector<vector<int>> matrixB){
	vector<vector<int>> matrixC(matrixA.size(), vector<int>(matrixB[0].size()));

	for (int i = 0; i < matrixA.size(); i++)
		for (int j = 0; j < matrixB[0].size(); j++)
			for (int k = 0; k < matrixA[0].size(); k++)
				matrixC[i][j] += matrixA[i][k] * matrixB[k][j];

	return matrixC;
}

int _tmain(int argc, _TCHAR* argv[]) {
	time_t begin, end;
	time(&begin);

	const char *fileName;
	fileName = "C:\\Sources\\itmo-c++labs\\omp\\files\\matrixA.in";
	vector<vector<int>> matrixA = ReadFile(fileName);
	fileName = "C:\\Sources\\itmo-c++labs\\omp\\files\\matrixB.in";
	vector<vector<int>> matrixB = ReadFile(fileName);

	printf("%d\n%d\n", matrixA.size(), matrixA[0].size());
	printf("%d\n%d\n", matrixB.size(), matrixB[0].size());

	if (matrixA[0].size() == matrixB.size()){
		vector<vector<int>> matrixC;
		matrixC = MultiplyMatrixes(matrixA, matrixB);

		for (int i = 0; i < matrixC.size(); i++){
			for (int j = 0; j < matrixC[0].size(); j++)
				cout << matrixC[i][j] << " ";
			cout << endl;
		}
	}

	time(&end);
	printf("%.2f\n", difftime(end, begin));
	system("pause");

	return 0;
}