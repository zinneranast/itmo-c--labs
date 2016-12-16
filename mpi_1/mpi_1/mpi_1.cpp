#define _CRT_SECURE_NO_WARNINGS
#pragma comment (lib, "C:\\Sources\\itmo-c++labs\\mpi_1\\msmpi.lib")
#include <Windows.h>
#include "mpi.h"
#include <stdio.h>
#include <fstream>
#include <vector>
#include <iostream>

#define ROOT (0)

using namespace std;

void GenerateAB(const char *fileName, int n) {
	ofstream file;

	try {
		file.open(fileName);
		file << n << " " << n + 1 << endl;
		for (int i = 0; i < n; i++){
			for (int j = 0; j < n + 1; j++) {
				if (i != j)
					file << (rand() % 200 - 100)*0.0001 << " ";
				else
					file << n << " ";
			}
			file << endl;
		}

	}
	catch (exception &e) {
		printf("Error: Invalid file path\n", fileName);
	}
}

void GenerateX(const char *fileName, int n) {
	ofstream file;

	try {
		file.open(fileName);
		file << n << endl;
		for (int i = 0; i < n; i++)
			file << (rand() % 200 - 100)*0.0001 << endl;
	}
	catch (exception &e) {
		printf("Error: Invalid file path\n", fileName);
	}
}


void ReadAB(const char *fileName, vector<double> &A, vector<double> &B) {
	int n, m;
	ifstream file;

	try {
		file.open(fileName);
		file >> n >> m;
		if (n < 0 || m != n + 1)
			throw new exception("Error: Invalid matrix dimensions\n");

		for (int i = 0; i < n; i++)
			for (int j = 0; j < m; j++) {
				if (j == m - 1) file >> B[i];
				else file >> A[i * n + j];
			}

		file.close();
	}
	catch (exception &e) {
		printf("Error: File %s does not exist\n", fileName);
	}
}

void ReadX(const char *fileName, vector<double> &X) {
	int n;
	ifstream file;

	try {
		file.open(fileName);
		file >> n;
		if (n < 0)
			throw new exception("Error: Invalid matrix dimensions\n");

		for (int i = 0; i < n; i++)
			file >> X[i];

		file.close();
	}
	catch (exception &e) {
		printf("Error: File %s does not exist\n", fileName);
	}
}

void WriteX(const char *fileName, vector<double> &X) {
	int n = X.size();
	ofstream file;

	try {
		file.open(fileName);
		file << n;

		for (int i = 0; i < n; i++)
			file << X[i] << endl;

		file.close();
	}
	catch (exception &e) {
		printf("Error: Cannot open or close file %s\n", fileName);
	}
}

void PJacobi(vector<double> &A, vector<double> &B, vector<double> &X, double eps) {
	int n = B.size();
	double norm, maxnorm;
	double start, finish;
	int rank, size, numberOfRows = 0;
	int *sendNumberOfRows, *displsRows, *sendNumberOfElements, *displsElements;

	MPI_Init(NULL, NULL);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	rank == ROOT ? numberOfRows = n / size + n % size : numberOfRows = n / size;

	MPI_Bcast(&numberOfRows, 1, MPI_INT, 0, MPI_COMM_WORLD);

	sendNumberOfRows = new int[size];
	displsRows = new int[size];
	displsRows[0] = 0;

	MPI_Allgather(&numberOfRows, 1, MPI_INT, sendNumberOfRows, 1, MPI_INT, MPI_COMM_WORLD);

	sendNumberOfElements = new int[size];
	displsElements = new int[size];
	int numberOfElements = numberOfRows*n;
	displsElements[0] = 0;

	MPI_Allgather(&numberOfElements, 1, MPI_INT, sendNumberOfElements, 1, MPI_INT, MPI_COMM_WORLD);

	for (int i = 1; i < size; i++) {
		displsRows[i] = displsRows[i - 1] + sendNumberOfRows[i - 1];
		displsElements[i] = displsElements[i - 1] + sendNumberOfElements[i - 1];
	}

	vector<double> newA(numberOfRows*n);
	vector<double> newB(numberOfRows);
	vector<double> tempx(numberOfRows);
	vector<double> newX(numberOfRows);

	MPI_Scatterv(B.data(), sendNumberOfRows, displsRows, MPI_DOUBLE, &newB.front(), numberOfRows, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Scatterv(X.data(), sendNumberOfRows, displsRows, MPI_DOUBLE, &newX.front(), numberOfRows, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Scatterv(A.data(), sendNumberOfElements, displsElements, MPI_DOUBLE, &newA.front(), numberOfRows*n, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	int shift = (rank == ROOT) ? numberOfRows*rank : (numberOfRows + n % size)*rank;

	start = MPI_Wtime();
	do {
		for (int i = 0; i < numberOfRows; i++) {
			tempx[i] = newB[i];
			for (int j = 0; j < n; j++) {
				if (i + shift != j) tempx[i] -= newA[i*n + j] * X[j];
			}
			tempx[i] /= newA[i*n + i + shift];
		}

		if (numberOfRows > 0) norm = fabs(newX[0] - tempx[0]);

		for (int i = 0; i < numberOfRows; i++) {
			if (fabs(newX[i] - tempx[i]) > norm) norm = fabs(newX[i] - tempx[i]);
			newX[i] = tempx[i];
		}

		MPI_Reduce(&norm, &maxnorm, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
		MPI_Bcast(&maxnorm, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Allgatherv(tempx.data(), numberOfRows, MPI_DOUBLE, X.data(), sendNumberOfRows, displsRows, MPI_DOUBLE, MPI_COMM_WORLD);

	} while (maxnorm > eps);

	finish = MPI_Wtime();
	if (rank == ROOT) printf("Time to solve equation: %f sec\n", finish - start);

	MPI_Finalize();
}

int main(int argc, char* argv[]) {

	const int n = 2000;
	double eps = 0.0001;

	GenerateAB("AB.in", n);
	GenerateX("X.in", n);

	vector<double> A(n*n);
	vector<double> B(n);
	vector<double> X(n);

	ReadAB("AB.in", A, B);
	ReadX("X.in", X);

	PJacobi(A, B, X, eps);
	WriteX("X.out", X);
}