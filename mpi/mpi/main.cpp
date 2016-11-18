#include <Windows.h>
#include <mpi.h>
#include <stdio.h>
#include <fstream>
#include <vector>
#include <iostream>
#include <chrono>


void read_syst(std::ifstream& input, long n, std::vector<double>& a, std::vector<double>& b) {
	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++) {
			input >> a[i*n + j];
			if (j == n - 1) input >> b[i];
		}
	return;
}

void read_vec(std::ifstream& input, long n, std::vector<double>& vec) {
	for (int i = 0; i < n; i++)
		input >> vec[i];
}

void print_vec(std::ofstream& output, std::vector<double> vec) {
	for (int i = 0; i < vec.size(); i++)
		output << vec[i] << std::endl;
}

boolean check_mat(std::vector<double> a, long n) {
	double sum = 0;

	//lets check for sufficient condition :-)

	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			if (i != j) sum += fabs(a[i*n + j]);
			if (i == j && a[i*n + j] == 0) return false;
		}
		if (a[i*n + i] < sum) return false;
		sum = 0;
	}

	return true;
}

void jacobi(std::vector<std::vector<double>>& a, std::vector<double>& b, std::vector<double>& x, double eps) {
	//ax=b !! 
	//vector x - начальное приближение
	long n = a.size();
	std::vector<double> tempx(n);
	double norm;

	do {
		for (int i = 0; i < n; i++) {
			tempx[i] = b[i];
			for (int j = 0; j < n; j++) {
				if (i != j)
					tempx[i] -= a[i][j] * x[j];
			}
			tempx[i] /= a[i][i];
		}
		norm = fabs(x[0] - tempx[0]);
		for (int i = 0; i < n; i++) {
			if (fabs(x[i] - tempx[i]) > norm)
				norm = fabs(x[i] - tempx[i]);
			x[i] = tempx[i];
		}
	} while (norm > eps);
}

void jacobi_p(std::vector<double>& a, std::vector<double>& b, std::vector<double>& x, double eps) {
	//ax=b !! 
	//vector x - начальное приближение

	MPI_Init(NULL, NULL);
	long n = b.size();
	double norm, maxnorm;
	int world_rank, row_num = 0, rows_left, world_size;
	MPI_Comm_size(MPI_COMM_WORLD, &world_size);
	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
	MPI_Status status;
	auto start = std::chrono::high_resolution_clock::now();
	rows_left = n%world_size;

	if (world_rank == 0) {
		row_num = n / world_size;
	}

	MPI_Bcast(&row_num, 1, MPI_INT, 0, MPI_COMM_WORLD);

	if (world_rank == 0) 
		row_num += rows_left;

	//preparation for gatherv and scatterv
	int *counts1 = (int*)malloc(world_size * sizeof(int));
	int *displs1 = (int*)malloc(world_size * sizeof(int));
	int *counts2 = (int*)malloc(world_size * sizeof(int));
	int *displs2 = (int*)malloc(world_size * sizeof(int));
	displs1[0] = 0;
	displs2[0] = 0;
	int num_el = row_num*n;
	MPI_Allgather(&row_num, 1, MPI_INT, counts1, 1, MPI_INT, MPI_COMM_WORLD);
	MPI_Allgather(&num_el, 1, MPI_INT, counts2, 1, MPI_INT, MPI_COMM_WORLD);
	for (int i = 1; i < world_size; i++) {
		displs1[i] = displs1[i - 1] + counts1[i - 1];
		displs2[i] = displs2[i - 1] + counts2[i - 1];
	}

	//std::vector<std::vector<double>> pa(row_num,std::vector<double>(n));
	std::vector<double> pa(row_num*n);
	std::vector<double> pb(row_num);
	std::vector<double> ptempx(row_num);
	std::vector<double> px(row_num);
	
	MPI_Scatterv(b.data(), counts1, displs1, MPI_DOUBLE, &pb.front(), row_num, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Scatterv(x.data(), counts1, displs1, MPI_DOUBLE, &px.front(), row_num, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Scatterv(a.data(), counts2, displs2, MPI_DOUBLE, &pa.front(), row_num*n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		
	int shift = (world_rank == 0) ? row_num*world_rank : (row_num + rows_left)*world_rank;

	do {
		for (int i = 0; i < row_num; i++) {
			ptempx[i] = pb[i];
			for (int j = 0; j < n; j++) {
				if (i + shift != j)
					ptempx[i] -= pa[i*n + j] * x[j];
			}
			ptempx[i] /= pa[i*n + i + shift];
			//pa[i*n+ i + world_rank*row_num] is an pa[i][i + world_rank*row_num] "converted" to one-dimensional vector.
		}

		if (row_num > 0)
			norm = fabs(px[0] - ptempx[0]);	

		for (int i = 0; i < row_num; i++) {
			if (fabs(px[i] - ptempx[i]) > norm)
				norm = fabs(px[i] - ptempx[i]);
			px[i] = ptempx[i];
		}
	
		MPI_Reduce(&norm, &maxnorm, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
		MPI_Bcast(&maxnorm, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Allgatherv(ptempx.data(), row_num, MPI_DOUBLE, x.data(), counts1, displs1, MPI_DOUBLE, MPI_COMM_WORLD);

	} while (maxnorm > eps);

	auto finish = std::chrono::high_resolution_clock::now();
	if (world_rank == 0) std::cout << "execution time " << double(std::chrono::duration_cast<std::chrono::milliseconds>(finish - start).count()) / 1000 << std::endl;

	MPI_Finalize();
}

int main(int argc, char* argv[]) {
	std::ifstream mat_input(argv[1]);
	std::ifstream pribl_input(argv[2]);
	std::ofstream res_output(argv[4]);

	if (!mat_input.is_open() || !pribl_input.is_open()) {
		std::cout << "files were not open" << std::endl;
		return -1;
	}

	double eps = atof(argv[3]);
	if (eps <= 0.0) {
		std::cout << "eps value was negative / zero / not valid" << std::endl;
		return -1;
	}

	long m, n;
	mat_input >> m >> n;
	if (n != m + 1) {
		std::cout << "n!=m+1" << std::endl;
		return -1;
	}

	std::vector<double> a(m*m);
	std::vector<double> b(m);
	std::vector<double> result(m);

	read_syst(mat_input, m, a, b);

	if (!check_mat(a, m)) {
		std::cout << "this system probably cannot be solved by Jacobi method" << std::endl;
		//	return -1;
	}

	read_vec(pribl_input, m, result);
	jacobi_p(a, b, result, eps);
	print_vec(res_output, result);


}