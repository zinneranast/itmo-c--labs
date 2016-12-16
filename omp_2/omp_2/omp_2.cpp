#include <omp.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <utility> 
#include <time.h>
#include <algorithm>

using namespace std;

const int INF = INT_MAX;

void GenerateAM(const char *fileName, int n) {
	ofstream file;

	try {
		file.open(fileName);
		file << n << endl;

		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				if (i == j) file << 0 << " ";
				else file << rand() % 10 << " ";
			}
			file << endl;
		}
	}
	catch (exception &e) {
		printf("Error: Invalid file path\n", fileName);
	}
}

void ReadAM(const char *fileName, vector<vector<pair<int, int>>> &M) {
	int n, w;
	ifstream file;

	try {
		file.open(fileName);
		file >> n;
		if (n < 0)
			throw new exception("Error: Invalid matrix dimensions\n");

		for (int i = 0; i < n; i++)
			for (int j = 0; j < n; j++) {
				file >> w;
				if (w > 0) M[i].push_back(make_pair(j, w));
			}

		file.close();
	}
	catch (exception &e) {
		printf("Error: File %s does not exist\n", fileName);
	}
}

void WriteWM(const char *fileName, vector<int> &M) {
	int n = M.size();
	ofstream file;

	try {
		file.open(fileName);
		file << n;

		for (int i = 0; i < n; i++)
			M[i] != INF ? file << M[i] << " " : file << "INF ";

		file.close();
	}
	catch (exception &e) {
		printf("Error: Cannot open or close file %s\n", fileName);
	}
}

void PDijkstra(vector<vector<pair<int, int>>> &WM, vector<int> &AM) {

	int n = AM.size();
	vector<int> parents(n);
	vector<int> processed(n, 0);
	double start, finish;

	AM[0] = 0;
	int v = -1, local_v, thread_num;
	int max_thread = omp_get_max_threads();

	#ifdef _OPENMP

	start = omp_get_wtime();

	for (int i = 0; i < n; i++) {
		vector<int> local_mins(max_thread, -1);

		#pragma omp parallel private(local_v, thread_num) 
		{
			thread_num = omp_get_thread_num();
			local_v = -1;

			#pragma omp for schedule(static)
			for (int j = 0; j < n; j++)
				if (!processed[j] && (local_v == -1 || AM[j] < AM[local_v])) {
					local_v = j;
					local_mins[thread_num] = local_v;
				}

		}

		int j = 0;
		while (local_mins[j] == -1) j++;

		v = local_mins[j];

		for (j; j < max_thread; j++) {
			if (AM[local_mins[j]] < AM[v] && local_mins[j] != -1) v = local_mins[j];
		}

		if (AM[v] == INF) {
			break;
		}

		processed[v] = 1;

		#pragma omp parallel for schedule(static)
		for (int j = 0; j < WM[v].size(); j++) {
			int to = WM[v][j].first,
				len = WM[v][j].second;
			if (AM[v] + len < AM[to]) {
				AM[to] = AM[v] + len;
			}
		}
	}

	finish = omp_get_wtime();

	#endif

	printf("Time: %f", finish - start);
}

int main(int argc, char* argv[]) {

	const int n = 300;

	GenerateAM("M.in", n);

	vector<int> AM(n, INF);
	vector<vector<pair<int, int>>> WM(n);

	ReadAM("M.in", WM);

	PDijkstra(WM, AM);
	WriteWM("M.out", AM);

	return 0;
}
