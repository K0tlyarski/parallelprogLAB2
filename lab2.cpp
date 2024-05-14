
#include <iostream>
#include <fstream>
#include <tuple>
#include <chrono>
#include <omp.h>

std::string MATRIX1_PATH = "matrix1.txt";
std::string MATRIX2_PATH = "matrix2.txt";
std::string RESULT_PATH = "result.txt";
std::string STATS_PATH = "stats.txt";


struct stats {
	int num_threads;
	int size;
	double time;
};


void writeStats(stats st, std::string path) {
	std::ofstream file(path, std::ios_base::app);
	file << st.num_threads << " " << st.size << " " << st.time << std::endl;
	file.close();
}


void clearFile(std::string path) {
	std::ofstream file(path, std::ios_base::trunc);
	file.close();
}


int** generateMatrix(int rows, int cols, int low, int up) {
	srand(time(NULL));
	int** matrix = new int* [rows];
	for (int i = 0; i < rows; ++i) {
		matrix[i] = new int[cols];
	}
	for (int i = 0; i < rows; ++i) {
		for (int j = 0; j < cols; ++j) {
			matrix[i][j] = rand() % (up - low + 1) + low;
		}
	}
	return matrix;
}

void printMatrix(const int* const* matrix, int rows, int cols) {
	for (int i = 0; i < rows; ++i) {
		for (int j = 0; j < cols; ++j) {
			std::cout << matrix[i][j] << "  ";
		}
		std::cout << "\n";
	}
}

int** mulMatrix(const int* const* matrix1, const int* const* matrix2, int rows1, int cols1, int rows2, int cols2, int threads_num) {
	if (cols1 != rows2) {
		return NULL;
	}
	int** result = new int* [rows1];
	for (int i = 0; i < rows1; ++i) {
		result[i] = new int[cols2]();
	}
#pragma omp parallel num_threads(threads_num)
	{
#pragma omp for 
		for (int i = 0; i < rows1; ++i) {
			for (int j = 0; j < cols2; ++j) {
				for (int k = 0; k < rows2; ++k) {
					result[i][j] += matrix1[i][k] * matrix2[k][j];
				}
			}
		}
	}

	return result;
}

int** mulMatrix(const int* const* matrix1, const int* const* matrix2, int size, int threads_num) {
	int** result = new int* [size];
	for (int i = 0; i < size; ++i) {
		result[i] = new int[size]();
	}
#pragma omp parallel num_threads(threads_num) 
	{
#pragma omp for 
		for (int i = 0; i < size; ++i) {
			for (int j = 0; j < size; ++j) {
				for (int k = 0; k < size; ++k) {
					result[i][j] += matrix1[i][k] * matrix2[k][j];
				}
			}
		}
	}
	return result;
}


void writeMatrix(const int* const* matrix, int rows, int cols, const std::string path) {
	std::ofstream file(path);
	file << rows << " " << cols << std::endl;
	for (int i = 0; i < rows; ++i) {
		for (int j = 0; j < cols; ++j) {
			file << matrix[i][j] << " ";
		}
		file << std::endl;
	}
	file.close();
}

std::tuple<int**, int, int> readMatrix(std::string path) {
	std::ifstream file(path);
	int rows, cols;
	file >> rows >> cols;
	int** result = new int* [rows];
	for (int i = 0; i < rows; ++i) {
		result[i] = new int[cols]();
	}
	for (int i = 0; i < rows; ++i) {
		for (int j = 0; j < cols; ++j) {
			file >> result[i][j];
		}
	}
	file.close();
	return std::make_tuple(result, rows, cols);
}

int** mulMatrix(std::string mat1Path, std::string mat2Path, int threads_num) {
	std::tuple<int**, int, int> m1 = readMatrix(mat1Path);
	std::tuple<int**, int, int> m2 = readMatrix(mat2Path);
	int** mat1 = std::get<0>(m1);
	int** mat2 = std::get<0>(m2);
	int size = std::get<1>(m2);
	int** result = new int* [size];
	for (int i = 0; i < size; ++i) {
		result[i] = new int[size]();
	}
#pragma omp parallel num_threads(threads_num)
	{
#pragma omp for 
		for (int i = 0; i < size; ++i) {
			for (int j = 0; j < size; ++j) {
				for (int k = 0; k < size; ++k) {
					result[i][j] += mat1[i][k] * mat2[k][j];
				}
			}
		}
	}
	return result;
}

void deleteMatrix(int** matrix, int size) {
	for (int i = 0; i < size; ++i) {
		free(matrix[i]);
	}
	free(matrix);
}

int main() {
	clearFile(MATRIX1_PATH);
	clearFile(MATRIX2_PATH);
	clearFile(RESULT_PATH);
	clearFile(STATS_PATH);
	int size = 500;
	stats result;
	int p = 0;
	int tmp = 2;
	int range = 1000000;
	int** mat1 = generateMatrix(2, 2, -2, 2);
	int** mat2 = generateMatrix(2, 2, -2, 2);
	int** res = generateMatrix(2, 2, -2, 2);
	int max_threads = omp_get_max_threads();
	int thr_num = pow(tmp, p);
	while (thr_num <= max_threads) {
		while (size <= 700) {
			result.num_threads = thr_num;
			result.size = size;
			result.time = 0;
			std::cout << "SIZE: " << size << "\n";
			for (size_t i = 0; i < 10; ++i) {
				std::cout << "---------- Iteration " << i << " ----------\n";
				mat1 = generateMatrix(size, size, -range, range);
				mat2 = generateMatrix(size, size, -range, range);
				std::cout << "Matrices are generated" << "\n";
				auto start = std::chrono::high_resolution_clock::now();
				res = mulMatrix(mat1, mat2, size, thr_num);
				auto stop = std::chrono::high_resolution_clock::now();
				std::cout << "Result is calculated" << "\n";
				auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
				result.time += (double)duration.count();
				writeMatrix(mat1, size, size, MATRIX1_PATH);
				writeMatrix(mat2, size, size, MATRIX2_PATH);
				std::cout << "Matrices are written" << "\n";
				writeMatrix(res, size, size, RESULT_PATH);
				std::cout << "Result is written" << "\n";
				deleteMatrix(mat1, size);
				deleteMatrix(mat2, size);
				deleteMatrix(res, size);
			}
			result.time /= 10;
			std::cout << result.time << "\n";
			writeStats(result, STATS_PATH);
			std::cout << "Stats is written" << "\n\n\n";
			size += 100;
		}
		size = 500;
		p++;
		thr_num = pow(tmp, p);
	}
	return 0;
}
