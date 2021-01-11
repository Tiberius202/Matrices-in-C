#pragma once
#include <vector>
#include <string>
#include <cstdarg>
#include <algorithm>
#include <numeric>
#define DEBUG 0
#if DEBUG
#include <iostream>
#define LOG(x) std::cout << x << std::endl
#else
#include <ostream>
#define LOG(x)
#endif
class Matrix {
private:
		//the only data that this class uses. 
		//Uses std::vector so matrix size can be inputted during runtime.
		std::vector<std::vector<double>> matrix;
public:
		Matrix() :
				matrix({ {1, 0}, {0, 1} })
		{
				LOG("default const");
		}
		Matrix(const size_t rows, const size_t columns) :
				matrix(rows, std::vector<double>(columns, 0))
		{
				LOG("identity constructed");
				for (size_t r = 0; r < matrix.size(); r++) {
						matrix[r][r] = 1;
				}
		}
		Matrix(const std::vector<std::vector<double>> mat) :
				matrix(mat)
		{
				LOG("manual constructed");
				size_t largestRow = 0;
				for (std::vector<double> row : mat) {
						largestRow = std::max(row.size(), largestRow);
				}
				/*std::for_each(mat.begin(), mat.end(),
					[this, largestRow](std::vector<double> curr, size_t col) {
						curr.resize(largestRow, 0);
					});*/
				for (size_t i = 0; i < matrix.size(); i++) {
						matrix[i].resize(largestRow, 0.0);
				}
		}
		Matrix(const Matrix& mat) {
				LOG("copied");
				matrix = mat.matrix;
		};
		Matrix(Matrix&& mat) noexcept {
				LOG("moved");
				matrix = std::move(mat.matrix);
		};
		Matrix& operator=(const Matrix& other) {
				LOG("assignment coperator");
				this->matrix = other.matrix;
				return *this;
		};
		Matrix& operator=(Matrix&& other) noexcept {
				LOG("move assignement");
				this->matrix = std::move(other.matrix);
				return *this;
		}
		~Matrix() {}
		//this is nothrow as long as allocators are the same???
		friend void swap(Matrix& first, Matrix& second) {
				std::swap(first.matrix, second.matrix);
		}

		size_t rows() const {
				return matrix.size();
		}
		size_t cols() const {
				return matrix[0].size();
		}
		void replace(const size_t sourceRow, const size_t recieveRow) {
				matrix[sourceRow].swap(matrix[recieveRow]);
		}
		void multiply(const size_t targetRow, const double scalarMultiple) {
				std::transform(matrix[targetRow].begin(), matrix[targetRow].end(),
						matrix[targetRow].begin(), [=](double x) {return x * scalarMultiple; });
		}
		static std::vector<double> multiply(std::vector<double> targetRow, const double scalarMultiple) {
				std::transform(targetRow.begin(), targetRow.end(),
						targetRow.begin(), [=](double x) {return x * scalarMultiple; });
				return targetRow;
		}
		void addRowToRow(const size_t sourceRow, const size_t recieveRow) {
				std::transform(matrix[sourceRow].begin(), matrix[sourceRow].end(),
						matrix[recieveRow].begin(), matrix[recieveRow].begin(), std::plus<double>());
		}
		static std::vector<double> addRowToRow(const std::vector<double> sourceRow, std::vector<double> recieveRow) {
				std::transform(sourceRow.begin(), sourceRow.end(),
						recieveRow.begin(), recieveRow.begin(), std::plus<double>());
				return recieveRow;
		}
		Matrix transpose() const {
				std::vector<std::vector<double>> newDimensions(matrix[0].size());
				for (size_t newRow = 0; newRow < matrix[0].size(); newRow++) {
						newDimensions[newRow].reserve(matrix.size());
						for (size_t newCol = 0; newCol < matrix.size(); newCol++) {
								newDimensions[newRow].push_back(matrix[newCol][newRow]);
						}
				}
				Matrix newMatrix(newDimensions);
				return newMatrix;
		}

		Matrix subMatrix(size_t top, size_t bottom, size_t left, size_t right) const {
				if (top > bottom) { std::swap(top, bottom); }
				//so left > right is going to have a wierd rollover effect to accomadate laplace expansion.
				if (right > left) {
						Matrix ret(bottom - top, right - left);
						for (size_t row = 0; top + row < bottom; row++) {
								for (size_t col = 0; left + col < right; col++) {
										ret.matrix[row][col] = matrix[top + row][left + col];
								}
						}
						return ret;
				}
				else {
						Matrix ret(bottom - top, right - left + this->cols());
						for (size_t row = 0; row < ret.rows(); row++) {
								for (size_t col = 0; col < right; col++) {
										ret.matrix[row][col] = matrix[top + row][col];
								}
								for (size_t col = 0; col < ret.cols() - right; col++) {
										ret.matrix[row][right + col] = matrix[top + row][left + col];
								}
						}
						return ret;
				}
		}
		double determinant() const {
				const size_t rows = matrix.size();
				if (rows != matrix[0].size()) {
						throw("Matrix must be a square matrix");
				}
				return rowEchelonForm(*this);
		}
		static double laPlaceExp(Matrix mat) {
				LOG(mat);
				if (mat.rows() == 1) {
						return mat.matrix[0][0];
				}
				double laPlaceRet = mat.matrix[0][0]
						* laPlaceExp(mat.subMatrix(1, mat.rows(), 1, mat.cols()));
				for (size_t row1 = 1; row1 < mat.cols() - 1; row1++) {
						LOG("bruh");
						laPlaceRet += (1 - (long long(row1) % 2) * 2) * mat.matrix[0][row1]
								* laPlaceExp(mat.subMatrix(1, mat.rows(), row1 + 1, row1)
								);
				}
				laPlaceRet += (1 - (long long(mat.cols() - 1) % 2) * 2) * mat.matrix[0][mat.cols() - 1]
						* laPlaceExp(mat.subMatrix(1, mat.rows(), 0, mat.cols() - 1)
						);
				return laPlaceRet;
		}
		//throws with vectors of different sizes
		static std::vector<double> multiplyTerms(std::vector<double> row1, std::vector<double> row2) {
				if (row1.size() != row2.size()) {
						throw std::invalid_argument("Trying to multiplyterms of different size vectors.");
				}
				std::transform(&row1[0], &row1[row1.size() - 1], &row2[0], &row1[0],
						std::multiplies<double>());
				return row1;
		}
		static double dotProduct(std::vector<double> row1, std::vector<double> row2) {
				if (row1.size() != row2.size()) {
						throw std::invalid_argument("Trying to dotproduct of different size vectors.");
				}
				return std::inner_product(row1.begin(), row1.end(), row2.begin(), 0.0);
		}
		//sign of determinant may be wrong
		double rowEchelonForm() {
				bool negative = false;
				if (rows() > 1) {
						auto partPlace = matrix.begin();
						//This is all to do the std::partition but to add the sign to the determinant. Very broken ATM.
						{
								auto backwardIterator = (matrix.end() - 1);
								while (partPlace < backwardIterator) {
										if ((*partPlace)[0] == 0.0) {
												std::swap(*partPlace, *(backwardIterator));
												negative = !negative;
												backwardIterator--;
												LOG(*this);
										}
										else {
												partPlace++;
												LOG(*this);
										}
								}
								if ((*partPlace)[0] != 0.0) {
										partPlace++;
										LOG(*this);
								}
								size_t distOver2 = (matrix.end() - partPlace) / 2;
								for (size_t i = 0; i < distOver2; i++) {
										std::swap(*(partPlace + i), *(matrix.end() - 1 - i));
										negative = !negative;
								}
						}
						LOG(*this);
						if (partPlace != matrix.begin()) {
								for (auto it = ++(matrix.begin()); it != partPlace; it++) {
										*it = addRowToRow(multiply(matrix[0], -1.0 * (*it)[0] / matrix[0][0]), *it);
								}
						}
						std::stable_sort(matrix.begin(), matrix.end(),
								[](std::vector<double> first, std::vector<double> second) {
										return std::find_if(first.begin(), first.end(), [](double x) {return x != 0.0; }) - first.begin()
												< std::find_if(second.begin(), second.end(), [](double x) {return x != 0.0; }) - second.begin();
								});
						Matrix small = subMatrix(1, rows(), 1, cols());
						small.rowEchelonForm();
						for (size_t smallRows = 0; smallRows < small.rows(); smallRows++) {
								for (size_t smallCols = 0; smallCols < small.cols(); smallCols++) {
										matrix[1 + smallRows][1 + smallCols] = small.matrix[smallRows][smallCols];
								}
						}
				}
				double determinant = matrix[0][0];
				for (size_t diagonalIterator = 1;
						diagonalIterator < std::min(cols(), rows());
						diagonalIterator++) {
						determinant *= matrix[diagonalIterator][diagonalIterator];
				}
				return ((1.0 - 2.0 * negative) * determinant);
		}
		static double rowEchelonForm(Matrix mat) {
				return mat.rowEchelonForm();
		}
		void reduceFromRowEchelon() {
				for (size_t revIt = matrix.size() - 1; revIt > 0; revIt--) {
						auto pivot = std::find_if(matrix[revIt].begin(), matrix[revIt].end(),
								[](double x) {return x != 0.0; });
						if (pivot != matrix[revIt].end()) {
								LOG(*pivot);
								multiply(revIt, 1.0 / (*pivot));
								for (size_t i = 0; i < revIt; i++) {
										matrix[i] = addRowToRow(multiply(matrix[revIt], -1.0 * (matrix[i][pivot - matrix[revIt].begin()])), matrix[i]);
								}
						}
				}
				auto pivot = std::find_if(matrix[0].begin(), matrix[0].end(),
						[](double x) {return x != 0.0; });
				if (pivot != matrix[0].end()) {
						multiply(0, 1.0 / (*pivot));
				}
				std::stable_sort(matrix.begin(), matrix.end(),
						[](std::vector<double> first, std::vector<double> second) {
								return std::find_if(first.begin(), first.end(), [](double x) {return x != 0.0; }) - first.begin()
										< std::find_if(second.begin(), second.end(), [](double x) {return x != 0.0; }) - second.begin();
						});
		}

		Matrix operator+(Matrix otherMatrix) const {
				if (matrix.size() != otherMatrix.matrix.size() || matrix[0].size() != matrix[0].size()) {
						throw std::invalid_argument("Matricies of different sizes");
				}
				for (size_t row = 0; row < matrix.size(); row++) {
						std::transform(otherMatrix.matrix[row].begin(), otherMatrix.matrix[row].end(),
								matrix[row].begin(), otherMatrix.matrix[row].begin(), std::plus<double>());
				}
				return otherMatrix;
		}
		void operator+=(Matrix otherMatrix) {
				if (matrix.size() != otherMatrix.matrix.size() || matrix[0].size() != matrix[0].size()) {
						throw std::invalid_argument("Matricies of different sizes");
				}
				for (size_t row = 0; row < matrix.size(); row++) {
						std::transform(matrix[row].begin(), matrix[row].end(),
								otherMatrix.matrix[row].begin(), matrix[row].begin(), std::plus<double>());
				}
		}
		//makes a friggin copy for no good reason. I will not pass by reference even if it kills me
		Matrix operator*(const Matrix otherMatrix) const {
				if (matrix[0].size() != otherMatrix.matrix.size()) {
						throw std::invalid_argument("Dimension mismatch (the columns of the first do not equal the rows of the second)");
				}
				std::vector<std::vector<double>> resultant(matrix.size());
				Matrix transpose = otherMatrix.transpose();
				for (size_t newRow = 0; newRow < matrix.size(); newRow++) {
						resultant[newRow].reserve(otherMatrix.matrix[0].size());
						for (size_t newCol = 0; newCol < otherMatrix.matrix[0].size(); newCol++) {
								resultant[newRow].push_back(dotProduct(matrix[newRow], transpose.matrix[newCol]));
						}
				}
				return resultant;
		}
		operator std::vector<std::vector<double>>() const {
				return matrix;
		}
		bool operator==(std::vector<std::vector<double>> otherMatrix) const {
				if (matrix.size() != otherMatrix.size()) {
						return false;
				}
				for (size_t r = 0; r < matrix.size(); r++) {
						if (matrix[r].size() != otherMatrix[r].size()) {
								return false;
						}
						for (size_t c = 0; c < matrix[r].size(); c++) {
								if (matrix[r][c] != otherMatrix[r][c]) {
										return false;
								}
						}
				}
				return true;
		}
		bool operator!=(std::vector<std::vector<double>> otherMatrix) const {
				return !(operator==(otherMatrix));
		}
		friend std::ostream& operator<<(std::ostream& strm, const Matrix& mat) {
				return strm << (std::string)mat;
		}
		operator std::string() const {
				std::string retValue = "\n";
				for (std::vector<double> rowVector : matrix) {
						retValue += "{ ";
						for (double term : rowVector) {
								retValue += std::to_string(term) + " ";
						}
						retValue += "}\n";
				}
				return retValue;
		}
};