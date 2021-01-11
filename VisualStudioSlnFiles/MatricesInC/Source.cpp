#include "Matrix.h"
#include "Timer.h"
#include <iostream>
#include <string>
#include <random>
int main() {
		{
				Timer time{};
				std::srand((unsigned int)time().time_since_epoch().count());
				Matrix mattie;
				std::cout << (std::string)mattie << "\n";
				mattie.replace(0, 1);
				std::cout << mattie << "\n";
				mattie.multiply(0, 3.1415);
				std::cout << mattie << "\n";
				Matrix mathew({ {1,2,3}, {4,5,6} });
				std::cout << mathew << "\n";
				mathew = mathew.transpose();
				std::cout << mathew << "\n";
				std::cout << mathew * mattie << "\n";
				Matrix barthelomew = mathew.subMatrix(0, 2, 0, 2);
				std::cout << barthelomew << "\n";
				Matrix bigBertha({ {1,2,3},{4,5},{7,8,9} });
				std::cout << bigBertha << "\n";
				std::cout << bigBertha.determinant() << "\n";
				Matrix biggerBertha({ {6,4,42,7}, {15,98,12,14}, {64,8,21,114}, {43,92,56,22} });
				std::cout << biggerBertha.determinant() << "\n";
				std::vector<std::vector<double>> building{};
				for (int r = 0; r < 4; r++) {
						std::vector<double> gurdur{};
						for (int c = 0; c < 4; c++) {
								gurdur.push_back(std::rand() >> 14);
						}
						building.push_back(gurdur);
				}
				Matrix built(building);
				double det;
				{
						Timer timer{};
						det = built.determinant();
				}
				std::cout << built << "\n"
						<< "determinant is " << det << "\n";
				double det2;
				{
						Timer timer{};
						det2 = built.rowEchelonForm();
				}
				std::cout << "actually it's " << det2 << std::endl;
				std::cout << built << "\n";
				built.reduceFromRowEchelon();
				std::cout << built << "\n";
				Matrix b({ {4,1}, {2,3} });
				std::cout << b * b << std::endl;
		}
		std::string end;
		std::cin >> end;
		return 0;
}