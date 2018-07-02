#include "Matrix.h"


int main() {
	
	Matrix<double> obj1(ROWS, COLUMNS);
	obj1.IenterElements(x, b, c);
	obj1.ItoDual();
	obj1.Isimplex();

	cout << "------------------------------------------\n";
	cout << "------------------------------------------\n";
	cout << "------------------------------------------\n";

	Matrix<double> obj2(ROWS, COLUMNS);
	obj2.IenterElements(x, b, c);
	obj2.Isimplex();
	
	cout << "------------------------------------------\n";
	cout << "------------------------------------------\n";
	cout << "------------------------------------------\n";

	obj1.IprintOptimalStrategy();
	obj2.IprintOptimalStrategy();
	system("pause");
	return 0;
}