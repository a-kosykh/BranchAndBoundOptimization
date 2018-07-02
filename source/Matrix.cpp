#include "Matrix.h"

template <typename T>
void Matrix<T>::IenterElements(double varX[][COLUMNS - 1], double freeB[], double objFunc[])
{
	for (unsigned int i = 0; i < _rows; ++i) {
		for (unsigned int j = 0; j < _columns; ++j) {
			if (i < (_rows - 1) && j == 0) {
				_array[i][j] = freeB[i];
			}
			else if (i == _rows - 1 && j > 0) {
				if (!_isMax) {
					_array[i][j] = objFunc[j - 1] * (-1);
				}
				else {
					_array[i][j] = objFunc[j - 1];
				}
			}
			else if (i == _rows - 1 && j == 0) {
				_array[i][j] = 0;
			}
			else {
				_array[i][j] = varX[i][j - 1];
			}
		}
	}

}

template<typename T>
void Matrix<T>::ItoDual()
{
	_toDual();
}

template<typename T>
int Matrix<T>::IfirstStep()
{
	return _stepFirst();
}

template<typename T>
void Matrix<T>::_toDual()
{
	// создание временной транспонированной матрицы
	Matrix<double> temp(COLUMNS, ROWS);
	// замена стремления
	temp._isMax = !(_isMax);

	// транспонирование элементов
	/*for (unsigned int i = 0; i < _rows; ++i) {
		for (unsigned int j = 0; j < _columns; ++j) {
			if (j > 0 && i != _rows - 1) {
				temp._array[j - 1][i + 1] = _array[i][j] * (-1);
			}
			if (j == 0 && i < _rows - 1) {
				temp._array[i][j] = _array[i][j] * (-1);
			}
			if (i == _rows - 1 && j > 0) {
				temp._array[j - 1][0] = _array[i][j] * (-1);
			}
			if (i == _rows - 1 && j == 0) {
				temp._array[temp._rows - 1][0] = 0;
			}
		}
	}*/

	for (unsigned int i = 0; i < temp._rows; ++i) {
		for (unsigned int j = 0; j < temp._columns; ++j) {
			if (i != temp._rows - 1 && j > 0) {
				temp._array[i][j] = _array[j - 1][i + 1] * (-1);
			}
			if (j == 0 && i < temp._rows - 1) {
				temp._array[i][j] = -1;
			}
			if (i == temp._rows - 1 && j > 0) {
				temp._array[i][j] = -1;
			}
			if (i == temp._rows - 1 && j == 0) {
				temp._array[i][j] = 0;
			}
		}
	}
	// возвращение элементов и полей в текущую матрицу
	*this = temp;
	_varBase = temp._varBase;
	_varFree = temp._varFree;
	_isDual = true;
	_isMax = temp._isMax;
}

template <typename T>
int Matrix<T>::_stepFirst()
{
	bool isSorted = false;

	while (!isSorted) { // проверка на то, что все свободные члены больше нуля
		for (unsigned int i = 0; i < _rows - 1; ++i) {
			if (_array[i][0] < 0) {
				// если нет, то запускается функция с жордановым преобразованием (номер строки, номер шага)
				if (_wrongCi(i, 1) == 0) // если нет, то запускается функция с жордановым преобразованием (номер строки, номер шага)
					break;
				else
					return 1;
			}
		}
		// сумма модулей каждого числа должна быть равна сумме всех чисел
		double sum = 0;
		double sumAbs = 0;
		for (unsigned int i = 0; i < _rows - 1; ++i) {
			sum += _array[i][0];
			sumAbs += abs(_array[i][0]);
		}
		if (sumAbs == sum) {
			isSorted = true;
		}
	}

	Iprint();
	_printSolution(0);
	cout << "----------------------------------\n";
	return 0;
}

template <typename T>
int Matrix<T>::_stepSecond()
{
	bool isSorted = false;

	while (!isSorted) { // проверка на то, что все коэфициенты целевой функции меньше нуля
		for (unsigned int j = 1; j < _columns; ++j) {
			if (_array[_rows - 1][j] > 0) {
				// если нет, то запускается функция с жордановым преобразованием (номер разрешающего столбца, номер шага)
				if (_wrongCi(j, 2) == 0) {
					break;
				}
				else {
					return 1;
				}
			}
		}
		// сумма модулей каждого числа должна быть равна сумме всех чисел со знаком минус
		double sum = 0;
		double sumAbs = 0;
		for (unsigned int j = 1; j < _columns; ++j) {
			sum += _array[_rows - 1][j];
			sumAbs += abs(_array[_rows - 1][j]);
		}
		if (sumAbs == -sum) {
			isSorted = true;
		}
	}

	Iprint();
	_printSolution(1);
	cout << "----------------------------------\n";
	return 0;
}

template <typename T>
int Matrix<T>::_wrongCi(unsigned int tempPos, int step)
{
	unsigned int acRow = 0;
	unsigned int acColumn = 0;
	Matrix<T> temp = *this; // cоздание временной матрицы

	if (step == 1) { // При поиске опорного решения определяется разрешающий столбец
		for (unsigned int j = 1; j < _columns; ++j) {
			if (_array[tempPos][j] < 0) {
				acColumn = j;
				break;
			}
			else if (j == _columns - 1 && _array[tempPos][j] >= 0) {
				cout << "NO SOLUTIONS TO THE PROBLEM!\n";
				return 1;
			}
		}
		acRow = tempPos;
	}
	else if (step == 2) {
		// при поиске оптимального решения разрешающий столбец определен на этапе обхода коэф. целевой функции
		acColumn = tempPos;
	}

	cout << "Getting ready to use Jordan Exception\n";
	cout << "Finding ratios of elements in column S to elements in column " << acColumn + 1 << "\n\n";
	// Создание временного массива с отношениями 
	// свободного члена к соответствующему элементу 
	// из разрешающего столбца
	T *tempArray = new T[_rows - 2];
	for (unsigned int i = 0; i < _rows - 1; ++i) {
		if (_array[i][0] == 0 && _array[i][acColumn] < 0) {
			tempArray[i] = numeric_limits<T>::max() - 1;
		}
		else {
			tempArray[i] = _array[i][0] / _array[i][acColumn];
		}
	}

	// поиск минимального соотношения
	T smallest = numeric_limits<T>::max();
	for (unsigned int i = 0; i < _rows - 1; ++i) {
		// если выполняется поиск оптимального решения -- отношение может быть неотрицательным
		if (tempArray[i] < smallest && tempArray[i] >= 0) {
			smallest = tempArray[i];
		}
	}

	// вывод соотношений на экран
	cout << "Ratios:\n";
	for (unsigned int i = 0; i < _rows - 1; ++i) {
		// если выполняется поиск оптимального решения -- отношение может быть только положительным
		if (step == 2 && tempArray[i] < 0) {
			cout << "Row: " << i + 1 << " *** " << tempArray[i] << setprecision(5) << " *** Not acceptable in 2nd step" << "\n";
		}
		else {
			cout << "Row: " << i + 1 << " *** " << tempArray[i] << setprecision(5) << "\n";
		}
	}


	// выбор разрешающего столбца по минимальному отношению
	for (unsigned int i = 0; i < _rows - 1; ++i) {
		if (smallest == tempArray[i]) {
			cout << "Minimal not negative ratio in row " << i + 1 << ": " << smallest << setprecision(5) << "\n\n";
			acRow = i;
			break;
		}
	}

	cout << "Ready to exec Jordan process\n";
	cout << "accepted row: " << acRow + 1 << "\n";
	cout << "accepted column: " << acColumn + 1 << "\n\n";

	// Жорданово исключение по таблице разрешающему столбцу и строке
	_jordanProcess(temp, acRow, acColumn);

	*this = temp;

	if (step == 2) {
		Iprint();
		_printSolution(0);
	}
	if (step == 1) {
		cout << "Trying to find solution\n";
		Iprint();
	}
	cout << "----------------------------------\n\n";
	return 0;
}

template<typename T>
void Matrix<T>::_jordanProcess(Matrix & temp, unsigned int acRow, unsigned int acColumn)
{
	static unsigned int callCount = 0;
	callCount++;
	//cout << "Jordan call counter: " << callCount << "\n";
	// жорданово исключение по матрице, разрешающему строке и столбцу (4 условия)
	for (unsigned int i = 0; i < _rows; ++i) {
		for (unsigned int j = 0; j < _columns; ++j) {
			if (i != acRow && j != acColumn) {
				temp._array[i][j] = _array[i][j] - ((_array[i][acColumn] * _array[acRow][j]) / _array[acRow][acColumn]);
			}
			else if (i == acRow && j != acColumn) {
				temp._array[i][j] = _array[acRow][j] / _array[acRow][acColumn];
			}
			else if (i != acRow && j == acColumn) {
				temp._array[i][j] = (-_array[i][j]) / _array[acRow][acColumn];
			}
			else if (i == acRow && j == acColumn) {
				temp._array[i][j] = 1 / _array[i][j];
			}
		}
	}
	// перестановка переменных местами
	swap(_varBase[acRow], _varFree[acColumn - 1]);
}

template<typename T>
void Matrix<T>::_printCanon()
{
	// выбор x или y в зависимости от задачи (двойственная или нет)
	char varCHAR = 'v';
	if (_isDual) {
		varCHAR = 'u';
	}

	cout << "\nCanon:\n";
	// Вывод целевой функции
	if (_isMax) {
		cout << "-";
	}
	if (_isDual) {
		cout << "W";
	}
	else {
		cout << "Z";
	}
	cout << " = ";
	for (unsigned int j = 1; j < _columns; ++j) {
		if (_isMax) {
			if (_array[_rows - 1][j] * (-1) > 0 && j != 1) {
				cout << " + " << _array[_rows - 1][j] * (-1) << varCHAR << j;
			}
			else {
				cout << _array[_rows - 1][j] * (-1) << varCHAR << j << " ";
			}
		}
		else {
			if (_array[_rows - 1][j] * (-1) > 0 && j != 1) {
				cout << " + " << _array[_rows - 1][j] * (-1) << varCHAR << j;
			}
			else {
				cout << _array[_rows - 1][j] * (-1) << varCHAR << j;
			}
		}
	}
	cout << " --> ";
	if (_isMax) {
		cout << "max";
	}
	else {
		cout << "min";
	}
	cout << "\n\n";
	// Вывод канонического вида
	for (unsigned int i = 0; i < _rows - 1; ++i) {
		for (unsigned int j = 1; j < _columns; ++j) {
			if (_array[i][j] == 0) {
				continue;
			}
			else if (_array[i][j] == 1 && j != 1) {
				cout << "+";
			}
			else if (_array[i][j] == -1) {
				cout << "-";
			}
			else if (_array[i][j] > 0 && j != 1) {
				cout << "+" << _array[i][j];
			}
			else if (_array[i][j] < 0) {
				cout << _array[i][j];
			}
			cout << varCHAR << j;
		}
		cout << "+" << varCHAR << _columns + i;
		cout << " = " << _array[i][0];
		cout << "\n";
	}
	cout << "\n" << varCHAR << "i>=0\n";

	cout << "\n\n";
}

template<typename T>
void Matrix<T>::_printSolution(bool optimal)
{
	char varChar = 'v';
	if (_isDual) {
		varChar = 'u';
	}

	// вывод решения
	cout << "Solutions:\n";
	if (optimal != 1) {
		cout << "NOT ";
	}
	if (!_isMax) {
		cout << "Optimal: F=" << _array[_rows - 1][0] << "\n";
	}
	else {
		cout << "Optimal: F=" << (-1) *_array[_rows - 1][0] << "\n";
	}
	for (unsigned int i = 0; i < _rows - 1; ++i) {
		for (unsigned int j = 0; j < _columns; ++j) {
			if (i == 0 && j > 0) {
				cout << varChar << _varFree[j - 1] << " = 0\n";
			}
			else if (j == 0) {
				cout << varChar << _varBase[i] << " = " << _array[i][j] << "\n";
			}
		}
	}
	cout << "\n";
}

template<typename T>
void Matrix<T>::_checkSolution(double x[][COLUMNS - 1], double freeB[], double objFunc[])
{
	vector<unsigned int> solNum;
	vector<T> solutions;
	for (unsigned int i = 0; i < _rows - 1; ++i) {
		if (_varBase.at(i) <= _columns - 1) {
			solNum.push_back(_varBase.at(i));
			solutions.push_back(_array[i][0]);
		}
	}
	cout << "CHECK:\n";
	Matrix<T> tempMatrix(ROWS, COLUMNS);
	tempMatrix.IenterElements(x, freeB, objFunc);
	if (this->_isDual) {
		tempMatrix.ItoDual();
	}
	for (unsigned int i = 0; i < tempMatrix._rows - 1; ++i) {
		T sum = 0;
		for (unsigned int j = 0; j < solNum.size(); ++j) {
			sum += tempMatrix._array[i][solNum.at(j)] * solutions.at(j);
		}
		if (_isDual) {
			sum *= (-1);
		}
		cout << sum;
		if ((sum >= 1 && _isDual)||(sum<=1 && !_isDual)) {
			cout << " True\n";
		}
	}
}

template<typename T>
Matrix<T>::Matrix()
{
	_rows = 0;
	_columns = 0;
	_array = nullptr;
}

template<typename T>
Matrix<T>::Matrix(unsigned int rows, unsigned int columns)
{
	_rows = rows;
	_columns = columns;

	_array = new T*[_rows];
	for (unsigned int i = 0; i < _rows; ++i) {
		_array[i] = new T[_columns];
	}

	unsigned int tempArgNumber = 0;
	for (unsigned int j = 0; j < _columns - 1; ++j) {
		tempArgNumber++;
		_varFree.push_back(tempArgNumber);
	}
	for (unsigned int i = 0; i < _rows - 1; ++i) {
		tempArgNumber++;
		_varBase.push_back(tempArgNumber);
	}
	_isDual = false;
	_isMax = true;
}

template<typename T>
Matrix<T>::Matrix(const Matrix & get)
{
	_rows = get._rows;
	_columns = get._columns;

	_array = new T*[_rows];
	for (unsigned int i = 0; i < _rows; ++i) {
		_array[i] = new T[_columns];
		for (unsigned int j = 0; j < _columns; ++j) {
			_array[i][j] = get._array[i][j];
		}
	}
}

template<typename T>
Matrix<T>::~Matrix()
{
	for (unsigned int i = 0; i < _rows; ++i) {
		delete[] _array[i];
	}
	delete[] _array;
}

template<typename T>
int Matrix<T>::IsecondStep()
{
	return _stepSecond();
}

template<typename T>
void Matrix<T>::Isimplex()
{	
	_printCanon();
	cout << "----------------------------------\n";
	Iprint();
	cout << "----------------------------------\n";

	if (_stepFirst() == 0) {
		_stepSecond();
	}

	_checkSolution(x, b, c);
}

template<typename T>
void Matrix<T>::Iprint()
{
	char varChar = 'v';
	if (_isDual) {
		varChar = 'u';
	}
	string fChar = "Z";
	if (_isDual) {
		fChar = 'W';
	}
	unsigned int tempRows = _rows + 1;
	unsigned int tempColumns = _columns + 1;

	Matrix<string> temp(tempRows, tempColumns);
	temp.IsetElement(0, 0, "       ");
	temp.IsetElement(0, 1, "S");
	temp.IsetElement(tempRows - 1, 0, fChar);

	for (unsigned int i = 1; i < tempRows - 1; ++i) {
		temp.IsetElement(i, 0, (varChar + to_string(_varBase[i - 1])));
	}

	for (unsigned int j = 2; j < tempColumns; ++j) {
		temp.IsetElement(0, j, (varChar + to_string(_varFree[j - 2])));
	}

	for (unsigned int i = 0; i < tempRows; ++i) {
		for (unsigned int j = 0; j < tempColumns; ++j) {
			if (i != 0 && j != 0) {
				temp.IsetElement(i, j, (to_string(_array[i - 1][j - 1])));
			}
		}
	}
	cout << "Simplex:\n";
	for (unsigned int i = 0; i < tempRows; ++i) {
		for (unsigned int j = 0; j < tempColumns; ++j) {
			string tempStr = temp.IgetElement(i, j);
			if (tempStr.at(0) != fChar.at(0) && tempStr[0] != varChar && tempStr[0] != 'S' && tempStr[0] != ' ') {
				size_t findDot = tempStr.find('.', 0);
				tempStr.erase(findDot + 4);
			}
			std::cout << tempStr << setw(10);
		}
		std::cout << std::endl;
	}
	cout << "\n";
}

template<typename T>
void Matrix<T>::IprintSol(bool error)
{
	if (error) {
		cout << "NO SOLUTIONS!!\n\n";
	}
	else {
		_printSolution(true);
	}
}

template<typename T>
void Matrix<T>::IprintCanon()
{
	_printCanon();
}

template<typename T>
void Matrix<T>::IcheckSolutions()
{
	_checkSolution(x, b, c);
}

template<typename T>
void Matrix<T>::IprintOptimalStrategy()
{
	char varChar = 'y';
	char fChar = 'h';
	if (_isDual) {
		varChar = 'x';
		fChar = 'g';
	}
	cout << "\n";
	T backFunc = 1 / _array[_rows - 1][0];
	if (_isMax) {
		backFunc *= -1;
	}
	cout << fChar << "(" << varChar << ")" << " = " << backFunc << "\n";
	T sum = 0;
	for (unsigned int i = 0; i < _rows - 1; ++i) {
		if (_varBase.at(i) < _columns) {
			cout << varChar << _varBase.at(i) << " = " << _array[i][0] * backFunc << "\n";
			sum += (_array[i][0] * backFunc);
		}
	}
	cout << "Other " << varChar << "'s = 0\n";
	cout << "Check sum = " << sum << "\n";
}

template<typename T>
void Matrix<T>::setDual(bool set)
{
	_isDual = set;
}

template<typename T>
bool Matrix<T>::getDual()
{
	return this->_isDual;
}

template<typename T>
vector<unsigned int> Matrix<T>::getVarBase()
{
	return this->_varBase;
}

template<typename T>
vector<unsigned int> Matrix<T>::getVarFree()
{
	return this->_varFree;
}

template<typename T>
void Matrix<T>::setBaseFree(vector<unsigned int> vb, vector<unsigned int> vf)
{
	this->_varBase = vb;
	this->_varFree = vf;
}

template<typename T>
void Matrix<T>::IsetElementsManual()
{
	for (unsigned int i = 0; i < _rows; ++i) {
		for (unsigned int j = 0; j < _columns; ++j) {
			if (j == 0 && i != _rows - 1) {
				double key = 0;
				cout << "s" << i << ": ";
				cin >> key;
				_array[i][j] = key;
			}
			else if (i == _rows - 1) {
				double key = 0;
				cout << "F[" << j << "]: ";
				cin >> key;
				if (!_isMax) {
					_array[i][j] = key * (-1);
				}
				else {
					_array[i][j] = key;
				}
			}
			else {
				double key = 0;
				cout << "x[" << i << "]" << "[" << j << "]: ";
				cin >> key;
				_array[i][j] = key;
			}
		}
	}
}

template<typename T>
T Matrix<T>::IgetElement(unsigned int row, unsigned int column)
{
	return _array[row][column];
}

template<typename T>
void Matrix<T>::IsetElement(unsigned int row, unsigned int column, T key)
{
	_array[row][column] = key;
}

template<typename T>
T ** Matrix<T>::getArray()
{
	return _array;
}

template<typename T>
unsigned int Matrix<T>::getRows()
{
	return this->_rows;
}

template<typename T>
unsigned int Matrix<T>::getColumns()
{
	return this->_columns;
}

template<typename T>
Matrix<T>& Matrix<T>::operator=(const Matrix & right)
{
	if (this == &right) {
		return *this;
	}
	for (unsigned int i = 0; i < _rows; ++i) {
		delete[]_array[i];
	}
	delete[]_array;
	_array = new T*[right._rows];
	for (unsigned int i = 0; i < right._rows; ++i) {
		_array[i] = new T[right._columns];
		for (unsigned int j = 0; j < right._columns; ++j)
			_array[i][j] = right._array[i][j];
	}
	_rows = right._rows;
	_columns = right._columns;
	//_varBase = right._varBase;
	//_varFree = right._varFree;
	return *this;
}

template<typename T>
bool Matrix<T>::operator==(const Matrix & right) const
{
	if ((_rows == right._rows) && (_columns == right._columns)) {
		for (unsigned int i = 0; i < _rows; ++i)
			for (unsigned int j = 0; j < _columns; ++j) {
				if (_array[i][j] != right._array[i][j]) { return false; }
			}
	}
	else return false;
}

template class Matrix<double>;