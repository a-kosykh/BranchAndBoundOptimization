#pragma once
#include "Elements.h"
#include <iostream>
#include <math.h>
#include <limits>
#include <vector>
#include <string>
#include <algorithm>
#include <iomanip>

using namespace std;

template <typename T>
class Matrix {
private:
	bool _isMax; // стремится ли целевая функция к максимуму
	bool _isDual; // двойственная ли задача

	unsigned int _rows; // количество строк
	unsigned int _columns; // количество столбцов

	T**	_array; // двумерный массив

	vector<unsigned int> _varFree; // вектор номеров свободных переменных
	vector<unsigned int> _varBase; // вектор номеров базисных переменных

								   // ввод элементов
	void _toDual(); // приведение задачи к двойственной
	int _stepFirst(); // поиск опорного решения
	int _stepSecond(); // поиск оптимального решения
	int _wrongCi(unsigned int tempPos, int step); // функция, вызываемая при 
	void _jordanProcess(Matrix &temp, unsigned int acRow, unsigned int acColumn); // жорданово исключение

	void _printCanon(); // вывод канонического вида
	void _printSolution(bool optimal); // вывод решений
	void _checkSolution(double x[][COLUMNS - 1], double freeB[], double objFunc[]); // проверка решения

public:

	//конструкторы
	Matrix();
	Matrix(unsigned int rows, unsigned int columns);
	Matrix(const Matrix& get);
	~Matrix();

	void IenterElements(double varX[][COLUMNS - 1], double freeB[], double objFunc[]); // заполнение элементов
	void ItoDual();
	int IfirstStep(); // первый шаг
	int IsecondStep(); // второй шаг
	void Isimplex(); // вызов функций симплекс метода
	void Iprint(); // вывод
	void IprintSol(bool error); // вывод решений
	void IprintCanon(); // вывод канонического вида
	void IcheckSolutions(); // проверка решений
	void IprintOptimalStrategy();

	void setDual(bool set); // установка двойственной задачи
	bool getDual();

	vector<unsigned int> getVarBase(); // геттер индексов базисных переменных
	vector<unsigned int> getVarFree(); // геттер индексов свободных переменных
	void setBaseFree(vector<unsigned int> vb, vector<unsigned int> vf); // установка индексов

	void IsetElementsManual(); // ручной ввод

	T IgetElement(unsigned int row, unsigned int column); // получение определения элемента из таблицы
	void IsetElement(unsigned int row, unsigned int column, T key); // установка элемента

	T** getArray(); // геттер массива

	unsigned int getRows(); // геттер количества строк
	unsigned int getColumns(); // кол-ва столбцов

	Matrix<T>& operator= (const Matrix& right);
	bool operator==(const Matrix& right) const;
};