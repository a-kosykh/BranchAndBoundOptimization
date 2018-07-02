#pragma once
static const int ROWS = 6;
static const int COLUMNS = 5;

static double x[ROWS - 1][COLUMNS - 1] = {
	{ 1, 3, 4, 16 },
	{ 11, 3, 7, 18 },
	{ 6, 15, 16, 7 },
	{ 15, 10, 0, 7 },
	{ 10, 7, 10, 5}
};


static double  c[COLUMNS - 1] = { 1,1,1,1 };
static double  b[ROWS - 1] = { 1,1,1,1,1 };