#include "matrix.h"
#include <iostream>
#include <assert.h>
#include <math.h>

using namespace std;

matrix::matrix(float** data, int numrows, int numcolumns)
{
	entries=data;
	rows=numrows;
	columns=numcolumns;
}
matrix::matrix()
{
	rows = 0;
	columns = 0;
	entries = new float*[rows];
}
matrix::matrix(int numrows, int numcolumns)
{
	rows=numrows;
	columns=numcolumns;
	entries = new float*[rows];
	for(int i = 0; i < rows; i++)
	{
		entries[i] = new float[columns];
		for(int j = 0; j < columns; j++)
		{
			entries[i][j]=0.0;
		}
	}
}

void matrix::setentry(int row, int col, float val)
{
	entries[row][col]=val;
}

float matrix::getentry(int row, int col)
{
	float val = entries[row][col];
	return val;
}
void matrix::print()
{
	for(int i = 0; i < rows; i++)
	{
		for(int j = 0; j < columns; j++)
			cout << entries[i][j] << "\t";
		cout << "\n";
	}
}

int matrix::getcolumns()
{
	return columns;
}

int matrix::getrows()
{
	return rows;
}

matrix matrix::operator +(const matrix &toadd)
{
	assert(rows == toadd.rows);
	assert(columns == toadd.columns);
	matrix result(rows, columns);
	for(int i = 0; i < rows; i++)
	{
		for(int j = 0; j < columns; j++)
		{
			result.setentry(i,j, entries[i][j]+toadd.entries[i][j] );
		}
	}
	return result;
}
matrix matrix::operator -(const matrix &tosub)
{
	assert(rows == tosub.rows);
	assert(columns == tosub.columns);
	matrix result(rows, columns);
	for(int i = 0; i < rows; i++)
	{
		for(int j = 0; j < columns; j++)
		{
			result.setentry(i,j, entries[i][j]-tosub.entries[i][j] );
		}
	}
	return result;
}
matrix matrix::operator *(const matrix &tomult)
{
	assert(columns == tomult.rows);
	matrix result(rows, tomult.columns);
	for(int i = 0; i < rows; i++)
	{
		for(int j = 0; j < tomult.columns; j++)
		{
			for(int k = 0; k < columns; k++)
				result.entries[i][j]+=entries[i][k]*tomult.entries[k][j];
		}
	}
	return result;
}

//Only works for vectors! (1 x n matrices)
float matrix::norm()
{
	assert(columns==1 || rows==1);
	matrix result = this->tr() * *this;
	return sqrt(result.entries[0][0]);
}

matrix matrix::tr(void)
{
	matrix result(columns, rows);
	for(int i = 0; i < rows; i++)
	{
		for(int j = 0; j < columns; j++)
		{
			result.entries[j][i] = entries[i][j];
		}
	}
	return result;
}

bool matrix::operator ==(const matrix &tocheck)
{
	if(columns != tocheck.columns || rows != tocheck.rows )
		return false;
	for(int i = 0; i < rows; i++)
	{
		for(int j = 0; j < columns; j++)
		{
			if(entries[i][j] != *tocheck.entries[i,j])
				return false;
		}
	}
	return true;
}

void matrix::operator =(const matrix &tocopy)
{

	for(int i = 0; i < rows; i++)
	{
		delete[] entries[i];
	}
	delete[] entries;
	rows = tocopy.rows;
	columns = tocopy.columns;
	entries = new float*[rows];
	for(int i = 0; i < rows; i++)
	{
		entries[i] = new float[columns];
		for(int j = 0; j < columns; j++)
		{
			entries[i][j]=tocopy.entries[i][j];
		}
	}
	
}

matrix::~matrix(void)
{
	for(int i = 0; i < rows; i++)
	{
		delete[] entries[i];
	}
	delete[] entries;
}
matrix matrix::getrow(int r)
{
	assert(r < rows);
	matrix row(1, columns);
	for(int i = 0; i < columns; i++)
	{
		row.entries[0][i] = entries[r][i];
	}
	return row;
}

matrix matrix::getcolumn(int c)
{
	assert(c < columns);
	matrix column(rows, 1);
	for(int i = 0; i < rows; i++)
	{
		column.entries[i][0] = entries[i][c];
	}
	column.print();
	return column;
}
void matrix::elementmultiply(float coeff)
{
	for(int i = 0; i < rows; i++)
	{
		for(int j = 0; j < columns; j++)
			entries[i][j]=coeff*entries[i][j];
	}
}
matrix matrix::orthonormalized(const matrix &b)
{
	float dotprodij;
	float dotprodii;
	float dotprodjj;
	matrix a = b;
	cout << a.columns << "\n";
	matrix v;
	matrix u;
	a.print();
	v=a.getcolumn(0);
	v.print();
	for(int j = 0; j < a.columns; j++)
	{
		for(int i = 0; i < (j-1) ; i++)
		{

			u = a.getcolumn(i);
			v = a.getcolumn(j);
			u.print();
			v.print();

			dotprodij = ((a.getcolumn(j)).tr() * (a.getcolumn(i))).entries[0][0];
			dotprodii = ((a.getcolumn(i)).tr() * (a.getcolumn(i))).entries[0][0];

			for(int k = 0; k < a.rows; k++)
			{
				a.entries[j][k] =  a.entries[j][k] - dotprodij/dotprodii*a.entries[j][k];
			}
		}
		dotprodjj = (a.getcolumn(j).tr() * a.getcolumn(j)).entries[0][0];
		for(int k = 0; k < b.rows; k++)
		{
			a.entries[j][k] = a.entries[j][k]/(a.getcolumn(j).norm());
		}
	}
	cout << "\n";
	a.print();
	return a;

}