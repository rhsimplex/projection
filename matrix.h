#pragma once
class matrix
{
public:
	matrix(float** data, int numrows, int numcolumns);
	matrix(int numrows, int numcolumns);
	matrix();
    int getcolumns();
	int getrows();
	matrix getrow(int r);
	matrix getcolumn(int c);
	static matrix orthonormalized(const matrix &b);
	void virtual print();
	void setentry(int row, int col, float val);
	float getentry(int row, int col);
	virtual ~matrix(void);
	float norm(void);
	matrix tr();
	void elementmultiply(float coeff);
	matrix operator +(const matrix &toadd);
	matrix operator -(const matrix &tosub);
	matrix operator *(const matrix &tomult);
	bool operator ==(const matrix &tocheck);
	void operator =(const matrix &tocopy);

protected:
	float** entries;	
	int rows;
	int columns;
private:

};

