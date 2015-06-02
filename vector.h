#pragma once
#include "matrix.h"
class vector :
	public matrix
{
public:
	vector(int dim);
	float norm();
	~vector(void);
};

