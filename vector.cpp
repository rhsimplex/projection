#include "vector.h"
#include <math.h>


vector::vector(int dim):matrix(dim, 1)
{
}

float vector::norm(void)
{
	matrix result = (this->tr())*(*this);
	return sqrt(result.entries[0][0]);
}


vector::~vector(void)
{
}
