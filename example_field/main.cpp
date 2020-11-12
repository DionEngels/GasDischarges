#include "field.h"
#include <iostream>

int main()
{
	Field field(10);
	for (unsigned i=0; i<field.size(); ++i)
	{
		field[i] = i*i;
	}
	field.write(std::cout);
	return 0;
}
