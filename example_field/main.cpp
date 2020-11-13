#include "field.h"
#include "grid.h"
#include <iostream>

int main() {

  Grid cart_grid(10, 0.0e-2, 1.0e-2, Grid::Cartesian);

  Field field(11);
  field = 1.0 / 3;

  cart_grid.write(std::cout, field);
  cart_grid.plot(field, "position", "value", "example field");
  return 0;
}
