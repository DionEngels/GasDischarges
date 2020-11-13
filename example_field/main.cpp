#include "field.h"
#include <iostream>

int main() {
  Field field(10);
  field = 1.0 / 3;
  field.write(std::cerr);
  field.plot("position", "value", "example field");
  return 0;
}
