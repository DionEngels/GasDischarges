#include "field.h"
#include <iostream>

int main() {
  Field field(10);
  for (unsigned i = 0; i < field.size(); ++i) {
    field[i] = i * i;
  }
  field.write(std::cerr);
  field.plot("position", "value", "example field");
  return 0;
}
