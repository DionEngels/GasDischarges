// include for input/output
#include <iostream>
//#include <vector>

// set the namespace
using namespace std;

// the main program function
int main() {
  // define variable
  int n1;
  // set value of variable
  n1 = 5;
  // define and initialize variable (constructor initialization)
  int n2(6);
  // define and initialize to sum of variables
  int n3 = n1 + n2;
  // write values to screen
  cout << n1 << "\t" << n2 << "\t" << n3 << endl;

  // now for a class
  /*
          // define variable
          vector<double> v1;
          // make a vector of 4 values
          v1.push_back(5.0);
          v1.push_back(5.0);
          v1.push_back(5.0);
          v1.push_back(5.0);
          // define and initialize variable
          vector<double> v2(4, 6.0);
          // define and initialize to sum of variables
          vector<double> v3(4);
          for (unsigned i = 0; i < v3.size(); i++)
          {
                  v3[i] = v1[i] + v2[i];
          }
          // write values to screen
          for (unsigned i = 0; i < v3.size(); i++)
          {
                  cout << i << "\t" << v1[i] << "\t" << v2[i] << "\t" << v3[i]
     << endl;
          }
  */

  // now the same for self defined a class
  /*
          // define variable
          Field f1(4);
          // set value of variable
          f1 = 5.0;
          // define and initialize variable
          Field f2(4, 6.0);
          // define and initialize to sum of variables
          Field f3 = f1 + f2;
          // write values to screen (using alternative loop)
          unsigned i = 0;
          while(i < f3.size())
          {
                  cout << i << "\t" << f1[i] << "\t" << f2[i] << "\t" << f3[i]
     << endl; i++;
          }

          // alternative initialization, per cell
          i = 0;
          while(i < f3.size())
          {
                  f3[i] = i * i;
                  i++;
          }
          // write to screen using member function
          f3.write(cout);
  */
  // exit program
  return 0;
}
