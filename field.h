/** \file
 *
 *  This file contains the declaration of the Field class.
 */

#ifndef H_FIELD_H
#define H_FIELD_H

#include "basics.h"
#include <cassert>
#include <iostream>
#include <vector>


/** A representation of (one-dimensional) real-valued field variables.
 *  Fields can be constructed by specifying the number of field points
 *  and the initial value. The number of grid points can be obtained by
 *  calling the member size(). Values of a field f can be read and written
 *  using the syntax f[i], where i is in the range [0,size()). Furthermore,
 *  members assigning a value to all field elements and for writing fields
 *  to streams have been provided.
 *
 *  \author Jan van Dijk
 */
class Field {
public:
  /** Constructs a field of \a size points with initial value \a value.
   *  The argument \a value may be omitted, in which case the default
   *  value of 0 is assumed.
   */
  Field(unsigned size, double value = 0.0) : m_data(size, value) {}
  /** This desctructor is called when a field variable goes out of scope.
   *  It does nothing by itself. However, a virtual _destructor_ is required
   *  when we derive classes with polymorphic functions from this class.
   */
  virtual ~Field() {}
  /** This accessor provides access to element \a i of a field variable f
   *  with the syntax 'f[i]'. This version returns a (constant) copy of that
   *  element. Since that does not modify the field, this has been made a
   *  constant member.
   */
  const double operator[](unsigned i) const {
    assert(i < m_data.size());
    return m_data[i];
  }
  /** This accessor provides access to element \a i of a field variable f
   *  with the syntax 'f[i]'. This version returns a reference to that
   *  element. Since that can be used to modify the field (for example,
   *  by writing f[i]=1) this possibly modifies the field. Therefore this
   *  member is _not_ constant.
   */
  double &operator[](unsigned i) {
    assert(i < m_data.size());
    return m_data[i];
  }
  /** this member returns the number of field points.
   *  Valid field indices are in the range [0,size()).
   */
  unsigned size() const { return m_data.size(); }
  /** Assigns the value \a val to all points of the field.
   *  As an example, to set all elements of a Field-variable
   *  f to 1.0, you may just write 'f=1.0'.
   */
  const Field &operator=(double val) {
    for (unsigned i = 0; i < m_data.size(); i++) {
      m_data[i] = val;
    }
    return *this;
  }
  /** Addition operator
   */
  const Field operator+(Field &f_in) {
    if (f_in.size() != m_data.size()) {
      npfFatalError("Fields not equal in size");
    }
    Field f_temp(m_data.size());
    for (unsigned i = 0; i < m_data.size(); i++) {
      f_temp[i] = m_data[i] + f_in[i];
    }
    return f_temp;
  }
  /** This writes the field to the output stream \a os.
   *  The first column of the output contains the field point indices,
   *  the second the field values.
   *  As an example, to write the Field-variable f to the
   *  standard output stream, you may write 'f.write(std::cout);'.
   */
  void write(std::ostream &os) const {
    for (unsigned i = 0; i < m_data.size(); i++) {
      os << i << '\t' << m_data[i] << '\n';
    }
  }
  /** This writes the field to stdout with some extra commands so the
   *  field data can be plotted by piping to gnuplot.
   */
  void plot(std::string xlabel = "", std::string ylabel = "",
            std::string title = "") const {
    std::cout << "unset key" << std::endl;
    std::cout << "set xlabel \"" << xlabel << "\"" << std::endl;
    std::cout << "set ylabel \"" << ylabel << "\"" << std::endl;
    std::cout << "set title \"" << title << "\"" << std::endl;
    std::cout << "plot '-' w l" << std::endl;
    write(std::cout);
    std::cout << "e" << std::endl;
  }

private:
  /** The internal representation of this class uses the standard
   *  vector template, std::vector. This is an implementation detail,
   *  users of the Field class have access to the field data only
   *  through the limited set of accessors in the public section.
   */
  std::vector<double> m_data;
};

#endif // H_FIELD_H
