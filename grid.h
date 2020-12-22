/** \file
 *
 *  Declarations of grid-related classes for the numerical plasma physics
 *  course.
 *
 *  \author Bart Hartgers and Jan van Dijk
 */
#ifndef H_GRID_H
#define H_GRID_H

#include "field.h"
#include "physconst.h"
#include <vector>

/** A base class for the grid.
 *  It defines various properties on nodal and cell boundary (ew) points.
 *  All coordinates, surfaces and volumes are expressed in m, m^2 and m^3,
 *  respectively.
 *
 *  The class has a protected constructor which requires the number of
 *  grid cells, the lower and upper boundaries of the grid coordinate,
 *  and a pointer to a free function which returns the perpendicular area
 *  as a function of the grid coordinate. That function (alone) determines
 *  the nature of the grid coordinate (Cartesian, Cylindrical, Spherical)
 *
 *  In the center of each cell, a nodal point
 *  is located. In addition, two nodal points are present on the left
 *  boundary of the first, and on the right boundary on the last cell.
 *  All cells have the same size. Note that as a result of this layout,
 *  the distances between the first two and between the last two points
 *  is only half the cell size.
 *
 *  \author Bart Hartgers and Jan van Dijk
 */
class Grid {
public:
  /// number of nodal points
  unsigned num_np() const { return m_pos_np.size(); }
  /// physical coordinates of all nodal points [m]
  const Field &pos_np() const { return m_pos_np; }
  /// physical coordinate of the nodal point \a i [m]
  double pos_np(unsigned i) const { return m_pos_np[i]; }
  /// length of control volume \a i [m]
  double dx_np(unsigned i) const { return m_dx[i]; }
  /// perpendicular area at nodal point \a i [m^2]
  double area_np(unsigned i) const { return m_area_np[i]; }
  /// volume of control volume \a i [m^3]
  double vol_np(unsigned i) const { return m_vol_np[i]; }

  /// number of cell boundary points
  unsigned num_ew() const { return m_pos_ew.size(); }
  /// physical coordinates of all cell boundary points [m]
  const Field &pos_ew() const { return m_pos_ew; }
  /// physical coordinate of the cell boundary point \a i [m]
  double pos_ew(unsigned i) const { return m_pos_ew[i]; }
  /// perpendicular area at cell boundary point \a i [m^2]
  double area_ew(unsigned i) const { return m_area_ew[i]; }

  /// mesh size [m] (size of an internal control volume)
  double del() const { return m_del; }
  /// total grid volume [m^3]
  double volume() const { return m_vol; }

  /** Writes the field variable \a f to the stream \a os.
   *  The field must be a nodal or ew field defined on this grid,
   *  the function uses the size of the field to determine what
   *  kind of field it is. The output consists of two columns
   *  that represent the coordinate-values and the field values,
   *  respectively.
   */
  void write(std::ostream &os, const Field &f) const;
  void write(FILE *pipe, const Field &f) const;
  /** Writes the field variables \a f1 and \a f2 to the stream \a os.
   *  Both fields must be defined on the same points, either
   *  the nodal or ew points defined for this grid.
   *  The function uses the size of the field to determine what
   *  kind of fields were provided. The output consists of three columns
   *  that represent the coordinate-values and the values of the
   *  two fields, respectively.
   */

  void write(std::ostream &os, const Field &f, const Field &f2) const;
  void write(FILE *pipe, const Field &f, const Field &f2) const;
  /** Writes the field variable \a f to stdout for piping to gnuplot.
   *  Optional axis labels \a xlabel and \a ylabel can be provided.
   *  The field must be defined on either the nodal or the cell boundary
   *  points of this grid.
   */
  void plot(const Field &f, std::string xlabel = "", std::string ylabel = "",
            std::string title = "", double ymin = 0, double ymax = -1);
  void plot_freeze(const Field &f, std::string xlabel = "",
                   std::string ylabel = "", std::string title = "",
                   double ymin = 0, double ymax = -1);
  /** Writes the field variables \a f1 and \a f2 to stdout for piping to
   *  gnuplot. Optional axis labels \a xlabel and \a ylabel can be provided.
   *  Both fields must be defined on the same mesh, either the nodal or the
   *  cell boundary points of this grid.
   */
  void plot(const Field &f1, const Field &f2, std::string xlabel = "",
            std::string ylabel = "", std::string y2label = "",
            std::string title = "") const;
  /// enum values for the supported grids
  enum GridType { Cartesian, Cylindrical, Spherical };
  /** This constructor creates a control volume grid consisting
   *  of \a n grid cells (\a n+2 nodal points). The coordinate
   *  is in the range [\a cmin, \a cmax]. The function \a area_func
   *  provides the area of the grid perpendicular to the coordinate
   *  direction as a function of the grid coordinate.
   */
  Grid(unsigned n, double cmin, double cmax, GridType grid_type);

  FILE *pipe = popen("gnuplot", "w");
  bool first_plot = true;

private:
  const GridType
      m_grid_type; ///< type of the grid (Cartesian, Cylindrical, Spherical)
  double m_cmin;   ///< left physical coordinate
  double m_cmax;   ///< right physical coordinate
  double m_del;    ///< the mesh size (equal to the length of interior cells)

  Field m_pos_np; ///< physical positions of the nodal points [m]
  Field m_pos_ew; ///< physical positions of the cell boundary points [m]
  Field m_dx;     ///< length of the control volumes [m]

  Field m_area_np; ///< perpendicular area at the nodal points [m^2]
  Field m_area_ew; ///< perpendicular area at the cell boundary points [m^2]

  Field m_vol_np; ///< volumes of the grid cells [m^3]

  double m_vol; ///< total volume of the grid

  /// returns the grid type-dependent perpendicular area [m^2] at coordinate
  /// value \a c.
  double area_func(double c) const;
};

#endif // H_GRID_H
