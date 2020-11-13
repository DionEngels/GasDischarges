#include "grid.h"
#include "basics.h"

Grid::Grid(unsigned n, double cmin, double cmax, GridType grid_type)
    : m_grid_type(grid_type), m_cmin(cmin), m_cmax(cmax),
      m_del((cmax - cmin) / n), m_pos_np(n + 2), m_pos_ew(n + 1), m_dx(n + 2),
      m_area_np(n + 2), m_area_ew(n + 1), m_vol_np(n + 2) {
  if (m_cmin >= m_cmax) {
    npfFatalError("Grid cmin>=cmax");
  }
  // in the code below we use the area function which was passed to the
  // constructor to calculate the perpendicular surface. That (and only
  // that) determines the nature of the grid (Cartesian, cylindrical or
  // spherical).

  // the nodal points
  m_vol = 0;
  for (unsigned i = 0; i < m_pos_np.size(); ++i) {
    if (i == 0) {
      m_pos_np[i] = m_cmin;
      m_dx[i] = 0;
    } else if (i == m_pos_np.size() - 1) {
      m_pos_np[i] = m_cmax;
      m_dx[i] = 0;
    } else {
      m_pos_np[i] = m_cmin + (i - 0.5) * m_del;
      m_dx[i] = m_del;
    }
    m_area_np[i] = area_func(m_pos_np[i]);
    m_vol_np[i] = m_dx[i] * m_area_np[i];
    // add the contribution of this cell to the grid volume
    m_vol += m_vol_np[i];
  }

  // the control volume boundary points
  for (unsigned i = 0; i < m_pos_ew.size(); ++i) {
    m_pos_ew[i] = m_cmin + m_del * i;
    m_area_ew[i] = area_func(m_pos_ew[i]);
  }
}

void Grid::write(std::ostream &os, const Field &f) const {
  if (f.size() == m_pos_np.size()) {
    for (unsigned i = 0; i < f.size(); i++) {
      os << pos_np(i) << '\t' << f[i] << '\n';
    }
  } else if (f.size() == m_pos_ew.size()) {
    for (unsigned i = 0; i < f.size(); i++) {
      os << pos_ew(i) << '\t' << f[i] << '\n';
    }
  } else {
    npfFatalError("Grid::write: field is not defined on this grid.");
  }
}

void Grid::write(FILE *pipe, const Field &f) const {
  if (f.size() == m_pos_np.size()) {
    for (unsigned i = 0; i < f.size(); i++) {
      fprintf(pipe, "%f\t%f\n", pos_np(i), f[i]);
    }
  } else if (f.size() == m_pos_ew.size()) {
    for (unsigned i = 0; i < f.size(); i++) {
      fprintf(pipe, "%f\t%f\n", pos_ew(i), f[i]);
    }
  } else {
    npfFatalError("Grid::write: field is not defined on this grid.");
  }
}

double Grid::area_func(double c) const {
  switch (m_grid_type) {
  case Cartesian:
    return 1.0;
    break;
  case Cylindrical:
    return 2 * PhysConst::pi * c;
    break;
  case Spherical:
    return 4 * PhysConst::pi * c * c;
    break;
  }
  npfFatalError(
      "Grid: illegal grid type. Must be Cartesian, Cylindrical or Spherical.");
  abort();
}

void Grid::write(std::ostream &os, const Field &f1, const Field &f2) const {
  if (f1.size() != f2.size()) {
    npfFatalError("Grid::write: fields must be defined on the same mesh.");
  }
  if (f1.size() == m_pos_np.size()) {
    for (unsigned i = 0; i < f1.size(); i++) {
      os << pos_np(i) << '\t' << f1[i] << '\t' << f2[i] << '\n';
    }
  } else if (f1.size() == m_pos_ew.size()) {
    for (unsigned i = 0; i < f1.size(); i++) {
      os << pos_ew(i) << '\t' << f1[i] << '\t' << f2[i] << '\n';
    }
  } else {
    npfFatalError("Grid::write: field is not defined on this grid.");
  }
}

void Grid::write(FILE *pipe, const Field &f1, const Field &f2) const {
  if (f1.size() != f2.size()) {
    npfFatalError("Grid::write: fields must be defined on the same mesh.");
  }
  if (f1.size() == m_pos_np.size()) {
    for (unsigned i = 0; i < f1.size(); i++) {
      fprintf(pipe, "%f\t%f\t%f\n", pos_np(i), f1[i], f2[i]);
    }
  } else if (f1.size() == m_pos_ew.size()) {
    for (unsigned i = 0; i < f1.size(); i++) {
      fprintf(pipe, "%f\t%f\t%f\n", pos_ew(i), f1[i], f2[i]);
    }
  } else {
    npfFatalError("Grid::write: field is not defined on this grid.");
  }
}

void Grid::plot(const Field &f, std::string xlabel, std::string ylabel,
                std::string title) const {

  FILE *pipe = popen("gnuplot", "w");
  if (pipe != NULL) {

    fprintf(pipe, "unset key\n");
    fprintf(pipe, "set xlabel \"%s\"\n", xlabel.c_str());
    fprintf(pipe, "set ylabel \"%s\"\n", ylabel.c_str());
    fprintf(pipe, "set xrange [%f:%f]\n", pos_np(0), pos_np(num_np() - 1));
    fprintf(pipe, "set title \"%s\"\n", title.c_str());
    fprintf(pipe, "plot '-' w l\n");
    write(pipe, f);
    fprintf(pipe, "%s\n", "e");
    fflush(pipe);

    char output[5];

    // wait for command line input
    std::cin.clear();
    std::cin.ignore(std::cin.rdbuf()->in_avail());
    std::cin.get(output, 5);

    if (strcmp(output, "save") == 0) {
      fprintf(pipe, "set terminal pdf\n");
      fprintf(pipe, "set output \"figs/%s.pdf\"\n", title.c_str());
      fprintf(pipe, "replot\n");
      fprintf(pipe, "unset output\n");
      fprintf(pipe, "unset terminal\n");
    }

    pclose(pipe);
  } else {
    std::cerr << "Creating plot failed!" << std::endl;
  }
}

void Grid::plot(const Field &f1, const Field &f2, std::string xlabel,
                std::string ylabel, std::string y2label,
                std::string title) const {
  FILE *pipe = popen("gnuplot", "w");
  if (pipe != NULL) {

    fprintf(pipe, "unset key\n");
    fprintf(pipe, "set ytics nomirror\n");
    fprintf(pipe, "set y2tics nomirror\n");
    fprintf(pipe, "set xlabel \"%s\"\n", xlabel.c_str());
    fprintf(pipe, "set ylabel \"%s\"\n", ylabel.c_str());
    fprintf(pipe, "set y2label \"%s\"\n", y2label.c_str());
    fprintf(pipe, "set xrange [%f:%f]\n", pos_np(0), pos_np(num_np() - 1));
    fprintf(pipe, "set title \"%s\"\n", title.c_str());
    fprintf(pipe, "plot '-' u 1:2 axis x1y1 w l, '' u 1:2 axis x1y2 w l\n");
    write(pipe, f1);
    fprintf(pipe, "%s\n", "e");
    write(pipe, f2);
    fprintf(pipe, "%s\n", "e");
    fflush(pipe);

    char output[5];

    // wait for command line input
    std::cin.clear();
    std::cin.ignore(std::cin.rdbuf()->in_avail());
    std::cin.get(output, 5);

    if (strcmp(output, "save") == 0) {
      fprintf(pipe, "set terminal pdf\n");
      fprintf(pipe, "set output \"figs/%s.pdf\"\n", title.c_str());
      fprintf(pipe, "replot\n");
      fprintf(pipe, "unset output\n");
      fprintf(pipe, "unset terminal\n");
    }

    pclose(pipe);
  } else {
    std::cerr << "Creating plot failed!" << std::endl;
  }
}
