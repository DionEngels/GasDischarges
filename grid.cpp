#include "grid.h"
#include "basics.h"

Grid::Grid( unsigned n, double cmin, double cmax, GridType grid_type)
 : m_grid_type(grid_type),
   m_cmin(cmin),
   m_cmax(cmax),
   m_del((cmax-cmin)/n),
   m_pos_np( n+2 ),
   m_pos_ew( n+1 ),
   m_dx( n+2 ),
   m_area_np( n+2 ),
   m_area_ew( n+1 ),
   m_vol_np( n+2 )
{
	if (m_cmin>=m_cmax)
	{
		npfFatalError("Grid cmin>=cmax");	
	}
	// in the code below we use the area function which was passed to the
	// constructor to calculate the perpendicular surface. That (and only
	// that) determines the nature of the grid (Cartesian, cylindrical or
	// spherical).

	// the nodal points
	m_vol=0;
	for( unsigned i=0; i<m_pos_np.size(); ++i)
	{
		if (i==0)
		{
			m_pos_np[i] = m_cmin;
			m_dx[i]     = 0;
		}
		else if (i==m_pos_np.size()-1)
		{
			m_pos_np[i] = m_cmax;
			m_dx[i]     = 0;
		}
		else
		{
			m_pos_np[i] = m_cmin + (i-0.5)*m_del;
			m_dx[i]     = m_del;
		}
		m_area_np[i] = area_func(m_pos_np[i]);
		m_vol_np[i]  = m_dx[i]*m_area_np[i];
		// add the contribution of this cell to the grid volume
		m_vol +=m_vol_np[i];
	}

	// the control volume boundary points
	for( unsigned i=0; i<m_pos_ew.size(); ++i)
	{
		m_pos_ew[i]  = m_cmin + m_del*i;
		m_area_ew[i] = area_func(m_pos_ew[i]);
	}	
}

void Grid::write(std::ostream& os, const Field& f) const
{
	if (f.size()==m_pos_np.size())
	{
		for (unsigned i=0; i<f.size(); i++)
		{
			os << pos_np(i) << '\t' << f[i] << '\n';
		}
	}
	else if (f.size()==m_pos_ew.size())
	{
		for (unsigned i=0; i<f.size(); i++)
		{
			os << pos_ew(i) << '\t' << f[i] << '\n';
		}
	}
	else
	{
		npfFatalError("Grid::write: field is not defined on this grid.");
	}
}

double Grid::area_func(double c) const
{
	switch (m_grid_type)
	{
		case Cartesian:   return 1.0;                 break;
		case Cylindrical: return 2*PhysConst::pi*c;   break;
		case Spherical:   return 4*PhysConst::pi*c*c; break;
	}
	npfFatalError("Grid: illegal grid type. Must be Cartesian, Cylindrical or Spherical.");
	abort();
}

void Grid::write(std::ostream& os, const Field& f1, const Field& f2) const
{
	if (f1.size()!=f2.size())
	{
		npfFatalError("Grid::write: fields must be defined on the same mesh.");
	}
	if (f1.size()==m_pos_np.size())
	{
		for (unsigned i=0; i<f1.size(); i++)
		{
			os << pos_np(i) << '\t' << f1[i] << '\t' << f2[i] << '\n';
		}
	}
	else if (f1.size()==m_pos_ew.size())
	{
		for (unsigned i=0; i<f1.size(); i++)
		{
			os << pos_ew(i) << '\t' << f1[i] << '\t' << f2[i] << '\n';
		}
	}
	else
	{
		npfFatalError("Grid::write: field is not defined on this grid.");
	}
}

void Grid::plot(const Field& f, std::string xlabel, std::string ylabel, std::string title) const
{
	std::cout << "unset key" << std::endl;
	std::cout << "set xlabel \"" << xlabel << "\"" << std::endl;
	std::cout << "set ylabel \"" << ylabel << "\"" << std::endl;
	std::cout << "set title \"" << title << "\"" << std::endl;
	std::cout << "set xrange [" << pos_np(0) << ":" << pos_np(num_np() - 1) << "]" << std::endl;
	std::cout << "plot '-' w l" << std::endl;
	write(std::cout, f);
	std::cout << "e" << std::endl;
}

void Grid::plot(const Field& f1, const Field& f2, std::string xlabel, std::string ylabel, std::string y2label, std::string title) const
{
	std::cout << "unset key" << std::endl;
	std::cout << "set ytics nomirror" << std::endl;
	std::cout << "set y2tics nomirror" << std::endl;
	std::cout << "set xlabel \"" << xlabel << "\"" << std::endl;
	std::cout << "set ylabel \"" << ylabel << "\"" << std::endl;
	std::cout << "set y2label \"" << y2label << "\"" <<std::endl;
	std::cout << "set title \"" << title << "\"" << std::endl;
	std::cout << "set xrange [" << pos_np(0) << ":" << pos_np(num_np() - 1) << "]" << std::endl;
	std::cout << "plot '-' u 1:2 axis x1y1 w l, '' u 1:2 axis x1y2 w l" << std::endl;
	write(std::cout, f1);
	std::cout << "e" << std::endl;
	write(std::cout, f2);
	std::cout << "e" << std::endl;
}
