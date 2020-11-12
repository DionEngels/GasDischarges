/** \file lutable.h 
 *  Declarations of the two-dimensional look-up table structure mdLookupTable.
 */
#ifndef H_LUTABLE_H
#define H_LUTABLE_H

#include <string>
#include <vector>

/// a lookup table. Assumes X - SINGLE Y.
struct mdLookupTable
{
	/** Reads a lookup table from the node \a xynode and translates
	 *  the units that it reads from the [XY]Multiplicator leafs to
	 *  the units \a x_unit and \a y_unit the table should be in.
	 */
	mdLookupTable(const std::string& fname, double xscale, double yscale);
	/** Dump table to \a os. */
	void dumptab(std::ostream& os);
	// some explanation...
	void transformX(const mdLookupTable& lut, bool inv=false);
	/// \todo - JvD Move to base class ???
	void Set(unsigned n, double x0, double y0);

	/// Returns inv ? x(y) : y(x)
	double lookup(double xy0, bool inv=false) const;
	/// returns \a y and the derivative \a dydx in the point \a x
	void lookup2(double &y, double &dydx, double x, bool inv=false) const;

	const std::string& Name() const { return m_name; }
	const std::vector<double>& X() const { return m_x; }
	const std::vector<double>& Y() const { return m_y; }
	const double& X(unsigned ndx) const { return m_x[ndx]; }
	const double& Y(unsigned ndx) const { return m_y[ndx]; }
	unsigned npoints() const { return m_x.size(); }
private:
	const std::string m_name;
	double do_lookup(const std::vector<double>& args, 
			const std::vector<double>& vals, 
			double arg0) const;
	void dolookup2(
		double &res_val, double &res_deriv, 
		const std::vector<double>& args, 
		const std::vector<double>& vals,
		double arg0) const;
	std::vector<double> m_x;
	std::vector<double> m_y;
};

#endif // H_LUTABLE_H
