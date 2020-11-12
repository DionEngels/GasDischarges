#include "lutable.h"
#include "basics.h"
#include <cassert>
#include <fstream>

#define WARN_CLIPPING

void mdLookupTable::transformX(const mdLookupTable& lut, bool inv)
{
	if (lut.npoints()==0)
	{
		npfFatalError("required transformation table empty!");
	}
	for (unsigned i=0;i<npoints();i++)
	{
		m_x[i] = lut.lookup(X(i),inv);
	}
}

void mdLookupTable::Set(unsigned n, double x0, double y0)
{
	m_x[n]=x0;
	m_y[n]=y0;
}

double mdLookupTable::do_lookup(const std::vector<double>& args, const std::vector<double>& vals, double arg0) const
{
	assert( args.size() == vals.size() );
	const int ndx_last=npoints()-1;
	int ndx=1;

#ifdef WARN_CLIPPING
	if ((arg0<args[0])||(arg0>args[args.size()-1])) 
	{
		std::cout << "value (" << arg0 << ") is outside range ("
			<< args[0] << "," << args[args.size()-1] << ") of " 
			<< m_name << " lutable" << std::endl;
	}
#endif
	
	if (arg0>args[ndx_last])
	{
		return vals[ndx_last];
	}
	while ( (arg0>args[ndx]) && (ndx<ndx_last) )
	{
		ndx++;
	}
	const double d1=(arg0-args[ndx-1])/(args[ndx]-args[ndx-1]);
	return d1*vals[ndx]+(1.0-d1)*vals[ndx-1];
}

double mdLookupTable::lookup(double xy0, bool inv) const
{
	return inv ? do_lookup( Y(), X(), xy0)
		: do_lookup( X(), Y(), xy0);
}

void mdLookupTable::dolookup2(double &res_val, double &res_deriv, 
	const std::vector<double>& args, 
	const std::vector<double>& vals,
	double arg0) const
{
	assert( args.size() == vals.size() );
#ifdef WARN_CLIPPING
	if ((arg0<args[0])||(arg0>args[args.size()-1])) 
	{
		std::cout << "value (" << arg0 << ") is outside range ("
			<< args[0] << "," << args[args.size()-1] << ") of " 
			<< m_name << " lutable" << std::endl;
	}
#endif
	

	if (arg0>args[npoints()-1])
	{
		res_deriv=0.0;
		res_val=vals[npoints()-1];
	}
	else {
		unsigned i=1;
		while ((arg0>args[i])&&(i<npoints()-1))
			i++;
		res_deriv=(vals[i]-vals[i-1])/(args[i]-args[i-1]);
		res_val=res_deriv*(arg0-args[i-1])+vals[i-1];
	}
}

void mdLookupTable::lookup2(double &res_val, double &res_deriv, double xy0, bool inv) const
{
	if (inv) 
		dolookup2( res_val, res_deriv, Y(), X(), xy0);
	else
		dolookup2( res_val, res_deriv, X(), Y(), xy0);
}

mdLookupTable::mdLookupTable(const std::string& fname, double xscale, double yscale)
 : m_name(fname)
{
	std::ifstream ifs(fname.c_str());
	if (!ifs)
	{
		const std::string msg("Could not open file '" + fname + "' for reading.");
		npfFatalError(msg.c_str());
	}
	std::cout << "Opened file " << fname << std::endl;
	std::string comment; std::getline(ifs,comment);
	std::cout << "Comment: " << comment << std::endl;
	while (1)
	{
		double xval,yval;
		ifs >> xval >> yval;
		if (ifs.eof())
		{
			break;
		}
		m_x.push_back(xval*xscale);
		m_y.push_back(yval*yscale);
	}
	std::cout << "Read " << npoints() << " values from file " << Name() << std::endl;
}

void mdLookupTable::dumptab(std::ostream& os)
{
	for (unsigned i=0;i<npoints();i++)
	{
		os << X(i) << "\t" << Y(i) << std::endl;
	}
}
