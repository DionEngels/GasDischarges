#ifndef H_NPF_CROSSSEC_H
#define H_NPF_CROSSSEC_H

class CrossSection
{
  public:
	CrossSection(double eps_th)
	 : m_eps_th(eps_th)
	{}
	virtual ~CrossSection() {}
	double operator()(double eps) const
	{
		return (eps<eps_th()) ? 0.0 : get_above_th(eps);
	}
	double eps_th() const { return m_eps_th; }
  protected:
	virtual double get_above_th(double eps) const=0;
  private:
	double m_eps_th;
};

// return a constant cross section above a certain threshold
class ConstCrossSec : public CrossSection
{
  public:
	ConstCrossSec(double eps_th, double value)
	 : CrossSection(eps_th), m_value(value)
	{}
  protected:
	virtual double get_above_th(double eps) const
	{ 
		return m_value; 
	}
  private:
	const double m_value;
};


#endif // H_NPF_CROSSSEC_H

