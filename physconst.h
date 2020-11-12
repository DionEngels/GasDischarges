#ifndef H_PHYSCONST_H
#define H_PHYSCONST_H

/** A set of fundamental constants and other physical values has been
 *  made available in the namespace PhysConst. Namespaces help to
 *  disambiguate identifiers with such short names as 'e'. In order
 *  to access these constants, prepend their names with PhysConst::,
 *  as in 'PhysConst::e'.
 */
namespace PhysConst {

	/// Our good-old pi
	const double pi = 3.1415926535897932384626433833;

	/// The permittivity of vacuum [F/m]
	const double epsilon0 = 8.854187e-12;

	/// The electron mass [kg]
	const double me = 9.10938189e-31; 

	/// The elementary charge unit [C]
	const double e = 1.6022e-19;

	/// The Boltzmann constant [J/K]
	const double k_b = 1.3806503e-23;
	
	//Planck's Constant
	const double h_planck =6.62606877e-34;
    
	/// The atomic mass unit [kg]
	const double AMU = 1.66053873e-27;

	/** The value of the energy unit eV [J].
	 *  This is obtained as e*1V; as a result, the numerical
	 *  value is equal to the elementary charge unit.
	 */
	const double eV = e;

} // end of namespace PhysConst

#endif // H_PHYSCONST_H

