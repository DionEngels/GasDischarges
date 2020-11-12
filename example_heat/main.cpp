#include <iostream>
#include <fstream>

#include "phivar.h"

int main( int argc, char **argv )
{
    const unsigned N = 16;
	const double xa = 0.0;
	const double xb = 10;
	const double Tr = 400.0;
	const double V = 0.0;
	const double lambda = 1.0;
	const double Sc = 10.0;
	const double Sp = 0.0;
	const double beta_cv = 1.0;
	
	Grid grid( N, xa, xb, Grid::Cartesian );

	Field w( grid.num_ew(), V ); 

	NeumannBndCond temp_left( 0,0 );
	DirichletBndCond temp_right( Tr );

	PhiVariable temp( grid, w, temp_left, temp_right );

	temp.beta_cv = beta_cv;
	for( unsigned i=0; i<grid.num_np(); ++i)
	{
		temp.lambda[i] = lambda;
		temp.sc[i] = Sc;
		temp.sp[i] = Sp;
	}
  
	temp.Update();

	grid.write( std::cerr, temp );
	std::ofstream ofs_T("T.dat");
	grid.write( ofs_T, temp);
	std::ofstream ofs_flux("heat_flux_density.dat");
	grid.write( ofs_flux, temp.flux_density());

	grid.plot( temp, "position (m)", "temperature (K)" );

	return 0;
} 
