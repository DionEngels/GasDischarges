#include <iostream>
#include <cmath>

#include "physconst.h"
#include "phivar.h"

//The physical radius of the grid in m.
const double grid_size = 0.0125;
//The maximum of the initialization function of Te
const double Te_max=12000;
//The minimum of the initialization function of Te
const double Te_min=11000;
//The maximum of the initialization function of Th
const double Th_max=400.0;
//The minimum of the initialization function of Th
const double Th_min=300.0;
//The background particle density
const double n0_dens=4e22;
//The ion density at the wall. Supposed to be low to
//represent recombination.
const double ion_wall_dens=1e12;
//The maximum of the initialization function of the ion density.
const double ion_dens_max=1e18;
//The minimum of the initialization function of the ion density.
const double ion_dens_min=1e12;
//The average density of the filling gas.
const double ave_dens=n0_dens;

//B.Broks, 9-1-04:
//A function which computes your local buffer density based on an average density of the plasma
//gas. This avearge density is the argument of your function call. It will require the input of
//following fields:
//An electron temperature Field called Te.
//A heavy particle temperature Field called Th.
//An electron density Field called ion_dens.
//A heavy particle density Field neutr_dens.
//It will modify the contents of neutr_dens so the neutr_dens becomes the needed buffer density.
void calculate_buffer_density(	double ave_dens,
				const Grid &grid,
				PhiVariable &neutr_dens,
				const PhiVariable &ion_dens,
				const PhiVariable &Te,
				const PhiVariable &Th)
{
	//Step zero: the volume of the grid
	//Iterate over the whole gird, adding the local volume tho the running total.
	double grid_volume=0.0;
	unsigned i=0;
	for(i=0; i<grid.num_np(); ++i)
		grid_volume +=grid.vol_np(i);

	//Step one: compute the point where we have the highest.
	//nonbuffer density. While we are at it, we also compute many
	//of our total particles we have used.

	double highest_local_nonbuffer_pressure=0.0;
	//The amount of particles used.
	double total_particles_used=0.0;
	for (i=0; i<grid.num_np(); ++i)
	{
		//p=n_e k_b T_e (elec) +n_e k_b T_h (ion)
		double nonbuffer_press_local =	PhysConst::k_b*Te[i]*ion_dens[i]+
						PhysConst::k_b*Th[i]*ion_dens[i];
		total_particles_used+=ion_dens[i]*grid.vol_np(i);
		if (nonbuffer_press_local>highest_local_nonbuffer_pressure)
			highest_local_nonbuffer_pressure=nonbuffer_press_local;
	}
	//Okay, now we know what the maximum pressure is.
	//Step two: We will now assign particles to each point to remove
	//the pressure difference.

	for (i=0; i<grid.num_np(); ++i)
	{
		double nonbuffer_press_local =	PhysConst::k_b*Te[i]*ion_dens[i]+
						PhysConst::k_b*Th[i]*ion_dens[i];
		double pressure_to_compensate=
			highest_local_nonbuffer_pressure-nonbuffer_press_local;
		double particles_used_to_compensate=
			pressure_to_compensate/(PhysConst::k_b*Th[i])*grid.vol_np(i);
		//Again, this costs us particles.
		total_particles_used+=particles_used_to_compensate;
		neutr_dens[i]=pressure_to_compensate/(PhysConst::k_b*Th[i]);

	}
	//Final step:Distribute the rest.
	//Compute the average value of 1/T_h. The particles are distributed proportional to 1/Th.
	double one_over_t_ave=0;
	for (i=0; i<grid.num_np(); ++i)
 		one_over_t_ave+=grid.vol_np(i)/grid_volume/Th[i];
	//Ditribute the paricles.
	for (i=0; i<grid.num_np(); ++i)
		neutr_dens[i]+=1.0/one_over_t_ave/Th[i]*(ave_dens-total_particles_used/grid_volume);
}

int main(int argc, char **argv)
{
	//The grid. The first argument is the amount of grid points, the
	//second argument is the position of the left edge, and the third
	//argument is the position of the right edge.
	Grid grid(16, 0.0, grid_size, Grid::Cylindrical);
	//the velocity is defined on the cv boundaries. The value is zero
	Field flow_vel(grid.num_ew(), 0.0);
	//Te. We have two Neumann conditions.
	NeumannBndCond Te_left(0,0);
	NeumannBndCond Te_right(0,0);
	PhiVariable Te(grid, flow_vel, Te_left, Te_right);
	//Th. We cool the right wall to 300 K.
	NeumannBndCond Th_left(0, 0);
	DirichletBndCond Th_right(Th_min);
	PhiVariable Th(grid, flow_vel, Th_left, Th_right);
	//Ion Density. We have recombination at the right wall, represented by the low
	//ion density value.
	NeumannBndCond ion_dens_left(0, 0);
	DirichletBndCond ion_dens_right(ion_wall_dens);
	PhiVariable ion_dens(grid, flow_vel, ion_dens_left, ion_dens_right);
	//neutr_dens
	NeumannBndCond neutr_dens_left(0, 0);
	NeumannBndCond neutr_dens_right(0, 0);
	PhiVariable neutr_dens(grid, flow_vel, neutr_dens_left, neutr_dens_right);

	// Ready
	return 0;
}
