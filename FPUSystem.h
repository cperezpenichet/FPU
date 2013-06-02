//==============================================================================
//Copyright (C) 2013 Carlos Pérez-Penichet <cperezpenichet@gmail.com>
//This program is free software: you can redistribute it and/or modify it 
//under the terms of the GNU General Public License version 3, as published 
//by the Free Software Foundation.
//
//This program is distributed in the hope that it will be useful, but 
//WITHOUT ANY WARRANTY; without even the implied warranties of 
//MERCHANTABILITY, SATISFACTORY QUALITY, or FITNESS FOR A PARTICULAR 
//PURPOSE.  See the GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License along 
//with this program.  If not, see <http://www.gnu.org/licenses/>.
//==============================================================================

#ifndef SYSTEM_H_
#define SYSTEM_H_

#include <fstream>
#include <math.h>
#include "Particle.h"

typedef struct {
	int 	NPart;
	float	Alpha;
	int		init_harm;
	float	tau;
} SIM_PARAMS;

typedef struct {
	int 		positions;
	int 		momenta;
	int 		energy;
	std::string	folder;
	short		ratio;
	short 		small_md;
	short		large_md;
} OUT_PARAMS;

class FPUSystem {
public:
	FPUSystem(SIM_PARAMS sim_params, OUT_PARAMS out_flags, Particle * parts = 0);
	~FPUSystem();

	void StepSim();

	//Generate initial conditions based on a given harmonic number
	static Particle * cond_init(int Npart, short harmonic);

private:
	void Initialize();
	void UpdateForces();
	void UpdatePos();
	void UpdateMomenta();
	void UpdateEnergy();
	void UpdateEperMode();

	void printTofile();
	static const double PI2 = (M_PIl* M_PIl);
	static const double PI4 = ((M_PIl* M_PIl)*(M_PIl* M_PIl));

	SIM_PARAMS sim_params;
	OUT_PARAMS out_params;

	Particle *Parts;
	Particle *prevParts;
	Particle *FourierParts;

	double	 *Forces;

	std::ofstream position_os;
	std::ofstream momenta_os;
	std::ofstream energy_os;
	std::ofstream EpM_os;

	double last_T;
	double last_V;
	double *EperMode;

	long currentStep;

	double V(double s, float alpha);
};

#endif /* SYSTEM_H_ */
