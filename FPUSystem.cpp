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

#include <sys/stat.h>
#include <string.h>
#include "FPUSystem.h"
#include <math.h>

FPUSystem::FPUSystem(SIM_PARAMS sim_params, OUT_PARAMS out_params, Particle * parts) {
	this->sim_params = sim_params;
	this->out_params = out_params;

	//Create output directory and files
	if (out_params.positions || out_params.momenta || out_params.energy) {
		mkdir(out_params.folder.c_str(), S_IRWXU|S_IRWXG|S_IROTH|S_IXOTH);
	}
	if (out_params.positions) {
		position_os.open((out_params.folder + "/positions.dat").c_str(), position_os.app|position_os.out);
		position_os.setf(position_os.scientific);
		position_os.precision(20);
	}
	if (out_params.momenta) {
		momenta_os.open((out_params.folder + "/momenta.dat").c_str(), momenta_os.app|momenta_os.out);
		momenta_os.setf(momenta_os.scientific);
		momenta_os.precision(20);
	}
	if (out_params.energy) {
		energy_os.open((out_params.folder + "/energy.dat").c_str(), energy_os.app|energy_os.out);
	}
	if (out_params.small_md > 0) {
		EpM_os.open((out_params.folder + "/EpM.dat").c_str(), energy_os.app|energy_os.out);
		EperMode = new double[out_params.large_md-out_params.small_md + 1];
//		FourierParts = new Particle[out_params.large_md - out_params.small_md + 1];
	}

	//Set or generate initial conditions
	if (parts != 0)
		Parts = parts;
	else
		Parts = cond_init(sim_params.NPart, sim_params.init_harm);

	Initialize();
}

FPUSystem::~FPUSystem(){
	if (out_params.positions)
		position_os.close();
	if (out_params.momenta)
		momenta_os.close();
	if (out_params.energy)
		energy_os.close();
	if (out_params.small_md)
		EpM_os.close();

	delete[] Parts;
	delete[] prevParts;
	delete[] Forces;
	delete[] EperMode;
	delete[] FourierParts;
}

void FPUSystem::Initialize(){
	currentStep = 0;
	UpdateForces();
	prevParts = new Particle[sim_params.NPart];
	for (int i = 0; i < sim_params.NPart; ++i) {
		prevParts[i].Stat.q = Parts[i].Stat.q -
							  sim_params.tau * Parts[i].Stat.p +
							  sim_params.tau * sim_params.tau * Forces[i] / 2;
	}
	if (out_params.energy)
		UpdateEnergy();
	if (out_params.small_md > 0)
		UpdateEperMode();
	printTofile();
}

void FPUSystem::UpdateForces(){
	delete[] Forces;
	Forces = new double[sim_params.NPart];
	//TODO Use second newton's law!!!
	for (int i = 1; i < sim_params.NPart-1; i++) {
		double a1 = Parts[i+1].Stat.q - Parts[i].Stat.q;
		double a2 = Parts[i].Stat.q - Parts[i-1].Stat.q;

		Forces[i] = 4 * PI2 * (Parts[i-1].Stat.q + Parts[i+1].Stat.q - 2 * Parts[i].Stat.q) +
					16 * PI4 * (sim_params.Alpha * a1 * a1 * (a1>0?1:-1) -
					sim_params.Alpha * a2 * a2 * (a1>0?1:-1));
	}
	Forces[0] = -Forces[1];
}

void FPUSystem::UpdatePos(){
	Particle *tmp = new Particle[sim_params.NPart];
	for (int i = 1; i < sim_params.NPart-1; ++i) {
		tmp[i].Stat.q = 2 * Parts[i].Stat.q - prevParts[i].Stat.q +
						sim_params.tau * sim_params.tau * Forces[i];
	}
	delete[] prevParts;
	prevParts = Parts;
	Parts = tmp;
}

void FPUSystem::UpdateMomenta(){
	for (int i = 1; i < sim_params.NPart-1; ++i) {
		Parts[i].Stat.p = (Parts[i].Stat.q - prevParts[i].Stat.q)/sim_params.tau +
						   sim_params.tau * Forces[i] / 2;
	}
}

void FPUSystem::UpdateEperMode() {
	delete[] FourierParts;
	FourierParts = new Particle[out_params.large_md - out_params.small_md + 1];
	float tmp = 0;
	int cursor = 0;
	while (out_params.small_md + cursor <= out_params.large_md) {
		int k = out_params.small_md + cursor;
		for (int i = 1; i < sim_params.NPart + 1; ++i) {
			tmp = sin(i * k * M_PI / (sim_params.NPart-1));
			FourierParts[cursor].Stat.q += Parts[i].Stat.q * tmp;
			FourierParts[cursor].Stat.p += Parts[i].Stat.p * tmp;
		}
		double omega = sin(k * M_PI / (sim_params.NPart - 1));
		EperMode[cursor] = (FourierParts[cursor].Stat.p * FourierParts[cursor].Stat.p)+
				4 * PI2 * omega * omega *	// <-------------- # aqui!!!
				(FourierParts[cursor].Stat.q * FourierParts[cursor].Stat.q);//(sim_params.NPart - 1)/(sim_params.NPart - 1);
		cursor++;
	}
}

double FPUSystem::V(double s, float alpha) {
	double s_2 = s * s;
	return 8 * PI4 * s_2 + s_2 * s * 64 / 3 * alpha * PI4 * PI4;
}

void FPUSystem::UpdateEnergy() {
	last_T = 0;
	last_V = 0;
	for (int i = 0; i < sim_params.NPart-1; ++i) {
			last_T += Parts[i].Stat.p * Parts[i].Stat.p * 2 * PI2;
			last_V += V(Parts[i+1].Stat.q - Parts[i].Stat.q, sim_params.Alpha);
		}
	}

void FPUSystem::printTofile() {
	if (out_params.positions) {
		for (int i = 0; i < sim_params.NPart; ++i) {
			position_os << Parts[i].Stat.q << "\t";
		}
		position_os << std::endl;
	}

	if (out_params.momenta) {
		for (int i = 0; i < sim_params.NPart; ++i) {
			momenta_os << Parts[i].Stat.p << "\t";
		}
		momenta_os << std::endl;
	}

	if (out_params.energy) {
		energy_os << last_T << "\t" << last_V << "\t";
		energy_os << last_T + last_V << std::endl;
	}

	if (out_params.small_md) {
		for (int k = 0; k < out_params.large_md - out_params.small_md + 1; ++k) {
			EpM_os << EperMode[k] << '\t';
		}
		EpM_os << std::endl;
	}
}

void FPUSystem::StepSim() {
	currentStep++;
	UpdatePos();
	UpdateForces();
	if ((currentStep % out_params.ratio) == 0) {
		if (out_params.momenta||out_params.energy||out_params.small_md)
			UpdateMomenta();
		if (out_params.energy)
			UpdateEnergy();
		if (out_params.small_md) {
			UpdateEperMode();
		}
		printTofile();
	}
}

Particle * FPUSystem::cond_init(int Npart, short harmonic){
	Particle *ret = new Particle[Npart];
	VECTOR st; st.p = 0;

	double scale = M_PI * harmonic / (Npart-1);
	for (int i = 1; i < Npart-1; ++i) {
		st.q = sin(scale * i) / 100;
		ret[i] = Particle(st);
	}
	ret[0].Stat.q = 0; ret[Npart-1].Stat.q = 0;
	return ret;
}
