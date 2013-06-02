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

#include <stdlib.h>
#include <iostream>
#include <sysexits.h>
#include <getopt.h>
#include "FPUSystem.h"

using namespace std;

int main(int argc, char **argv) {

	//Define default values for parameters
	SIM_PARAMS fpu_params;
	fpu_params.NPart 	=	3;
	fpu_params.init_harm=	1;
	fpu_params.tau		=	0.001;
	fpu_params.Alpha	=	0;

	OUT_PARAMS out_params;
	out_params.momenta	=	0;
	out_params.positions= 	0;
	out_params.energy	=	0;
	out_params.folder 	=	".";
	out_params.ratio	=	1;
	out_params.small_md =	0;
	out_params.large_md	=	out_params.small_md;

	long	STEPS		=	1000;

	//Define command line parameters to be parsed by getopt
	static struct option long_options[] =
	{
			{"N-particles", required_argument, 0,	'N'},
			{"initial-mode", required_argument, 0,	'i'},
			{"delta-t", required_argument, 0,		't'},
			{"steps", required_argument, 0,			'S'},
			{"alpha", required_argument, 0,			'a'},
			{"output-folder", required_argument, 0,	'O'},
			{"initial-cond", required_argument, 0,	'I'},
			{"saving-ratio", required_argument, 0,	'r'},
			{"smaller-mode", required_argument, 0,	'm'},
			{"larger-mode", required_argument, 0, 	'M'},

			{"save-positions", no_argument, &(out_params.positions), 1},
			{"save-momenta", no_argument, &(out_params.momenta), 1},
			{"save-energy", no_argument, &(out_params.energy), 1},
			{0,0,0,0}
	};

	//Do the actual command line parsing
	int opt_index = 0;
	int c;
	bool setIC = false;
	string ICfile = "";
	//loop through all given options reading parameters
	while ((c = getopt_long(argc, argv, "qpEN:m:t:S:a:O:i:r:I:M:", long_options, &opt_index)) != -1) {
		switch (c) {
		case 0:
			break;
		case 'N':
			fpu_params.NPart = (int)strtol(optarg,0,10) + 2;
			break;
		case 'i':
			fpu_params.init_harm = (int)strtol(optarg,0,10);
			break;
		case 't':
			fpu_params.tau = strtof(optarg, 0);
			break;
		case 'S':
			STEPS = strtol(optarg, 0, 10);
			break;
		case 'a':
			fpu_params.Alpha = strtof(optarg, 0);
			break;
		case 'O':
			out_params.folder = optarg;
			break;
		case 'q':
			out_params.positions = 1;
			break;
		case 'p':
			out_params.momenta = 1;
			break;
		case 'E':
			out_params.energy = 1;
			break;
		case 'I':
			setIC = true;
			ICfile = optarg;
			break;
		case 'r':
			out_params.ratio = strtol(optarg,0, 10);
			break;
		case 'm':
			out_params.small_md = strtol(optarg,0, 10);
			if (out_params.large_md < out_params.small_md) {
				out_params.large_md = out_params.small_md;
			}
			break;
		case 'M':
			out_params.large_md = strtol(optarg, 0, 10);
			if (out_params.small_md <= 0) {
				out_params.small_md = 1;
			}
			if (out_params.large_md < out_params.small_md) {
				out_params.small_md = out_params.large_md;
			}
			break;
		default:
			abort();
		}
	}

	//Create a System instance and read initial conditions
	//if asked.
	FPUSystem *FPU;
	if (!setIC)
		FPU = new FPUSystem(fpu_params, out_params);
	else {
		VECTOR st; st.p = 0;
		Particle *icParts = new Particle[fpu_params.NPart];

		ifstream in(ICfile.c_str());
		if (!in.good()) {
			cerr << "Error opening initial conditions file: " << ICfile << endl;
			abort();
		}
		for (int i = 0; i < fpu_params.NPart; ++i) {
			in >> st.q;
			icParts[i] = Particle(st);
		}
		for (int i = 0; i < fpu_params.NPart; ++i) {
			in >> icParts[i].Stat.p;
		}
		in.close();

		FPU = new FPUSystem(fpu_params, out_params, icParts);
	}

	//Main simulation loop
	//loop STEPS times
	for (int i = 0; i < STEPS; i++) {
		FPU->StepSim();

		cout << ((float)i / STEPS * 100) << "\t%\r";
	}
	FPU->~FPUSystem();

	return EX_OK;
}
