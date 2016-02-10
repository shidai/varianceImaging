// Dynamic spectrum simulator
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <fftw3.h>
//#include <gsl/gsl_rng.h>
//#include <gsl/gsl_randist.h>
#include "dynSpecSim.h"

int main (int argc, char* argv[])
{
	controlStruct control;
	acfStruct acfStructure;

	char fname[1024];   // read in parameter file
	//char oname[1024];   // output file name
	char pname[1024];   // file to plot
	char dname[1024];   // graphics device
	int i;
	int n = 1;  // number of dynamic spectrum to simulate

	int plotMode = 0;  // plot existing dynamic spectrum, off by default
	control.noplot = 0;  // show dynamic spectrum while simulationg, on by default

	// read options
	for (i=0;i<argc;i++)
	{
		if (strcmp(argv[i],"-f") == 0)
		{
			strcpy(fname,argv[++i]);
			printf ("Parameters are in %s\n", fname);
		}
		//else if (strcmp(argv[i],"-o")==0)
		//{
		//	strcpy(oname,argv[++i]);
			//printf ("Dynamic spectrum is output into %s\n", oname);
		//}
		else if (strcmp(argv[i],"-n")==0)
		{
			n = atoi (argv[++i]);
			printf ("Number of dynamic spectrum to simulate %d\n", n);
		}
		else if (strcmp(argv[i],"-p")==0) // just plot 
		{
			strcpy(pname,argv[++i]);
			plotMode = 1;
			printf ("Plotting exiting dynamic spectrum.\n");
		}
		else if (strcmp(argv[i],"-noplot")==0)
		{
			control.noplot = 1; // Don't show dynamic spetrum while simulationg
		}
		else if (strcmp(argv[i],"-dev")==0)
		{
			strcpy(dname,argv[++i]);
			printf ("Graphic device: %s\n", dname);
		}
	}

	if (plotMode==0)
	{
		// Simulate dynamic spectrum
		// read parameters
		initialiseControl(&control);
		readParams (fname,dname,n,&control);
		//readParams (fname,oname,dname,n,&control);
		printf ("Finished reading parameters.\n");

		// simulate dynamic spectra
		calculateScintScale (&acfStructure, &control);
		//calculateScintScale (&acfStructure, &control, seed);

		if (control.noplot==0 && n == 1)
		{
			// plot while simulating
			heatMap (&acfStructure, dname);
		}

		qualifyVar (&acfStructure, &control);

		// deallocate memory
		deallocateMemory (&acfStructure);
	}
	else
	{
		plotDynSpec(pname, dname);
	}

	return 0;
}

