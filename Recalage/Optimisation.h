//---------------------------------------------------------------------------
#ifndef OptimisationH
#define OptimisationH
//---------------------------------------------------------------------------

#include "CostFunction.h"
#include "Globdef.h"
#include "Simann.h"
#include "../CImg.h"

using namespace cimg_library;

void local1(int,double,int,double*,double*,int*,double*,double*,CostFunction&);
CImg<double> QuasiNewton( CostFunction& cost, CImg<double> &t, double &angle );

CImg<double> SimulatedAnnealing( CostFunction& CF, CImg<double> &t, double &angle );

int global(double*,double*,int,int,int,FILE*,int,TOMB_nx21,int*,double*,int*,CostFunction&);
CImg<double> StochasticClustering( CostFunction& CF, CImg<double> &t, double &angle, double delta_translation, double delta_rotation );

void update(double*,int,double*,double*,double*,int*,int,double);
void fun(double*,double*,int,double*,double*,CostFunction&);
void urdmn(TOMB_nx101,int);

//---------------------------------------------------------------------------
#endif