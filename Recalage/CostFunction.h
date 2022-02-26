//---------------------------------------------------------------------------

#ifndef CostFunctionH
#define CostFunctionH

//---------------------------------------------------------------------------

#include "../CImg.h"
#include <iostream>

using namespace cimg_library;

class CostFunction
{
 protected:
	CImg<unsigned char> target, source;
	int st;

 public:
  CostFunction(const CImg<unsigned char>& target, const CImg<unsigned char>& source);
  void set_Start(int start){st = start; }
  virtual double Evaluate(double* p)=0;
  CImg<unsigned char> get_Target() {return target;} 
};

class ChamferDistance : public CostFunction
{
 private:

 public:
  ChamferDistance(const CImg<unsigned char>& target, const CImg<unsigned char>& source):CostFunction(target, source){};
  virtual double Evaluate(double* p);
};

class MutualInformation : public CostFunction
{
 private:
 int bin;
 int subsample;

 public:
  MutualInformation(const CImg<unsigned char>& target, const CImg<unsigned char>& source, int bin, int subsample);
  virtual double Evaluate(double* p);
};

#endif
