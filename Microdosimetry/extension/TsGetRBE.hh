// Claculate RBE 
// Author: Hongyu Zhu
// Date: 04/19/2019

#ifndef TsGetRBE_hh
#define TsGetRBE_hh

#include <stdint.h>
#include <vector>
#include <string>


using namespace std;

class TsGetRBE 
{
public:
    TsGetRBE(double*BinLimit, double* hBinWidth, double* hfy, double* hdy, double hyF, double hyF_var, std::vector<double> hfy_var, std::vector<double>hdy_var, int SpecLength);
    ~TsGetRBE();

	void GetRBEWithMKModel();
    void GetRBEWithBioWeightFunction();
    void SetMKModel_alpha0(double num ){ MKModel_alpha0 = num; }
    void SetMKModel_beta  (double num ){ MKModel_beta   = num; }
    void SetMKModel_rho   (double num ){ MKModel_rho    = num; }
    void SetMKModel_rd    (double num ){ MKModel_rd     = num; }
    void SetMKModel_y0    (double num ){ MKModel_y0     = num; }
    void SetBioWeightFunctionDataFile(string fileName ){ BioWeightFunctionDataFile= fileName; }
    
private:

   double* fBinLimit;
   double* fBinWidth;
   double* fhy;
   double* fhfy;
   double* fhdy;
   double* fhydy;
   double  yF;
   double  yF_var;
   std::vector<double> fy_var;
   std::vector<double> dy_var;
   
   int fSpecLength;

   double MKModel_alpha0,MKModel_beta, MKModel_rho, MKModel_rd, MKModel_y0;
   string BioWeightFunctionDataFile;


};



#endif
