// Extra Class for TsYScorer

// Claculate RBE 
// Author: Hongyu Zhu
// Date: 04/19/2019

#include "TsGetRBE.hh"

#include <fstream>
#include <iostream>
#include <iomanip>
#include <cmath>

using namespace std;

TsGetRBE::TsGetRBE(double*hBinLimit, double* hBinWidth,  double* hfy,double* hdy, double hyF, double hyF_var, std::vector<double> hfy_var, std::vector<double> hdy_var, int SpecLength )
:fBinLimit(hBinLimit), fBinWidth(hBinWidth), fhfy(hfy), fhdy(hdy), yF(hyF),yF_var(hyF_var), fy_var(hfy_var), dy_var(hdy_var),fSpecLength(SpecLength)
{
    // default values of MK model
    // 10% cell survival relative to 200 kVp X-rays for HSG cells
    MKModel_alpha0 = 0.13; // Unit:Gy-1
    MKModel_beta   = 0.05; // Unit:Gy-2, 
    MKModel_rd     = 0.42; // Unit:um
    MKModel_rho    = 1.;    // Unit:g/cm3
    MKModel_y0     = 150;  // Unit:keV/um

    // default biological weighting function data file
    BioWeightFunctionDataFile = "BioWeightFuncData_interpolation.txt";


};
    
TsGetRBE::~TsGetRBE()
{};

// See reference:  H. Paganetti, et al.(1997). 
// Calculation of relative biological effectiveness for proton beams using biological weighting functions. 
// International Journal of Radiation Oncology• Biology• Physics, 37, 719-729.
 void TsGetRBE::GetRBEWithBioWeightFunction()
 {

    vector< vector<double> > BioWeightFuncData;
    
    //******************************************************************
    //                            ReadData
    //******************************************************************
    ifstream infile;
    infile.open(BioWeightFunctionDataFile);
    double x=0; // y (kev/um)
    double y=0; //r(y)


    if (!infile)    {
        std::cout<<"Error in open "<<BioWeightFunctionDataFile<<" !"<<std::endl;
        return;
    }
    

    while (!infile.eof())
    {
        infile >> x >>y;

 		vector<double> singleData;
        singleData.push_back(x);
        singleData.push_back(y);

        BioWeightFuncData.push_back(singleData);
     }


    //******************************************************************
    //                           Calculate RBE
    //******************************************************************
    double RBE = 0;
    double RBE_var =0;
    double RBE_std =0;
    double ry  = 1;

    for(int i=0; i<fSpecLength; i++)
    {
        if(fBinLimit[i]>1000 )                // if y(keV/um) set r(y)=0
            ry=0;   
        else
        {
            for(int j=0; j<BioWeightFuncData.size()-1; j++)
            { 
                if( fBinLimit[i]>=BioWeightFuncData[j][0] && fBinLimit[i]<BioWeightFuncData[j+1][0])
                    ry = BioWeightFuncData[j][1];
            }           
        }
        RBE += ry*fhdy[i]*fBinWidth[i];
        RBE_var += pow(ry*fBinWidth[i], 2)*dy_var[i];

    }
    RBE_std = sqrt(RBE_var);

    std::cout<<"******************** Get RBE with biological weight function **************************"<<std::endl;
    std::cout<<setiosflags(ios::fixed)<<setprecision(4)<<"RBE = "<<RBE;
    std::cout<<setiosflags(ios::fixed)<<setprecision(6)<<" ("<< RBE_std<<")"<<std::endl;
    std::cout<<"Default parameters:"<<endl;
    std::cout<<"1. RBE is calculated with endpoint of intestinal tolerance in mice."<<endl;
    std::cout<<"2. The value of biological weight function, r(y), was set as 0 when y> 1000 keV/um."<<endl;
    std::cout<<"***************************************************************************************\n"<<std::endl;


 }
 
 // See reference: Y. Kase, et al.(2006). 
 // Microdosimetric measurements and estimation of human cell survival for heavy-ion beams. 
 // Radiat Res, 166, 629-38.
 void TsGetRBE::GetRBEWithMKModel()
 {
    //Set parameters
	double alpha0 = MKModel_alpha0; // Unit:Gy-1
    double beta   = MKModel_beta;   // Unit:Gy-2
    double rd     = MKModel_rd;     // Unit:um
    double y0     = MKModel_y0;     // Unit:keV/um
    double rho    = MKModel_rho;    // Unit:g/cm3
    double pi     = 3.1415;
	 
	 // *************************** Calculate y_star *****************************
    //double len = sizeof(fhy) / sizeof(fhy[0]);
    double Integrate_y = 0; 
    for(int i=0; i<fSpecLength-1;i++) 
    {
        double yi = (fBinLimit[i]+fBinLimit[i+1])/2;  
        Integrate_y += (1-exp(-pow(yi,2)/pow(y0,2)))*fhfy[i]*fBinWidth[i];
    }
    double y_star = y0*y0*Integrate_y/yF;  // Unit:keV/um
 

    // ***************************** Calculate alphta *****************************
    double sub1 = y_star/(rho*pi*pow(rd,2)); // unit:keV*cm3/g/um3
    sub1 = sub1*0.16;                        // unit:J/kg = Gy
    double alpha = alpha0 + beta*sub1;       // unit:Gy-1
    //cout<<"alphta = "<<alpha<<" beta = "<<beta<<endl;

    double D10 = (sqrt(pow(alpha,2)-4*beta*log(0.1))-alpha)/(2*beta);

    // ***************************** Calculate D10_R  *****************************
    alpha = 0.19;  // Unit:Gy-1
    beta   = 0.05; // Unit:Gy-2
    double D10_R = (sqrt(pow(alpha,2)-4*beta*log(0.1))-alpha)/(2*beta);
    double RBE = D10_R/D10;

    // ************************* Calculate statistic error ************************
    double ystar_var =0;
    double alpha_var=0;
    double dose_var =0;
    double RBE_var =0;
    double RBE_std =0;

    for(int i=0; i<fSpecLength-1;i++) 
    {    
        double yi = (fBinLimit[i]+fBinLimit[i+1])/2;     
        double aa = (1-exp(-pow(yi,2)/pow(y0,2)));
        ystar_var += pow(y0*y0*aa*fBinWidth[i]/yF,2)*fy_var[i] +pow(y0*y0*aa*fhfy[i]*fBinWidth[i]/(yF*yF),2)*yF_var;
    }
    alpha_var = pow(0.16*beta/(rho*pi*pow(rd,2)),2)*ystar_var; // 0.16 convert keV*cm3/g/um3 into Gy
    double bb = pow(alpha*alpha-4*beta*log(0.1),-0.5);
    dose_var  = pow(bb*alpha/(2*beta),2)*alpha_var + pow(1/(2*beta),2)*alpha_var;
    RBE_var   = pow(D10_R/(D10*D10),2)*dose_var;
    RBE_std   = sqrt (RBE_var);

    std::cout<<"************************************ Get RBE with MK method ***************************"<<std::endl;
    std::cout<<setiosflags(ios::fixed)<<setprecision(4)<<"RBE = "<<RBE;
    std::cout<<setiosflags(ios::fixed)<<setprecision(6)<<" ("<< RBE_std<<")"<<std::endl;
    std::cout<<"Default parameters:"<<endl;
    std::cout<<"1.The reference radiation is X-ray(200 kVp) with alpha = 0.19 Gy-1 and beta = 0.05 Gy-2"<<std::endl;
    std::cout<<"2.THe bilogical end point is 10% survival of he human salivary gland (HSG) tumor cells."<<std::endl;
    std::cout<<"Parameter used in this calculation:"     <<std::endl;
    std::cout<<"alpha0="<<MKModel_alpha0<<" Gy-1; "<<"beta="<<MKModel_beta<<" Gy-2; "
             <<"rd="<<MKModel_rd<<" um; "<<"y0="<<MKModel_y0<<" keV/um; "
             <<"rho="<<MKModel_rho<<" g/cm3"  <<std::endl;
    std::cout<<"***************************************************************************************\n"<<std::endl;
	 
 }
 
