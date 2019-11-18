#ifndef TsYScorer_hh
#define TsYScorer_hh

#include "TsVNtupleScorer.hh"
#include "TsTrackerHit.hh"
#include "vector"

class TsYScorer : public TsVNtupleScorer
{
public:
    TsYScorer(TsParameterManager* pM, TsMaterialManager* mM, TsGeometryManager* gM, TsScoringManager* scM, TsExtensionManager* eM,
                G4String scorerName, G4String quantity, G4String outFileName, G4bool isSubScorer);
    
    virtual ~TsYScorer();


    G4bool ProcessHits(G4Step*,G4TouchableHistory*);
    void AbsorbResultsFromWorkerScorer(TsVScorer* workerScorer);
    void UserHookForEndOfRun();

    void PushBackData(G4double* Edep);
    void GetSpectrum();
    void Calculatefy(std::vector<G4double> dataVector);
    void GetStatisticInfo(G4int Binindex, G4double variable);
    void GetErrorPropagation();
    void InitializeMicrodosimetricSpectrum();
    void InitializeStatistic();
    void OutputResult();
    void GetRBE();

    

private:
    TsTrackerHitsCollection* fHitsCollection;
    G4bool IncludeYF, IncludeYD;
    G4double TEradius;
    G4double SVradius, SVheight;
    G4double MeanPathLength;
    G4double ylimitL, ylimitU;
    G4double CenterX, CenterY, CenterZ;
    G4int    GeoNo;
    G4int    NumberOfHistoriesInRun;
    G4int    fStatisticUpdateFrequency;
    G4double fMeanChordLength;
    G4bool fGetRBEWithBioWeightFunction;
	G4bool fGetRBEWithMKModel;
    G4bool fGetSecondariesContribution;
    G4bool fGetStatisticInfo;
    G4int  fSpectrumUpdateTimes;

    std::vector<G4double> yVector;
    std::vector<G4double * > yVector_Particle;



    // statisic
    std::vector<G4double> fFirstMomentMap;
    std::vector<G4double> fSecondMomentMap;
    std::vector<G4double> fCountMap;
    std::vector<G4double> fVariance;
    std::vector<G4double> fStandardDeviation;

    G4double yF, yD;
    G4double yF_var, yF_std;
    G4double yD_var, yD_std;
    std::vector<G4double> ydy_var, ydy_std;
    std::vector<G4double> yfy_var, yfy_std;
    std::vector<G4double>  dy_var,  dy_std;


protected:

    void AccumulateEvent();

    G4double fy;
    G4double fy_z0; // e-
    G4double fy_z1; // proton, deuteron, triton
    G4double fy_z2; // He3, alpha, He6
    G4double fy_z3; // Li7
    G4double fy_z4; // Be7, Be9, 
    G4double fy_z5; // B10, B11
    G4double fy_z6; // C11, C12
    G4double fy_z_; // other

    const G4double yBinMagnitude         =  5.0;
    const G4double yBinMagnitudeInterval = 20.0;
    const G4int    yBinNum               = 100; // yBinNum == yBinMagnitude*yBinMagnitudeInterval
    G4double *hfy;
    G4double *hdy;
    G4double *hyfy;
    G4double *hydy;
    G4double **hfy_particle;
    G4double *BinLimit;
    G4double *BinWidth;

};

#endif
