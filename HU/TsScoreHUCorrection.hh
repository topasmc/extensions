//
// ********************************************************************
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * TOPAS collaboration.                                             *
// * Use or redistribution of this code is not permitted without the  *
// * explicit approval of the TOPAS collaboration.                    *
// * Contact: Joseph Perl, perl@slac.stanford.edu                     *
// *                                                                  *
// ********************************************************************
//

#ifndef TsScoreHUCorrection_hh
#define TsScoreHUCorrection_hh

#include "TsVBinnedScorer.hh"

#include "G4EmCalculator.hh"

class G4Material;

class TsScoreHUCorrection : public TsVBinnedScorer
{
public:
	TsScoreHUCorrection(TsParameterManager* pM, TsMaterialManager* mM, TsGeometryManager* gM, TsScoringManager* scM, TsExtensionManager* eM,
					   G4String scorerName, G4String quantity, G4String outFileName, G4bool isSubScorer);
    
	
	virtual ~TsScoreHUCorrection();
	
	G4bool ProcessHits(G4Step*,G4TouchableHistory*);

private:
	G4EmCalculator fEmCalculator;
	G4Material* fWater;
	G4ParticleDefinition* fProton;
    G4int fFlag;
    G4String fFileName;
    std::ofstream fFile;
};
#endif
