// Scorer for DoseAlpha_Tabulated
//
// ********************************************************************
// *                                                                  *
// * This  code implementation is an extension to developments by the *
// * TOPAS collaboration.                                             *
// * This extension  is an freely  available in  accordance  with the *
// * freeBSD license:                                                 *
// * Copyright (c) <2015>, <Harald Paganetti>                         *
// * All rights reserved.                                             *
// * Redistribution    and   use in   source and   binary    forms,   *
// * with or without modification, are permitted provided that the    *
// * following conditions are met:                                    *
// *                                                                  *
// *                                                                  *
// * 1. Redistributions of source code must retain the above          *
// * copyright notice, this                                           *
// * list of conditions and the following disclaimer.                 *
// * 2. Redistributions in binary form must reproduce the above       *
// * copyright notice, this list of conditions and the following      *
// * disclaimer in the documentation and/or other materials provided  *
// * with the distribution.                                           *
// *                                                                  *
// * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND           *
// * CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES,      *
// * INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF         *
// * MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE         *
// * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR             *
// * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,     *
// * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT *
// * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF *
// * USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED  *
// * AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT      *
// * LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING   *
// * IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF   *
// * THE POSSIBILITY OF SUCH DAMAGE.                                  *
// *                                                                  *
// * The views and conclusions contained in the software and          *
// * documentation are those of the authors and should not be         *
// * interpreted as representing official policies, either expressed  *
// * or implied, of the FreeBSD Project.                              *
// *                                                                  *
// * Contacts: Jan Schuemann, jschuemann@mgh.harvard.edu              *
// *           Harald Paganetti, hpaganetti@mgh.harvard.edu           *
// *                                                                  *
// ********************************************************************
//

#include "TsScoreDoseAlpha_Tabulated.hh"

#include "TsParameterManager.hh"

TsScoreDoseAlpha_Tabulated::TsScoreDoseAlpha_Tabulated(TsParameterManager* pM, TsMaterialManager* mM, TsGeometryManager* gM, TsScoringManager* scM, TsExtensionManager* eM, G4String scorerName, G4String quantity, G4String outFileName, G4bool isSubScorer)
: TsVScoreBiologicalEffect(pM, mM, gM, scM, eM, scorerName, quantity, outFileName, isSubScorer)
{
    SetUnit("");
}


TsScoreDoseAlpha_Tabulated::~TsScoreDoseAlpha_Tabulated()
{;}


TsVModelBiologicalEffect* TsScoreDoseAlpha_Tabulated::ConstructModel(G4String cellLine)
{
    G4String modelName = fPm->GetStringParameter(GetFullParmName("ModelName"));

    return new TsModelAlpha_Tabulated(cellLine, modelName, fPm);
}


G4bool TsScoreDoseAlpha_Tabulated::ProcessHits(G4Step* aStep, G4TouchableHistory*)
{
    if (!fIsActive) {
        fSkippedWhileInactive++;
        return false;
    }

    G4double edep = aStep->GetTotalEnergyDeposit();
    if ( edep > 0. ) {
        ResolveSolid(aStep);
        G4int index = fComponent->GetIndex((G4Step*)aStep);
        TsModelAlpha_Tabulated* model = dynamic_cast<TsModelAlpha_Tabulated*>(GetModelForVoxel(index));

        G4int particleZ = aStep->GetTrack()->GetDefinition()->GetAtomicNumber();
        G4double kineticEnergyPerNucleon = aStep->GetPreStepPoint()->GetKineticEnergy() / aStep->GetTrack()->GetDefinition()->GetAtomicMass();

        G4double alpha = model->InterpolateAlpha(particleZ, kineticEnergyPerNucleon);
        if (alpha < 0)
            return false;

        G4double density = aStep->GetPreStepPoint()->GetMaterial()->GetDensity();
        G4double dose = edep / (density * fSolid->GetCubicVolume());
        G4double doseAlpha = dose * alpha;

        doseAlpha *= aStep->GetPreStepPoint()->GetWeight();
        AccumulateHit(aStep, doseAlpha);
        return true;
    }
    return false;
}


TsModelAlpha_Tabulated::TsModelAlpha_Tabulated(const G4String &cellLine, const G4String &modelName, TsParameterManager* pM)
{
    G4String paramPrefix = "Sc/" + cellLine + "/" + modelName + "/";

    G4String name = paramPrefix + "KineticEnergyPerNucleon";
    fKineticEnergyPerNucleon = pM->GetDoubleVector(name, "Energy");
    fNumberOfEnergyBins = pM->GetVectorLength(name);

    name = paramPrefix + "ParticleName";
    G4String* particleName = pM->GetStringVector(name);
    fNumberOfParticleNames = pM->GetVectorLength(name);

    name = paramPrefix + "ParticleZ";
    G4int* particleZ = pM->GetIntegerVector(name);
    if (pM->GetVectorLength(name) != fNumberOfParticleNames) {
        G4cerr << "Topas is exiting due to a serious error in scoring." << G4endl;
        G4cerr << paramPrefix << " has different vector lengths" << G4endl;
        G4cerr << "for ParticleName and ParticleZ parameters." << G4endl;
        exit(1);
    }

    for (G4int i=0; i<fNumberOfParticleNames; i++) {

        if (fAlpha.count(particleZ[i]) > 0) {
            G4cerr << "Topas is exiting due to a serious error in scoring." << G4endl;
            G4cerr << paramPrefix << "ParticleZ contains duplicate entries." << G4endl;
            exit(1);
        }
        fAlpha[particleZ[i]] = pM->GetDoubleVector(paramPrefix + particleName[i] + "/Alpha", "perDose");
    }
}


G4double TsModelAlpha_Tabulated::InterpolateAlpha(G4int particleZ, G4double kineticEnergyPerNucleon)
{
    if (fAlpha.count(particleZ)==0)
        return -1;

    G4int i=0;
    while (fKineticEnergyPerNucleon[i] < kineticEnergyPerNucleon && i < fNumberOfEnergyBins-1)
        i++;
    i--;

    return fAlpha[particleZ][i] + ( (kineticEnergyPerNucleon - fKineticEnergyPerNucleon[i]) / (fKineticEnergyPerNucleon[i+1] - fKineticEnergyPerNucleon[i]) * (fAlpha[particleZ][i+1] - fAlpha[particleZ][i]) );
}
