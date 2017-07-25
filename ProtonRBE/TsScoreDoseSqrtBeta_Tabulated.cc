// Scorer for DoseSqrtBeta_Tabulated
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

#include "TsScoreDoseSqrtBeta_Tabulated.hh"

#include "TsParameterManager.hh"
#include <math.h>

TsScoreDoseSqrtBeta_Tabulated::TsScoreDoseSqrtBeta_Tabulated(TsParameterManager* pM, TsMaterialManager* mM, TsGeometryManager* gM, TsScoringManager* scM, TsExtensionManager* eM, G4String scorerName, G4String quantity, G4String outFileName, G4bool isSubScorer)
: TsVScoreBiologicalEffect(pM, mM, gM, scM, eM, scorerName, quantity, outFileName, isSubScorer)
{
    SetUnit("/Gy");
}


TsScoreDoseSqrtBeta_Tabulated::~TsScoreDoseSqrtBeta_Tabulated()
{;}


TsVModelBiologicalEffect* TsScoreDoseSqrtBeta_Tabulated::ConstructModel(G4String cellLine)
{
    G4String modelName = fPm->GetStringParameter(GetFullParmName("ModelName"));

    return new TsModelBeta_Tabulated(cellLine, modelName, fPm);
}


G4bool TsScoreDoseSqrtBeta_Tabulated::ProcessHits(G4Step* aStep, G4TouchableHistory*)
{
    if (!fIsActive) {
        fSkippedWhileInactive++;
        return false;
    }

    G4double edep = aStep->GetTotalEnergyDeposit();
    if ( edep > 0. ) {
        ResolveSolid(aStep);
        G4int index = fComponent->GetIndex((G4Step*)aStep);
        TsModelBeta_Tabulated* model = dynamic_cast<TsModelBeta_Tabulated*>(GetModelForVoxel(index));

        G4int particleZ = aStep->GetTrack()->GetDefinition()->GetAtomicNumber();
        G4double kineticEnergyPerNucleon = aStep->GetPreStepPoint()->GetKineticEnergy() / aStep->GetTrack()->GetDefinition()->GetAtomicMass();

        G4double beta = model->InterpolateBeta(particleZ, kineticEnergyPerNucleon);
        if (beta < 0)
            return false;

        G4double density = aStep->GetPreStepPoint()->GetMaterial()->GetDensity();
        G4double dose = edep / (density * fSolid->GetCubicVolume());
        G4double doseSqrtBeta = dose * sqrt(beta);

        doseSqrtBeta *= aStep->GetPreStepPoint()->GetWeight();
        AccumulateHit(aStep, doseSqrtBeta);
        return true;
    }
    return false;
}


TsModelBeta_Tabulated::TsModelBeta_Tabulated(const G4String &cellLine, const G4String &modelName, TsParameterManager* pM)
{
    G4String paramPrefix = "Sc/" + cellLine + "/" + modelName + "/";

    G4String name = paramPrefix + "UseReferenceBeta";
    fUseReferenceBeta = pM->ParameterExists(name) && pM->GetBooleanParameter(name);

    if (fUseReferenceBeta) {
        name = "Sc/" + cellLine + "/Betax";
        fBetax = pM->GetDoubleParameter(name, "perDoseSquare");
    }
    else {
        name = paramPrefix + "KineticEnergyPerNucleon";
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

            if (fBeta.count(particleZ[i]) > 0) {
                G4cerr << "Topas is exiting due to a serious error in scoring." << G4endl;
                G4cerr << paramPrefix << "ParticleZ contains duplicate entries." << G4endl;
                exit(1);
            }
            fBeta[particleZ[i]] = pM->GetDoubleVector(paramPrefix + particleName[i] + "/Beta", "perDoseSquare");
        }
    }
}


G4double TsModelBeta_Tabulated::InterpolateBeta(G4int particleZ, G4double kineticEnergyPerNucleon)
{
    if (fUseReferenceBeta)
        return fBetax;

    if (fBeta.count(particleZ)==0)
        return -1;

    G4int i=0;
    while (fKineticEnergyPerNucleon[i] < kineticEnergyPerNucleon && i < fNumberOfEnergyBins-1)
        i++;
    i--;

    return fBeta[particleZ][i] + ( (kineticEnergyPerNucleon - fKineticEnergyPerNucleon[i]) / (fKineticEnergyPerNucleon[i+1] - fKineticEnergyPerNucleon[i]) * (fBeta[particleZ][i+1] - fBeta[particleZ][i]) );
}
