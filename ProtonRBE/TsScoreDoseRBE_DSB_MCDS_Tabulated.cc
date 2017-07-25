// Scorer for DoseRBE_DSB_MCDS_Tabulated
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

#include "TsParameterManager.hh"
#include "TsScoreDoseRBE_DSB_MCDS_Tabulated.hh"

#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include <math.h>

TsScoreDoseRBE_DSB_MCDS_Tabulated::TsScoreDoseRBE_DSB_MCDS_Tabulated(TsParameterManager* pM, TsMaterialManager* mM, TsGeometryManager* gM, TsScoringManager* scM, TsExtensionManager* eM, G4String scorerName, G4String quantity, G4String outFileName, G4bool isSubScorer)
: TsVScoreBiologicalEffect(pM, mM, gM, scM, eM, scorerName, quantity, outFileName, isSubScorer)
{
    SetUnit("Gy");

    G4String particleName = fPm->GetStringParameter(GetFullParmName("GetRBERMFForParticleNamed"));
    fParticleDefinition = G4ParticleTable::GetParticleTable()->FindParticle(particleName);

    // If didn't find the particle name, see if it exists in lower case
    if (!fParticleDefinition) {
        G4String particleNameLower = particleName;
        particleNameLower.toLower();
        fParticleDefinition = G4ParticleTable::GetParticleTable()->FindParticle(particleNameLower);

        if (!fParticleDefinition) {
            G4cerr << "Topas is exiting due to a serious error in scoring setup." << G4endl;
            G4cerr << GetFullParmName("GetLETForParticleNamed") << " refers to an unknown particle name: " << particleName << G4endl;
            exit(1);
        }
    }
}


TsScoreDoseRBE_DSB_MCDS_Tabulated::~TsScoreDoseRBE_DSB_MCDS_Tabulated()
{;}


TsVModelBiologicalEffect* TsScoreDoseRBE_DSB_MCDS_Tabulated::ConstructModel(G4String cellLine)
{
    return new TsModelRBE_DSB_MCDS_Tabulated(cellLine, fPm);
}


G4bool TsScoreDoseRBE_DSB_MCDS_Tabulated::ProcessHits(G4Step* aStep, G4TouchableHistory*)
{
    if (!fIsActive) {
        fSkippedWhileInactive++;
        return false;
    }

    if (aStep->GetTrack()->GetParticleDefinition()!=fParticleDefinition)
        return false;

    G4double edep = aStep->GetTotalEnergyDeposit();
    if (edep > 0.) {
        ResolveSolid(aStep);
        G4int index = fComponent->GetIndex((G4Step*)aStep);
        TsModelRBE_DSB_MCDS_Tabulated* model = dynamic_cast<TsModelRBE_DSB_MCDS_Tabulated*>(GetModelForVoxel(index));

        G4double density = aStep->GetPreStepPoint()->GetMaterial()->GetDensity();
        G4int particleMass = aStep->GetTrack()->GetDefinition()->GetAtomicMass();
        G4double kineticEnergyPerNucleon = aStep->GetPreStepPoint()->GetKineticEnergy() / particleMass;

        G4double RBE_DSB = model->InterpolateRBE_DSB(kineticEnergyPerNucleon);

        G4double dose = edep / (density * fSolid->GetCubicVolume());
        G4double doseRBE_DSB = dose * RBE_DSB;

        doseRBE_DSB *= aStep->GetPreStepPoint()->GetWeight();
        AccumulateHit(aStep, doseRBE_DSB);

        return true;
    }

    return false;
}


TsModelRBE_DSB_MCDS_Tabulated::TsModelRBE_DSB_MCDS_Tabulated(const G4String &cellLine, TsParameterManager* pM)
{
    G4String name = "Sc/" + cellLine + "/KE";
    fKineticEnergyPerNucleon = pM->GetDoubleVector(name, "Energy");
    fNumberOfEnergyBins = pM->GetVectorLength(name);

    name = "Sc/" + cellLine + "/DSBperGyPerCell";
    fDSBperGyPerCell = pM->GetDoubleVector(name, "perDose");

    name = "Sc/" + cellLine + "/ReferenceDSBperGyPerGbp";
    fReferenceDSBperGyPerGbp = pM->GetDoubleParameter(name, "perDose");

    name = "Sc/" + cellLine + "/GbpPerCell";
    fGbpPerCell = pM->GetUnitlessParameter(name);
    fDSBperGyPerCellx = fGbpPerCell * fReferenceDSBperGyPerGbp;
}


G4double TsModelRBE_DSB_MCDS_Tabulated::InterpolateRBE_DSB(G4double kineticEnergyPerNucleon)
{
    G4int i=0;
    while (fKineticEnergyPerNucleon[i] < kineticEnergyPerNucleon && i < fNumberOfEnergyBins-1)
        i++;
    i--;

    return (fDSBperGyPerCell[i] + (kineticEnergyPerNucleon - fKineticEnergyPerNucleon[i]) / (fKineticEnergyPerNucleon[i+1] - fKineticEnergyPerNucleon[i]) * (fDSBperGyPerCell[i+1] - fDSBperGyPerCell[i]) )  / fDSBperGyPerCellx;
}
