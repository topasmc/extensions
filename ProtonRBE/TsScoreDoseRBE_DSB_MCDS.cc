// Scorer for DoseRBE_DSB_MCDS
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

#include "TsScoreDoseRBE_DSB_MCDS.hh"

#include "TsParameterManager.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include <math.h>

TsScoreDoseRBE_DSB_MCDS::TsScoreDoseRBE_DSB_MCDS(TsParameterManager* pM, TsMaterialManager* mM, TsGeometryManager* gM, TsScoringManager* scM, TsExtensionManager* eM, G4String scorerName, G4String quantity, G4String outFileName, G4bool isSubScorer)
: TsVScoreBiologicalEffect(pM, mM, gM, scM, eM, scorerName, quantity, outFileName, isSubScorer)
{
    SetUnit("Gy");
}


TsScoreDoseRBE_DSB_MCDS::~TsScoreDoseRBE_DSB_MCDS()
{;}


TsVModelBiologicalEffect* TsScoreDoseRBE_DSB_MCDS::ConstructModel(G4String cellLine)
{
    return new TsModelRBE_DSB_MCDS(cellLine, fPm);
}


G4bool TsScoreDoseRBE_DSB_MCDS::ProcessHits(G4Step* aStep, G4TouchableHistory*)
{
    if (!fIsActive) {
        fSkippedWhileInactive++;
        return false;
    }

    G4double edep = aStep->GetTotalEnergyDeposit();
    if ( edep > 0. ) {
        G4int particleZ = aStep->GetTrack()->GetDefinition()->GetAtomicNumber();
//        if (particleZ>0) {
            ResolveSolid(aStep);
            G4int index = fComponent->GetIndex((G4Step*)aStep);
            TsModelRBE_DSB_MCDS* model = dynamic_cast<TsModelRBE_DSB_MCDS*>(GetModelForVoxel(index));

            G4double density = aStep->GetPreStepPoint()->GetMaterial()->GetDensity();

            G4double kineticEnergy = aStep->GetPreStepPoint()->GetKineticEnergy();
            G4double particleMass = aStep->GetTrack()->GetDefinition()->GetPDGMass();
            G4double beta = sqrt(1. - 1./(1.+kineticEnergy/particleMass)/(1.+kineticEnergy/particleMass));
            G4double Zeff = particleZ * (1. - exp(-125.*beta*pow(particleZ, -2./3.)));
            G4double x = Zeff*Zeff/beta/beta;

            G4double RBE_DSB = model->GetRBE_DSB(x);

            G4double dose = edep / (density * fSolid->GetCubicVolume());
            G4double doseRBE_DSB = dose * RBE_DSB;

            doseRBE_DSB *= aStep->GetPreStepPoint()->GetWeight();
            AccumulateHit(aStep, doseRBE_DSB);
            return true;
        }
//      return false;
//  }

    return false;
}


TsModelRBE_DSB_MCDS::TsModelRBE_DSB_MCDS(const G4String &cellLine, TsParameterManager* pM)
{
    // default = normoxic cells
    fIsAnoxic = false;
    fNormoxicA = 0.9902;
    fNormoxicB = 2.411;
    fNormoxicC = 7.32e-4;
    fNormoxicD = 1.539;

    G4String name = "Sc/" + cellLine + "/IsAnoxic";
    if (pM->ParameterExists(name)){
        fIsAnoxic = pM->GetBooleanParameter(name);
        // Anoxic cells
        if (fIsAnoxic){
            fAnoxicA = 1.502;
            fAnoxicB = 1.611;
            fAnoxicC = 1.037;
            fAnoxicD = -0.0115;
            fAnoxicE = 0.135;
            fAnoxicF = -6.096e-4;
            fAnoxicG = -8.230e-3;
            fAnoxicH = 3.047e-5;
            fAnoxicI = 3.077e-4;
        }
    }
}


G4double TsModelRBE_DSB_MCDS::GetRBE_DSB(G4double x)
{
    G4double RBE_DSB = 0.995;
    if (x>2){
        if (x<100000) {
            G4double sqrtx = sqrt(x);
            if (fIsAnoxic)
                RBE_DSB = ( fAnoxicA + sqrtx*(fAnoxicC+sqrtx*(fAnoxicE+sqrtx*(fAnoxicG+fAnoxicI*sqrtx))) ) / ( 1. + sqrtx*(fAnoxicB+sqrtx*(fAnoxicD+sqrtx*(fAnoxicF+fAnoxicH*sqrtx))) );
            else
                RBE_DSB = fNormoxicA + fNormoxicB - pow(pow(fNormoxicB, 1.-fNormoxicD) + fNormoxicC*x*(fNormoxicD-1.), 1./(1.-fNormoxicD));
        }
        else {
            if (fIsAnoxic)
                RBE_DSB = 9.93;
            else
                RBE_DSB = 3.41;
        }
    }

    return RBE_DSB;
}
