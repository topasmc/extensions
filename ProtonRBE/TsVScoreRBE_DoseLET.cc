// Extra Class for scorers
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


#include "TsVScoreRBE_DoseLET.hh"

#include "TsParameterManager.hh"

TsVScoreRBE_DoseLET::TsVScoreRBE_DoseLET(TsParameterManager* pM, TsMaterialManager* mM, TsGeometryManager* gM, TsScoringManager* scM, TsExtensionManager* eM,
                         G4String scorerName, G4String quantity, G4String outFileName, G4bool isSubScorer)
: TsVScoreRBE(pM, mM, gM, scM, eM, scorerName, quantity, outFileName, isSubScorer)
{
    if (fOutputQuantity == "rbe")
        SetUnit("");
    else if (fOutputQuantity == "alpha")
        SetUnit("/Gy");
    else if (fOutputQuantity == "beta")
        SetUnit("/Gy2");
    else if (fOutputQuantity == "survivalfraction")
        SetUnit("");
    else if (fOutputQuantity == "rbe_x_dose")
        SetUnit("Gy");
    else {
        G4cerr << "Topas is exiting due to a serious error in scoring setup." << G4endl;
        G4cerr << "No output quantity " << fOutputQuantity << " defined." << G4endl;
        G4cerr << "Valid output quantities are: RBE, Alpha, Beta, SurvivalFraction and RBE_x_Dose (BiologicalDose)." << G4endl;
        exit(1);
    }

    if (!fPm->ParameterExists(GetFullParmName("PreCalculateStoppingPowerRatios")))
        fPm->AddParameter("b:" + GetFullParmName("PreCalculateStoppingPowerRatios"), "\"True\"");
    InstantiateSubScorer("DoseToWater", outFileName, "Dose");
    InstantiateSubScorer("ProtonLET", outFileName, "LET");
}


TsVScoreRBE_DoseLET::~TsVScoreRBE_DoseLET()
{;}


G4int TsVScoreRBE_DoseLET::CombineSubScorers()
{
    TsVBinnedScorer* doseScorer = dynamic_cast<TsVBinnedScorer*>(GetSubScorer("Dose"));
    TsVBinnedScorer* letScorer = dynamic_cast<TsVBinnedScorer*>(GetSubScorer("LET"));
    G4double density = 1.0 * g/cm3;  // ProtonLET scores LET per unit density

    std::vector<G4double> normalizedDose = NormalizeDose(doseScorer->fFirstMomentMap);

    if (fOutputQuantity == "rbe") {
        for (unsigned index = 0; index<fFirstMomentMap.size(); index++) {
            TsVModelRBE_DoseLET* model = dynamic_cast<TsVModelRBE_DoseLET*>(GetModelForVoxel(index));
            fFirstMomentMap[index] = model->GetRBE(normalizedDose[index], density*letScorer->fFirstMomentMap[index]);
        }
    }
    else if (fOutputQuantity == "alpha") {
        for (unsigned index = 0; index<fFirstMomentMap.size(); index++) {
            TsVModelRBE_DoseLET* model = dynamic_cast<TsVModelRBE_DoseLET*>(GetModelForVoxel(index));
            fFirstMomentMap[index] = model->GetAlpha(density*letScorer->fFirstMomentMap[index]);
        }
    }
    else if (fOutputQuantity == "beta") {
        for (unsigned index = 0; index<fFirstMomentMap.size(); index++) {
            TsVModelRBE_DoseLET* model = dynamic_cast<TsVModelRBE_DoseLET*>(GetModelForVoxel(index));
            fFirstMomentMap[index] = model->GetBeta(density*letScorer->fFirstMomentMap[index]);
        }
    }
    else if (fOutputQuantity == "survivalfraction") {
        for (unsigned index = 0; index<fFirstMomentMap.size(); index++) {
            TsVModelRBE_DoseLET* model = dynamic_cast<TsVModelRBE_DoseLET*>(GetModelForVoxel(index));
            fFirstMomentMap[index] = model->GetSurvivalFraction(normalizedDose[index], model->GetAlpha(density*letScorer->fFirstMomentMap[index]), model->GetBeta(density*letScorer->fFirstMomentMap[index]));
        }
    }
    else if (fOutputQuantity == "rbe_x_dose") {
        for (unsigned index = 0; index<fFirstMomentMap.size(); index++) {
            TsVModelRBE_DoseLET* model = dynamic_cast<TsVModelRBE_DoseLET*>(GetModelForVoxel(index));
            fFirstMomentMap[index] = normalizedDose[index] * model->GetRBE(normalizedDose[index], density*letScorer->fFirstMomentMap[index]);
        }
    }

    return 0;
}


G4double TsVModelRBE_DoseLET::GetSurvivalFraction(G4double dose, G4double alpha, G4double beta)
{
    return exp(-alpha*dose - beta*dose*dose);
}
