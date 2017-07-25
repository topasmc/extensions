// Scorer for RBE_TabulatedAlphaBeta
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


#include "TsScoreRBE_TabulatedAlphaBeta.hh"

#include "TsParameterManager.hh"

TsScoreRBE_TabulatedAlphaBeta::TsScoreRBE_TabulatedAlphaBeta(TsParameterManager* pM, TsMaterialManager* mM, TsGeometryManager* gM, TsScoringManager* scM, TsExtensionManager* eM,
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
    InstantiateSubScorer("DoseAlpha_Tabulated", outFileName, "DoseAlpha");
    InstantiateSubScorer("DoseSqrtBeta_Tabulated", outFileName, "DoseSqrtBeta");
}


TsScoreRBE_TabulatedAlphaBeta::~TsScoreRBE_TabulatedAlphaBeta()
{;}


TsVModelBiologicalEffect* TsScoreRBE_TabulatedAlphaBeta::ConstructModel(G4String cellLine)
{
    return new TsModelRBE_TabulatedAlphaBeta(cellLine, fPm, fOutputQuantity);
}


G4int TsScoreRBE_TabulatedAlphaBeta::CombineSubScorers()
{
    TsVBinnedScorer* doseScorer = dynamic_cast<TsVBinnedScorer*>(GetSubScorer("Dose"));
    TsVBinnedScorer* doseAlphaScorer = dynamic_cast<TsVBinnedScorer*>(GetSubScorer("DoseAlpha"));
    TsVBinnedScorer* doseSqrtBetaScorer = dynamic_cast<TsVBinnedScorer*>(GetSubScorer("DoseSqrtBeta"));

    std::vector<G4double> normalizedDose = NormalizeDose(doseScorer->fFirstMomentMap);

    if (fOutputQuantity == "rbe") {
        for (unsigned index = 0; index<fFirstMomentMap.size(); index++) {
            TsModelRBE_TabulatedAlphaBeta* model = dynamic_cast<TsModelRBE_TabulatedAlphaBeta*>(GetModelForVoxel(index));
            G4double alpha = model->GetAlpha(doseScorer->fFirstMomentMap[index], doseAlphaScorer->fFirstMomentMap[index]);
            G4double beta = model->GetBeta(doseScorer->fFirstMomentMap[index], doseSqrtBetaScorer->fFirstMomentMap[index]);
            fFirstMomentMap[index] = model->GetRBE(normalizedDose[index], alpha, beta);
        }
    }
    else if (fOutputQuantity == "alpha") {
        for (unsigned index = 0; index<fFirstMomentMap.size(); index++) {
            TsModelRBE_TabulatedAlphaBeta* model = dynamic_cast<TsModelRBE_TabulatedAlphaBeta*>(GetModelForVoxel(index));
            fFirstMomentMap[index] = model->GetAlpha(doseScorer->fFirstMomentMap[index], doseAlphaScorer->fFirstMomentMap[index]);
        }
    }
    else if (fOutputQuantity == "beta") {
        for (unsigned index = 0; index<fFirstMomentMap.size(); index++) {
            TsModelRBE_TabulatedAlphaBeta* model = dynamic_cast<TsModelRBE_TabulatedAlphaBeta*>(GetModelForVoxel(index));
            fFirstMomentMap[index] = model->GetBeta(doseScorer->fFirstMomentMap[index], doseSqrtBetaScorer->fFirstMomentMap[index]);
        }
    }
    else if (fOutputQuantity == "survivalfraction") {
        for (unsigned index = 0; index<fFirstMomentMap.size(); index++) {
            TsModelRBE_TabulatedAlphaBeta* model = dynamic_cast<TsModelRBE_TabulatedAlphaBeta*>(GetModelForVoxel(index));
            G4double alpha = model->GetAlpha(doseScorer->fFirstMomentMap[index], doseAlphaScorer->fFirstMomentMap[index]);
            G4double beta = model->GetBeta(doseScorer->fFirstMomentMap[index], doseSqrtBetaScorer->fFirstMomentMap[index]);
            fFirstMomentMap[index] = model->GetSurvivalFraction(normalizedDose[index], alpha, beta);
        }
    }
    else if (fOutputQuantity == "rbe_x_dose") {
        for (unsigned index = 0; index<fFirstMomentMap.size(); index++) {
            TsModelRBE_TabulatedAlphaBeta* model = dynamic_cast<TsModelRBE_TabulatedAlphaBeta*>(GetModelForVoxel(index));
            G4double alpha = model->GetAlpha(doseScorer->fFirstMomentMap[index], doseAlphaScorer->fFirstMomentMap[index]);
            G4double beta = model->GetBeta(doseScorer->fFirstMomentMap[index], doseSqrtBetaScorer->fFirstMomentMap[index]);
            fFirstMomentMap[index] = normalizedDose[index] * model->GetRBE(normalizedDose[index], alpha, beta);
        }
    }

    return 0;
}


TsModelRBE_TabulatedAlphaBeta::TsModelRBE_TabulatedAlphaBeta(const G4String &cellLine, TsParameterManager* pM, const G4String &outputQuantity)
{
    if (outputQuantity == "rbe" || outputQuantity == "rbe_x_dose") {
        G4String name = "Sc/" + cellLine + "/Alphax";
        fAlphax = pM->GetDoubleParameter(name, "perDose");

        name = "Sc/" + cellLine + "/Betax";
        fBetax = pM->GetDoubleParameter(name, "perDoseSquare");
    }
}


G4double TsModelRBE_TabulatedAlphaBeta::GetSurvivalFraction(G4double dose, G4double alpha, G4double beta)
{
    return exp(-alpha*dose - beta*dose*dose);
}


G4double TsModelRBE_TabulatedAlphaBeta::GetAlpha(G4double dose, G4double doseAlpha)
{
    if (dose == 0)
        return 0;

    return doseAlpha / dose;
}


G4double TsModelRBE_TabulatedAlphaBeta::GetBeta(G4double dose, G4double doseSqrtBeta)
{
    if (dose == 0)
        return 0;

    return doseSqrtBeta*doseSqrtBeta / (dose*dose);
}


G4double TsModelRBE_TabulatedAlphaBeta::GetRBE(G4double dose, G4double alpha, G4double beta)
{
    if (dose == 0)
        return 0;

    return (sqrt(fAlphax*fAlphax + 4*fBetax*dose * (alpha + beta*dose)) - fAlphax) / (2*fBetax*dose);
}
