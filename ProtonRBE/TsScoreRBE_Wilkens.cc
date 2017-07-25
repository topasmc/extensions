// Scorer for RBE_Wilkens
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

#include "TsScoreRBE_Wilkens.hh"

#include "TsParameterManager.hh"
#include <math.h>

TsScoreRBE_Wilkens::TsScoreRBE_Wilkens(TsParameterManager* pM, TsMaterialManager* mM, TsGeometryManager* gM, TsScoringManager* scM, TsExtensionManager* eM, G4String scorerName, G4String quantity, G4String outFileName, G4bool isSubScorer)
: TsVScoreRBE_DoseLET(pM, mM, gM, scM, eM, scorerName, quantity, outFileName, isSubScorer)
{;}


TsScoreRBE_Wilkens::~TsScoreRBE_Wilkens()
{;}


TsVModelBiologicalEffect* TsScoreRBE_Wilkens::ConstructModel(G4String cellLine)
{
    return new TsModelRBE_Wilkens(cellLine, fPm, fOutputQuantity);
}


TsModelRBE_Wilkens::TsModelRBE_Wilkens(const G4String &cellLine, TsParameterManager* pM, const G4String &outputQuantity)
{
    fAlpha0 = 0.1 * (1./gray);
    fLambda = 0.02 * um/keV/gray;

    if (outputQuantity == "rbe" || outputQuantity == "survivalfraction" || outputQuantity == "rbe_x_dose") {
        G4String name = "Sc/" + cellLine + "/Alphax";
        fAlphax = pM->GetDoubleParameter(name, "perDose");
    }

    if (outputQuantity == "rbe" || outputQuantity == "survivalfraction" || outputQuantity == "rbe_x_dose" || outputQuantity == "beta") {
        G4String name = "Sc/" + cellLine + "/Betax";
        fBetax = pM->GetDoubleParameter(name, "perDoseSquare");
    }
}


G4double TsModelRBE_Wilkens::GetRBE(G4double dose, G4double LETd)
{
    if (dose == 0)
        return 0;

    return (sqrt(fAlphax*fAlphax + 4*fBetax*dose * (GetAlpha(LETd) + GetBeta(LETd)*dose)) - fAlphax) / (2*fBetax*dose);
}


G4double TsModelRBE_Wilkens::GetAlpha(G4double LETd)
{
    return fAlpha0 + (fLambda * LETd);
}


G4double TsModelRBE_Wilkens::GetBeta(G4double)
{
    return fBetax;
}
