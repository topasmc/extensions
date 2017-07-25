// Scorer for RBE_MKMLET
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

#include "TsScoreRBE_MKMLET.hh"

#include "TsParameterManager.hh"
#include <math.h>
#include <iostream>
#include <fstream>

TsScoreRBE_MKMLET::TsScoreRBE_MKMLET(TsParameterManager* pM, TsMaterialManager* mM, TsGeometryManager* gM, TsScoringManager* scM, TsExtensionManager* eM, G4String scorerName, G4String quantity, G4String outFileName, G4bool isSubScorer)
: TsVScoreRBE_DoseLET(pM, mM, gM, scM, eM, scorerName, quantity, outFileName, isSubScorer)
{;}


TsScoreRBE_MKMLET::~TsScoreRBE_MKMLET()
{;}


TsVModelBiologicalEffect* TsScoreRBE_MKMLET::ConstructModel(G4String cellLine)
{
    return new TsModelRBE_MKMLET(cellLine, fPm, fOutputQuantity);
}


TsModelRBE_MKMLET::TsModelRBE_MKMLET(const G4String &cellLine, TsParameterManager* pM, const G4String &)
{
    fGamma0 = 0.229 * gray*um*um*um/keV;

    G4String name = "Sc/" + cellLine + "/AlphaBetaRatiox";
    fAlphaBetax = pM->GetDoubleParameter(name, "Dose");

    name = "Sc/" + cellLine + "/Alphax";
    fAlphax = pM->GetDoubleParameter(name, "perDose");

    name = "Sc/" + cellLine + "/Betax";
    fBetax = pM->GetDoubleParameter(name, "perDoseSquare");

    name = "Sc/" + cellLine + "/Alpha0MKM";
    fAlpha0 = pM->GetDoubleParameter(name, "perDose");

    name = "Sc/" + cellLine + "/DomainDiameter";
    fDomainDiameter = pM->GetDoubleParameter(name, "Length");
}


G4double TsModelRBE_MKMLET::GetRBE(G4double dose, G4double LETd)
{
    if (dose == 0)
        return 0;

    return sqrt(0.25*fAlphaBetax*fAlphaBetax + GetAlpha(LETd)/fAlphax*fAlphaBetax*dose + GetBeta(LETd)/fBetax*dose*dose) / dose - fAlphaBetax/(2*dose);
}


G4double TsModelRBE_MKMLET::GetAlpha(G4double LETd)
{
    G4double gamma = fGamma0/fDomainDiameter/fDomainDiameter * LETd;
    return fAlpha0 + gamma*fBetax;
}


G4double TsModelRBE_MKMLET::GetBeta(G4double)
{
    return fBetax;
}
