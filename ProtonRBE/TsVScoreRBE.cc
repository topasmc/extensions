// Extra Class for AllScorers
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

#include "TsVScoreRBE.hh"

#include "TsParameterManager.hh"

TsVScoreRBE::TsVScoreRBE(TsParameterManager* pM, TsMaterialManager* mM, TsGeometryManager* gM, TsScoringManager* scM, TsExtensionManager* eM,
                         G4String scorerName, G4String quantity, G4String outFileName, G4bool isSubScorer)
: TsVScoreBiologicalEffect(pM, mM, gM, scM, eM, scorerName, quantity, outFileName, isSubScorer)
{
    fOutputQuantity = fPm->GetStringParameter(GetFullParmName("OutputQuantity"));
    fOutputQuantity.toLower();
    // OutputQuantity aliases
    if (fOutputQuantity == "sf")
        fOutputQuantity = "survivalfraction";
    if (fOutputQuantity == "rwd") // RBE weighted dose (RWD)
        fOutputQuantity = "rbe_x_dose";


    fIsSimultaneousExposure = true;
    if (fPm->ParameterExists(GetFullParmName("SimultaneousExposure")))
        fIsSimultaneousExposure = fPm->GetBooleanParameter(GetFullParmName("SimultaneousExposure"));

    fPrescribedDose = fPm->GetDoubleParameter(GetFullParmName("PrescribedDose"), "Dose");

    fPrescribedDoseMetric = "Max";
    if (fPm->ParameterExists(GetFullParmName("PrescribedDoseMetric")))
        fPrescribedDoseMetric = fPm->GetStringParameter(GetFullParmName("PrescribedDoseMetric"));
    fPrescribedDoseMetric.toLower();
    if (fPrescribedDoseMetric != "max" && fPrescribedDoseMetric != "mean" && fPrescribedDoseMetric != "d90") {
        G4cerr << "Topas is exiting due to a serious error in scoring setup." << G4endl;
        G4cerr << "Invalid " << GetFullParmName("PrescribedDoseMetric") << " parameter value: " << fPrescribedDoseMetric << G4endl;
        G4cerr << "Valid values are: Max, Mean and D90." << G4endl;
        exit(1);
    }

    if (fPm->ParameterExists(GetFullParmName("PrescribedDoseStructure"))) {
        fPrescribedDoseStructureName = fPm->GetStringParameter(GetFullParmName("PrescribedDoseStructure"));

        if (!fIsSimultaneousExposure) {
            G4cerr << "Topas is exiting due to a serious error in scoring setup." << G4endl;
            G4cerr << "Scorer \"" << GetName() << "\" is using both SimultaneousExposure" << G4endl;
            G4cerr << "and PrescribedDoseStructure parameters, but these are mutually exclusive." << G4endl;
            exit(1);
        }
    }
}


TsVScoreRBE::~TsVScoreRBE()
{;}


std::vector<G4double> TsVScoreRBE::NormalizeDose(const std::vector<G4double> &dose)
{
    std::vector<G4double> normalizedDose;

    if (fIsSimultaneousExposure) {

        G4double scaleFactor = 1.0;

        if (fPrescribedDoseStructureName.isNull())
            scaleFactor = GetScaleFactor(dose);
        else {
            G4int prescribedDoseStructureID = fComponent->GetStructureID(fPrescribedDoseStructureName);

            std::vector<G4double> structureDose;
            for (unsigned index = 0; index<dose.size(); index++)
                if (fComponent->IsInNamedStructure(prescribedDoseStructureID, index))
                    structureDose.push_back(dose[index]);

            scaleFactor = GetScaleFactor(structureDose);
        }

        for (std::vector<G4double>::const_iterator it = dose.begin(); it != dose.end(); ++it)
            normalizedDose.push_back(*it * scaleFactor);
    }
    else {
        normalizedDose = std::vector<G4double>(dose.size(), fPrescribedDose);
    }

    return normalizedDose;
}


G4double TsVScoreRBE::GetScaleFactor(const std::vector<G4double> &dose)
{
    G4double refDose = 0;

    if (fPrescribedDoseMetric == "max") {
        refDose = *std::max_element(dose.begin(), dose.end());
    }
    else if (fPrescribedDoseMetric == "mean") {
        refDose = 0;
        for (std::vector<G4double>::const_iterator it = dose.begin(); it != dose.end(); ++it)
            refDose += *it;
        refDose /= dose.size();
    }
    else if (fPrescribedDoseMetric == "d90") {
        std::vector<G4double> dose_copy = dose;

        const size_t i90 = (1-0.9) * dose.size();
        std::nth_element(dose_copy.begin(), dose_copy.begin() + i90, dose_copy.end());
        refDose = dose_copy[i90];
    }

    if (refDose == 0) {
        G4cerr << "Topas is exiting due to a serious error in scoring." << G4endl;
        G4cerr << "Scorer \"" << GetName() << "\" attempted to use a null reference dose." << G4endl;
        exit(1);
    }

    return fPrescribedDose / refDose;
}
