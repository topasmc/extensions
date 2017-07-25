RBE Scorer extensions for TOPAS

Full documentation is found at http://topas.readthedocs.io/en/latest/extensions/rbe.html

These extension scorers require TOPAS version 3.1 or later.

The directory consists of the following file types:
Scorers:
TsScoreDose* : These score quantities in the ProcessHits function like normal scorers.

RBE Scorers: These do not have a ProcessHits function and instead combine scored
properties (dose, LET, etc) to RBE or biological dose, etc.

TsV*: Base classes for the scorers

Additionally, the example directory contains an example experiment irradiation
(experiment.txt) scoring each of the available RBE scorers (rbe_scorers.txt) for
V79 cells (CellLineV79.txt). V79 cells are used because they are one of the most
studied cells and biological parameters for all models were available.
The simulations can be run with "topas run.txt" and analyzed with the provided
python script.

In order to change the experimental setup edit experiment.txt.
In order to change the cell line, provide a new cell line file and change
sv:Sc/CellLines     = 1 "CellLineV79"
in run.txt.

run.txt also controls the Prescribed dose used to calculate RBE and the output
quantity. The output quantities available depend on the RBE model.

RBE scorers are defined in rbe_scorers.txt and can be edited there. Typically,
we recommend not to run too many scorers at once as that increases memory use.
In particular, the two parameters:
ReferencedSubScorer_Dose
and
ReferencedSubScorer_LET
should be set if a dose and LET scorer already exists, otherwise each RBE scorer
will create scorers for all properties it needs, potentially resulting in
multiple parallel dose scorers.
