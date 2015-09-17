# DEMENT
Decomposition Model of Enzymatic Traits maintained by Steven Allison

DEMENT is an agent-based numerical simulation model, open-source, coded in R, and available for all operating systems. The model uses microbial trait values and correlations to predict rates of organic matter decomposition. The most recent model version represents traits related to enzyme kinetics, enzyme production, growth efficiency, stoichiometry, and responses to moisture availability. During a model run, a large number (>100) of bacterial and fungal taxa are allowed to compete on a spatial grid representing the surface of decomposing organic material. Each taxon possesses a suite of physiological traits that are assigned based on trait correlations. Enzymes produced by the microbial taxa interact locally with substrates to generate monomers that are available for uptake. The model predicts microbial community composition by simulating the abundances of the initial taxa at a daily time step. Enzymatic degradation is a Michaelis-Menten process with Arrhenius temperature sensitivity functions built into Vmax and Km kinetic parameters.

Instructions to get started

The current set of files is configured to run in batch mode, for example on a computer cluster, but you can also run in batch mode on your desktop computer from the command line.
1) Download all of the files into the same directory.
2) Create a subdirectory called "params".
3) Move the "params1.txt" file into the params subdirectory.
4) Modify any parameter values you wish (parameter descriptions are in "ParamDescription.txt").
5) Open a terminal or command prompt in the directory you created with the model files.
6) At the prompt, type (without the quotes): "Rscript --no-restore DEMENTBatch.R 000000000000 1".
