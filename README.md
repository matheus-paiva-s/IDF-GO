# IDF_GO
Sample script of methodologies for constructing “PEU” and “PP” IDF curves for GEV and Gumbel distributions

This set of R scripts performs the calculation of
Intensity–Duration–Frequency (IDF) curves from annual maximum
precipitation series using two distinct methodologies: PP and PEU.

FILES

funcoes.R
Contains only auxiliary functions shared by the PP and PEU
methodologies. It includes IDF model functions, an objective function,
error metrics, and supporting routines. No code is executed directly
within this file.

PP.R
Execution script for the PP methodology. It is responsible for data
preprocessing, IDF model fitting, optimization, performance evaluation,
and result generation. It relies on the functions defined in funcoes.R.

PEU.R
Execution script for the PEU methodology, with its own model-fitting
logic. It shares the auxiliary functions from funcoes.R, allowing
direct comparison with the PP methodology.

DATA

Series Selecionadas
This folder contains the datasets used as input in the dissertation,
including treated historical rainfall series and annual maximum
precipitation intensities by duration.


ORGANIZATION

The separation between auxiliary functions and execution scripts is
intended to improve methodological clarity, reproducibility, and code
maintainability.

