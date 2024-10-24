# Lévy graphical models
R code related to the paper *Lévy graphical models*.

The file <code>functions.R</code> contains all necesary functions for simulating and estimating Lévy graphical models on trees.

The file <code>simulations.R</code> contains the code to reproduce the simulation study of the paper.  It simulates a multivariate Lévy process with Hüsler--Reiss type jumps. Note that running this code may take several minutes. The code relies on the file <code>functions.R</code>.

The file <code>data_analysis.R</code> performs the analysis of stock price data, which is described in the paper. It also reproduces all figures related to this application. The code relies on the file <code>functions.R</code>. The data is contained in the file <code>data/data_all.csv</code>.
