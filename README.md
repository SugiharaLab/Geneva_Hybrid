# lake_geneva_DO

Main and supplemental analyses for

A hybrid empirical and parametric approach for managing ecosystem complexity: water quality in Lake Geneva under nonstationary futures.
Ethan Deyle1,2,*,†, Damien Bouffard3,†, Victor Frossard4, Robert Schwefel5,6, John Melack7, George Sugihara2

1 Department of Biology, Boston University, Boston, MA 02215
2 Scripps Institution of Oceanography, University of California San Diego, La Jolla, CA 92093
3 Department of Surface Waters – Research and Management, Eawag - Swiss Federal Institute of Aquatic Science and Technology, Kastanienbaum, Switzerland 6047
4 UMR 42 CARRTEL, Université Savoie Mont Blanc, Le Bourget du Lac, France 73376
5 Ecology, Evolution, and Marine Biology, University of California Santa Barbara, Santa Barbara, CA 93106
6 Environmental Engineering Institute, École Polytechnique Fédérale de Lausanne, Lausanne, Switzerland
7 Bren School of Environmental Science & Management, University of California Santa Barbara, Santa Barbara, CA 93106
* Corresponding authors Ethan Deyle, Damien Bouffard, George Sugihara
† Contributed equally as primary authors

## File Structure

Project directory should have scripts in the root folder as well as folders labeled "DATA", "FUNCTIONS", "INPUTS", and "OUTPUTS. Analyses are all organized in Rmarkdown files in the main project directoy. Below each is described briefly.

### PNAS_SI_notebook.Rmd

This Rmarkdown contains the core calculations of the manuscript and recreates main text and supplement plots (without full annotation and manicure).

### Monad_S_map_demo.Rmd

This Rmarkdown contains demonstration of S-map approach using a Monad model.

### analysis_chl.Rmd

This Rmarkdown contains additional empirical dynamic modeling analysis of "chl" (Chlorophyll) to support core calculations in PNAS_SI_notebook.Rmd

### analysis_po4_epi.Rmd

This Rmarkdown contains additional empirical dynamic modeling analysis of "PO4_epi" (Orthophosphate in the epilimnion) to support core calculations in PNAS_SI_notebook.Rmd

### figures_for_manuscript.Rmd

This Rmarkdown reproduces plots in PNAS_SI_notebook.Rmd with added manicure for the use in the official manuscript.