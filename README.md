# Homeocurvature adaptation of phospholipids to pressure in deep-sea invertebrates

### Jacob R. Winnikoff, Daniel Milshteyn, Sasiri J. Vargas-Urbano, Miguel A. Pedraza, Aaron M. Armando, Oswald Quehenberger, Alexander Sodt, Richard E. Gillilan, Edward A. Dennis, Edward Lyman, Steven H. D. Haddock, and Itay Budin

**If you believe this repo is incomplete**, please message its owner or email the corresponding author.

### All files herein are distributed under the [CC BY-NC-SA (Attribution-NonCommercial-ShareAlike)](http://creativecommons.org/licenses/by-nc-sa/4.0/deed.en) license

## Guide to manuscript repo

0. `homeocurvature-2024.Rproj` - Run any of the included R code from within this project to have the relative paths set up by `here()` work correctly.

1. `01-rawdata` - Mostly SAXS profiles and lipidomics data in ASCII and spreadsheet formats. Some other input files for R scripts are included here.

2. `02-saxsanalysis` - Scripts for and results from parametric fitting of lamellar and inverted SAXS profiles. Data and Shiny code for the (ctenophore SAXS P-T eXplorer)[https://octopode.shinyapps.io/PTX_ctenos_overlay/] are also present here. 

3. `02-tidydata` - Intermediate "tidied" lipidomics data files.

4. `03-scripts` - R scripts to run all statistical analyses and visualize results. `lamellar_mcg` subdirectory contains a JScatter python implementation of the Modified Caillé Gaussian (MCG) model for multilamellar vesicle SAXS.

5. `04-mainfigs` - PDF and PNG versions of all main-text figures.

6. `04-suppfigs` - PDF and PNG versions of all SI figures.