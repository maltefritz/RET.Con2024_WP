# RET.Con2024_WP

# Overview

This repository contains additional information and data relating to the contribution to the "7. RET.Con" from the authors  M. Fritz, J. Freißmann and I. Tuschy. It serves to allow for reproduction of the obtained results.

# The conference

You can find information regarding the conference on the respective [homepage](https://www.hawk.de/de/hochschule/fakultaeten-und-standorte/fakultaet-ressourcenmanagement/profil/nwf](https://www.hs-nordhausen.de/forschung/inret/ret-con/)). There, you will also be able to find the book of abstracts and eventually the full paper.

# Description of contents

## Optimization

This folder contains the files for carrying out the combined investment and dispatch optimization as well as for the dispatch optimization. This requires auxiliary files such as those for the economic functions and helpers. It also contains a postprocessing file for processing the results. In addition to control files, the folders for the various district heating systems (primary network, sub network, and 4GDH network) are also stored.

### Input data

This folder contains the necessary input data for the combined investment and dispatch optimization. Generally, this includes a JSON file of the constant parameters (no variation throughout the observed period) and a CSV file of the time dependent data.

### Output data

This folder contains the output data of the combined investment and dispatch optimization. For each setup, there are the unit commitment time series, as well as key parameters and unit cost. Additionally, the results are visualized in plots, that allow an even more detailed analysis than provided in the paper.
