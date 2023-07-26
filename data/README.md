# NZ Wastewater Modelling: Data

We rely on two core sources of data for our modelling, these are:
- [ESR's wastewater sampling data](https://github.com/ESR-NZ/covid_in_wastewater)
- [MOH's reported case data](https://github.com/minhealthnz/nz-covid-data)

The relevant files are saved in */data/raw/*. *ww_data_all.csv*, *ww_national.csv*, and *sites.csv* were obtained from ESR, while *covid-case-counts-moh.csv* was obtained from MOH.

Additionally this folder contains .csv files describing the assumed onset-to-shedding distribution (*dist_ww_shedding.csv*, estimated by [Hewitt et al](https://doi.org/10.1016/j.watres.2021.118032)), and the assumed onset-to-reporting distribution (*dist_esr_onset.csv*, estimated directly from NZ COVID-19 line-data).

Finally, *borderWorkerInfections.csv* is a subset of the data available [here](https://github.com/michaelplanknz/modelling-ba5-in-nz). This is used in Figure 4 of the paper to show cumlative infections in a cohort of regularly tested border workers for comparison with model outputs.