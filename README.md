[![DOI](https://zenodo.org/badge/214175219.svg)](https://zenodo.org/badge/latestdoi/214175219)

# Pressure solution modelling using IODP data from Expedition 375, Site U1520
Code for recreating and exploring the figures produced for a study on shallow (<1 km below sea floor) stylolites and pressure solution seaward of the Hikurangi Margin (resubmitted to Tectonics, 2020).

There are two scripts, one produces stress estimates and performs the modelling of pressure solution ([u1520_ps_modelling.m](./u1520_ps_modelling.m)), the other calculates the coefficient of variation from the clustering of faults and stylolites ([u1520_styl_fault_clustering.m](./u1520_styl_fault_clustering.m)).

Data is available from IODP.org or the BGS repository No. 130545 (DOI: [10.5285/af50c4d2-fc87-4052-a41b-226d0590c12b](https://www.bgs.ac.uk/services/NGDC/citedData/catalogue/af50c4d2-fc87-4052-a41b-226d0590c12b.html) ). A zipped datafile is also available on request from the author ( LeahHR[at]cardiff.ac.uk ) for direct use with the scripts once unzipped to a folder named 'data' in the same directory.

By downloading this code, you accept that I am not liable for any effects these scripts may have on anything due to your running them.

To save figures, [export_fig](https://uk.mathworks.com/matlabcentral/fileexchange/23629-export_fig) is receommended.

If you find any errors or have any problems running these scripts, please feel free to raise an issue in the GitHub repository or email me at the address above.
