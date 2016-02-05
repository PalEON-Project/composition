composition
===========

Code for fitting composition models across the entire domain.  master.sh is the controller script that indicates how all processing is done.  Data files are downloaded from the Paleon Wiki to the 'dataDir' specified in the config file. Output is placed in the 'outputDir' specified in the config file.

If you want to reproduce the work as a non-Paleon participant you can find the input data files in the public domain as follows:
western subdomain: https://portal.lternet.edu/nis/mapbrowse?packageid=msb-paleon.2.0 is equivalent to western-0.6-2.csv
northeastern subdomain: SetTreeComp_Northeast_Level1_v1.0.csv is equivalent to 1372polygons_v0.9-1.csv and SetTreeComp_Northeast_Level1_v1.0.zip contains the shapefile info in 1372polygons_v0.9-1.*
ohio subdomain:  SetTreeComp_Ohio_Level1_v1.0.csv is equivalent to OH471polygons_v0.2-1.csv and SetTreeComp_Ohio_Level1_v1.0.zip contains the shapefile infor in OH471polygons_v0.2-1.*

Current master branch creates the v0.4 product posted on DataOne/NIS and used in the PLOS ONE paper. Eastern product was not refit as the input data did not change from v0.3. Simply copied eastern v0.3 to v0.4.

The version0.3 branch was used in early 2015 to create the v0.3 product.

The version0.2 branch is intended as static snapshot of the code used to create the version 0.2 of the statistical composition product, as of Apr. 1, 2014.

The version0.1 branch is intended as static snapshot of the code used to create the version 0.1 of the statistical composition product, as of Feb. 28, 2014.

lindgren_nu1_logsigma2 is branch that contains code for both CAR and Lindgren(nu = 1) as well as for cross-validation, with joint sampling of hyperparams and alphas. It was merged into master in fall 2014.

See *labbook* for notes on statistical model development and model runs.
