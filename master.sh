#!/bin/bash
# controller script for fitting the PalEON settlement era composition model
# both the western PLS data and the eastern township data
# note that these steps are intended for use on UNIX-like machines and will need to be modified for Windows (and possibly for Mac OS X)

# modify the contents of the config file to reflect the data versions to be used, relevant directories, and parameters of the MCMC
source config

# get cookie with Wiki authentication info
wget --post-data="u=${WIKI_USERNAME}&p=${WIKI_PASSWORD}&sectok=b7649cb05e87dc0f45d9d91851a7b097&id=start&r=1&do=login" --save-cookies=${dataDir}/my-cookies.txt --keep-session-cookies https://paleon.geography.wisc.edu/doku.php/dw__login

export cookieArgs="--load-cookies=${dataDir}/my-cookies.txt --save-cookies=${dataDir}/my-cookies.txt --keep-session-cookies"

# not working... problems with quotes I think
function wgetWiki() {
    wget --load-cookies=${dataDir}/my-cookies.txt --save-cookies=${dataDir}/my-cookies.txt --keep-session-cookies $1 -O $2
}

########################################################################
# download meta info    ------------------------------------------------
########################################################################

# get taxon conversion table
cd $projectDir
wget $cookieArgs "https://paleon.geography.wisc.edu/lib/exe/fetch.php/data_and_products%3Bpublic_data%3Blevel3s_v${productVersion}.csv" -O level3s_v${productVersion}.csv


########################################################################
# download western data ------------------------------------------------
########################################################################

cd $dataDir
# note this doesn't work because of authentication issues, so navigate to here via browser instead
wget $cookieArgs "https://paleon.geography.wisc.edu/lib/exe/fetch.php/data_and_products%3Bpublic_data%3Bwesterncompv${westernVersion}.csv" -O westerncomp-${westernVersion}.csv



mv westerncompv${westernVersion}.csv western-${westernVersion}.csv
rm -f western.csv
ln -s western-${westernVersion}.csv western.csv

cd $projectDir

########################################################################
# preprocess western data ----------------------------------------------
########################################################################

./build_western.R >& log.build_western_${productVersion}-${uniqueRunID} &
# this creates 'westernData.Rda'

########################################################################
# fit Bayesian composition model to western data -----------------------
########################################################################

./fit_western.R >& log.fit_western_${productVersion}-${uniqueRunID} &
# this creates 'PLScomposition_western_${productVersion}.nc'

########################################################################
# download eastern township data -------------------------------------
########################################################################

cd $dataDir
# note this doesn't work because of authentication issues, so navigate to here via browser instead
wget $cookieArgs "https://paleon.geography.wisc.edu/lib/exe/fetch.php/data_and_products%3Bpublic_data%3Boh${ohioVersionID}v${ohioVersion}.zip" -O ohiocomp-${ohioVersionID}.${ohioVersion}.zip
wget $cookieArgs "https://paleon.geography.wisc.edu/lib/exe/fetch.php/data_and_products%3Bpublic_data%3B${easternVersionID}centroid_polygonver${easternVersion}.zip" -O easterncomp-${easternVersionID}.${easternVersion}.zip

ohioDir=${ohioVersionID}.${ohioVersion}
mkdir ${ohioDir}
cp oh${ohioVersionID}centroid_polygonsver${ohioVersion}.zip ${ohioDir}
cd ${ohioDir}
unzip oh${ohioVersionID}centroid_polygonsver${ohioVersion}.zip

cd ..
easternDir=${easternVersionID}.${easternVersion}
mkdir ${easternDir}
cp ${easternVersionID}centroid_polygonsver${easternVersion}.zip ${easternDir}
cd ${easternDir}
unzip ${easternVersionID}centroid_polygonsver${easternVersion}.zip

cd $projectDir

########################################################################
# preprocess eastern township data -------------------------------------
########################################################################

./intersect_towns_cells.R  >& log.intersect_towns_cells &
# creates intersection.Rda

./build_eastern.R  >& log.build_eastern_${productVersion}${uniqueRunID} &
# this reads intersection.Rda and creates easternData.Rda


########################################################################
# fit Bayesian composition model to eastern data -----------------------
########################################################################

export OMP_NUM_THREADS=${numCoresToUse} # this seems the sweet spot
./fit_eastern.R >& log.fit_eastern_${productVersion}${uniqueRunID} &
# this creates 'PLScomposition_eastern_${productVersion}.nc'

########################################################################
# subset final output to burned-in samples
########################################################################

# eastern
burnin=25000
echo "burnin=$burnin" > tmp.config
echo "domain=\"eastern\"" >> tmp.config
./remove_burnin.R

# western
burnin=25000
echo "burnin=$burnin" > tmp.config
echo "domain=\"western\"" >> tmp.config
./remove_burnin.R 


########################################################################
# make plots -----------------------
########################################################################


# download us_alb.zip, water tiff, domain tiff from Wiki
# http://144.92.235.115/dokuwiki/doku.php/public_data%3Brasters
wget $cookieArgs "https://paleon.geography.wisc.edu/lib/exe/fetch.php/data_and_products%3Bpublic_data%3Bus_alb.zip" -O tmp_alb.zip
wget $cookieArgs "https://paleon.geography.wisc.edu/lib/exe/fetch.php/data_and_products%3Bpublic_data%3Bwater_pct_albv0.1.tif" -O tmp_water.tif
wget $cookieArgs "https://paleon.geography.wisc.edu/lib/exe/fetch.php/data_and_products%3Bpublic_data%3Bdomain%3Bpaleon_full_alb_v0.1.tif" -O tmp_domain.tif

cp tmp_domain.tif $dataDir/paleonDomain.tif
cp tmp_water.tif $dataDir/water.tif

cd $dataDir
unzip tmp_alb.zip

cd $projectDir
./plot_eastern.R
./plot_western.R

# this creates the mask (paleonMask.nc) for screening out water and non-paleon state cells
./create_mask.R >& log.create_mask &


./plot_full.R