#!/bin/bash
# controller script for fitting the PalEON settlement era composition model
# both the western PLS data and the eastern township data
# note that these steps are intended for use on UNIX-like machines and will need to be modified for Windows (and possibly for Mac OS X)

# this code is being run under R 3.1.2 and with package versioning controlled by packrat
# restore any packages that are not installed on the system
Rscript -e "require(packrat); packrat::restore()"

# modify the contents of the config file to reflect the data versions to be used, relevant directories, and parameters of the MCMC
source config

export OMP_NUM_THREADS=1

########################################################################
# create directories    ------------------------------------------------
########################################################################

if [ ! -e $plotDir ]; then
    mkdir $plotDir
fi
if [ ! -e $outputDir ]; then
    mkdir $outputDir
fi
if [ ! -e $dataDir ]; then
    mkdir $dataDir
fi
if [ ! -e $tmpDir ]; then
    mkdir $tmpDir
fi

########################################################################
# setup for down/uploading from Wiki   ---------------------------------
########################################################################

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
if [ ! -e level3s_v${productVersion}.csv ]; then
    cd $projectDir
    wget $cookieArgs "https://paleon.geography.wisc.edu/lib/exe/fetch.php/data_and_products%3Blevel3s_v${productVersion}.csv" -O level3s_v${productVersion}.csv
fi
if [ -e level3s.csv ]; then
    rm -f level3s.csv
fi

ln -s level3s_v${productVersion}.csv level3s.csv

########################################################################
# download western data ------------------------------------------------
########################################################################

cd $dataDir

if [ ! -e western-${westernVersion}.csv ]; then
    wget $cookieArgs "https://paleon.geography.wisc.edu/lib/exe/fetch.php/data_and_products%3Bwestern_comp_v${westernVersion}.csv" -O western-${westernVersion}.csv
fi

if [ -e western.csv ]; then
    rm -f western.csv
fi

ln -s western-${westernVersion}.csv western.csv

cd $projectDir

########################################################################
# preprocess western data ----------------------------------------------
########################################################################

# cp config_{version}-{runID} config

if [ ! -e data_western_${runID}.Rda ]; then
    ./build_western.R >& log.build_western_${runID} &
# this creates 'data_western_${runID}.Rda'
fi

########################################################################
# fit Bayesian composition model to western data -----------------------
########################################################################

./fit_western.R >& log.fit_western_${runID} &
# this creates 'PLScomposition_western_${runID}_full.nc'

########################################################################
# download eastern township data -------------------------------------
########################################################################


cd $dataDir

if [ ! -e ohio ]; then
    mkdir ohio
fi
if [ ! -e eastern ]; then
    mkdir eastern
fi

if [ ! -e oh${ohioVersionID}centroid_polygonsver${ohioVersion}.zip ]; then
    wget $cookieArgs "https://paleon.geography.wisc.edu/lib/exe/fetch.php/data_and_products%3Boh${ohioVersionID}centroid_polygonsver${ohioVersion}.zip" -O ohio/ohio.zip
    cd ohio
    unzip ohio.zip
    cd ..
fi

if [ ! -e ${easternVersionID}polygonsver${easternVersion}.zip ]; then
    wget $cookieArgs "https://paleon.geography.wisc.edu/lib/exe/fetch.php/data_and_products%3B${easternVersionID}polygonsver${easternVersion}.zip" -O eastern/eastern.zip
    cd eastern
    unzip eastern.zip
fi

cd $projectDir

########################################################################
# preprocess eastern township data -------------------------------------
########################################################################

if [ ! -e data_eastern_${runID}.Rda ]; then
    ./intersect_towns_cells.R  >& log.intersect_towns_cells 
# creates intersection_${runID}.Rda

    ./build_eastern.R  >& log.build_eastern_${runID} &
# this reads intersection_eastern_${runID}.Rda and creates 'data_eastern_${runID}.Rda
fi

########################################################################
# fit Bayesian composition model to eastern data -----------------------
########################################################################

export OMP_NUM_THREADS=${numCoresToUse} # this seems the sweet spot
./fit_eastern.R >& log.fit_eastern_${runID} &
# this creates 'PLScomposition_eastern_${runID}_full.nc'

########################################################################
# subset final output to burned-in samples
########################################################################

# eastern
burnin=25000
domain=eastern
./remove_burnin.R $burnin $domain 
# this creates 'PLScomposition_eastern_${runID}.nc'

# western
burnin=25000
domain=western
./remove_burnin.R $burnin $domain
# this creates 'PLScomposition_western_${runID}.nc'

########################################################################
# do cross-validation -----------------------
########################################################################

if [ $cv = "TRUE" ]; then
    ./calc_cv_western.R >& log.calc_cv_western_${runID} &
fi


########################################################################
# make plots -----------------------
########################################################################


# download us_alb.zip, water tiff, domain tiff from Wiki
# http://144.92.235.115/dokuwiki/doku.php/public_data%3Brasters
wget $cookieArgs "https://paleon.geography.wisc.edu/lib/exe/fetch.php/data_and_products%3Bus_alb.zip" -O tmp_alb.zip
wget $cookieArgs "https://paleon.geography.wisc.edu/lib/exe/fetch.php/data_and_products%3Bwater_pct_albv0.1.tif" -O tmp_water.tif
wget $cookieArgs "https://paleon.geography.wisc.edu/lib/exe/fetch.php/data_and_products%3Bdomain%3Bpaleon_full_alb_v0.1.tif" -O tmp_domain.tif

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