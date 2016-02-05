#!/bin/bash
# controller script for fitting the PalEON settlement era composition model
# both the western PLS data and the eastern township data
# note that these steps are intended for use on UNIX-like machines and will need to be modified for Windows (and possibly for Mac OS X)

# this is not intended to be run as a full script as later components depend on earlier ones having finished (and the fitting steps take many days); rather it is intended to allow one to run all of the steps of the model fitting/analysis

# this code was originally being run under R 3.1.2 (now R 3.2.0) and with package versioning controlled by packrat; note that use of BCE-2015-fall is very likely to work fine (it provides R 3.2.2)
# restore any packages that are not installed on the system
Rscript -e "require(packrat); packrat::restore()"

# modify the contents of the config file to reflect the data versions to be used, relevant directories, and parameters of the MCMC
# in general, it's good to create a version of config, say config_0.4-0, specific to each run and then copy that file to 'config'

\cp config_0.4-0 config

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
if [ ! -e level3s_v${taxonTranslationVersion}.csv ]; then
    cd $projectDir
    wget $cookieArgs "https://paleon.geography.wisc.edu/lib/exe/fetch.php/data_and_products%3Blevel3s_v${taxonTranslationVersion}.csv" -O level3s_v${taxonTranslationVersion}.csv
fi
if [ -e level3s.csv ]; then
    rm -f level3s.csv
fi

ln -s level3s_v${taxonTranslationVersion}.csv level3s.csv

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

if [ ! -e $dataDir/data_western_${runID}.Rda ]; then
    ./build_western.R >& log.build_western_${runID} &
# this creates 'data_western_${runID}.Rda'
fi

########################################################################
# fit Bayesian composition model to western data -----------------------
########################################################################

./fit_western.R >& log.fit_western_${runID} &
# this creates 'PLScomposition_western_${runID}_full.nc'
# note that this netCDF has y-values from N to S (contradicting the 
# dim info for the y dim)

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

if [ ! -e oh${ohioVersionID}polygons_v${ohioVersion}.zip ]; then
    wget $cookieArgs "https://paleon.geography.wisc.edu/lib/exe/fetch.php/data_and_products%3Boh${ohioVersionID}polygons_v${ohioVersion}.zip" -O ohio/ohio.zip
    cd ohio
    unzip ohio.zip
    cd ..
fi

if [ ! -e ${easternVersionID}polygonsver${easternVersion}.zip ]; then
    wget $cookieArgs "https://paleon.geography.wisc.edu/lib/exe/fetch.php/data_and_products%3B${easternVersionID}polygons_v${easternVersion}.zip" -O eastern/eastern.zip
    cd eastern
    unzip eastern.zip
fi

cd $projectDir

########################################################################
# preprocess eastern township data -------------------------------------
########################################################################

if [ ! -e $dataDir/data_eastern_${runID}.Rda ]; then
    ./intersect_towns_cells.R  >& log.intersect_towns_cells 
# creates intersection_${runID}.Rda

    ./build_eastern.R  >& log.build_eastern_${runID} &
# this reads intersection_eastern_${runID}.Rda and creates 'data_eastern_${runID}.Rda
fi

########################################################################
# fit Bayesian composition model to eastern data -----------------------
########################################################################

./fit_eastern.R >& log.fit_eastern_${runID} &
# this creates 'PLScomposition_eastern_${runID}_full.nc'
# note that this netCDF has y-values from N to S (contradicting the 
# dim info for the y dim)

########################################################################
# subset final output to burned-in samples
########################################################################

# this also flips the N->S orientation of netCDF files to their correct
# S->N to match the y dimension as given in the netCDF

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


cp $outputDir/PLScomposition_western_${runID}.nc /server/web/share/paciorek/paleon/composition_midwest_v${productVersion}.nc
cp $outputDir/PLScomposition_eastern_${runID}.nc /server/web/share/paciorek/paleon/composition_east_v${productVersion}.nc

# create merged full domain version

burnin=25000
./stitch_domains.R $burnin
# this gives errors for Atl Wh Ced and Chestnut in west - I presume this was happening before and is why the ncvar_get is wrapped in a try()

########################################################################
# make netCDFs with summary stats (posterior mean, sd, etc.) -----------
########################################################################

./summarize_posterior.R 

########################################################################
# do cross-validation -----------------------
########################################################################

if [ $cv = "TRUE" ]; then
# machs="scf-sm00 scf-sm01 scf-sm02 scf-sm03 scf-sm10 scf-sm11 scf-sm12 scf-sm13"
#for mach in $machs; do  ssh $mach mkdir -p /var/tmp/paciorek/paleon/comp/{data,output}; rsync -av /var/tmp/paciorek/paleon/comp/data/data_western* paciorek@$mach:/var/tmp/paciorek/paleon/comp/data/; done
    for (( i = 101; i <= 106; i++ )); do
        \cp config_${productVersion}-${i} config
        source config
        if [ ! -e $dataDir/data_western_${runID}.Rda ]; then
            ./build_western.R >& log.build_western_${runID} 
            # this creates 'data_western_${runID}.Rda'
        fi
        ./fit_western.R >& log.fit_western_${runID}
        ./calc_cv_western.R >& log.calc_cv_western_${runID} 
    done
fi


########################################################################
# make plots -----------------------
########################################################################


# download us_alb.zip, water tiff, domain tiff from Wiki
# http://144.92.235.115/dokuwiki/doku.php/public_data%3Brasters
wget $cookieArgs "https://paleon.geography.wisc.edu/lib/exe/fetch.php/data_and_products%3Bus_alb.zip" -O $dataDir/tmp_alb.zip
wget $cookieArgs "https://paleon.geography.wisc.edu/lib/exe/fetch.php/data_and_products%3Bwater_pct_alb_v0.1.tif" -O $dataDir/water.tif 
wget $cookieArgs "https://paleon.geography.wisc.edu/lib/exe/fetch.php/data_and_products%3Bdomain%3Bpaleon_full_alb_v0.1.tif" -O $dataDir/paleonDomain.tif

cd $dataDir
unzip tmp_alb.zip

# this is deprecated as we concentrate on plots of full domain below
cd $projectDir
./plot_eastern.R
./plot_western.R

# this creates the mask (paleonMask.nc) for screening out water and non-paleon state cells
./create_mask.R >& log.create_mask &

./plot_full.R
