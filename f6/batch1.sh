#!/bin/bash -l

#PBS -P v14
#PBS -q normal
#PBS -l walltime=01:00:00
#PBS -l mem=80gb
#PBS -l ncpus=12
#PBS -j oe
#PBS -l wd
#PBS -lstorage=gdata/v14+scratch/ux06

#module purge
pwd
ls -al

conda activate pangeo

ci=( nino3 nino4 nino34)
yr=( 2015 2016 2017 2018 )
yr=( 2014 2013 2012 2011 2010 2009  2008  )
yr=( 2020 )
for year in ${yr[*]}
do
for var in ${ci[*]}
do
echo ${year} ${var}
time python climate_indices.py ${year} ${var} 
done
done

mv *.pdf plots
mv *.nc data
