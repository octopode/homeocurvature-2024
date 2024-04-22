# will just use starting params hardcoded in script

for DAT in ../01-rawdata/saxs/JWL225*1.dat; do BN=$( basename $DAT ); ./lamellar_itfit.py $DAT > ../02-tidydata/mcg_params/${BN/.dat/.tsv}; done
for DAT in ../01-rawdata/saxs/JWL224*1.dat; do BN=$( basename $DAT ); ./lamellar_itfit.py $DAT > ../02-tidydata/mcg_params/${BN/.dat/.tsv}; done