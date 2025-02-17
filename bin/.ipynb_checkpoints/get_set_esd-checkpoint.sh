#!/bin/bash
# 2021-06-24
# M Lazecky 2020, 2021
# Nergizci for ESD 2023
# Script for calculating Earth tides

# Check if the required arguments are provided
if [ -z "$2" ]; then
    echo "Usage: $0 in_df out_SET"
    echo "e.g.: $0 frame.csv tides.csv"
    exit 1
fi

# Store input and output file paths and frame name
in_esds="$1"
out_SET="$2"
# Extract frame name without the extension
frame=$(basename "$in_esds" | cut -d '.' -f1)

homedir=$HOME
batchdir=$BATCH_CACHE_DIR

# Define input and output directory paths
input_dir="$batchdir/daz_esd"
output_dir="$batchdir/daz_esd"

# Check if the input file exists and the output file doesn't exist
if [ ! -f "$input_dir/$in_esds" ]; then
    echo "ERROR: Please check if the following files exist:"
    echo "in_esds: $input_dir/$in_esds"
    echo "Note that this file must NOT exist:"
    echo "out_SET: $output_dir/$out_SET"
    exit
fi
echo $output_dir/$out_SET
echo "epoch,lat,lon,E,N,U,daz_tide " > "$output_dir/$out_SET"
echo "Processing... It can take time"

for aline in `cat $input_dir/$in_esds | tail -n+2 `; do
    # Extract values from the current line
    heading=`echo $aline | cut -d ',' -f3`
    epochdate=`echo $aline | cut -d ',' -f1`
    masterdate=`echo $aline | cut -d ',' -f2`
    centertime=`echo $aline | cut -d ',' -f4`
    lat=`echo $aline | cut -d ',' -f6`
    lon=`echo $aline | cut -d ',' -f5`
    masterdt="$masterdate""T""$centertime"
    epochdt="$epochdate""T""$centertime"
    #echo "$lat,$lon,$masterdt,$heading"
    #echo "$lat,$lon,$epochdt,$heading"
    
    # Calculate earth tide values for master date and current epoch date
    mtide=`gmt earthtide -L$lon/$lat -T$masterdt 2>/dev/null | sed 's/\t/,/g'`
    NM=`echo $mtide | cut -d ',' -f2`
    EM=`echo $mtide | cut -d ',' -f3`
    VM=`echo $mtide | cut -d ',' -f4`
    # Convert scientific notation to decimal notation
    NM=$(printf "%.12f" "$NM")
    EM=$(printf "%.12f" "$EM")
    VM=$(printf "%.12f" "$VM")
    
    etide=`gmt earthtide -L$lon/$lat -T$epochdt 2>/dev/null | sed 's/\t/,/g'`
    NE=`echo $etide | cut -d ',' -f2`
    EE=`echo $etide | cut -d ',' -f3`
    VE=`echo $etide | cut -d ',' -f4`
    # Convert scientific notation to decimal notation
    NE=$(printf "%.12f" "$NE")
    EE=$(printf "%.12f" "$EE")
    VE=$(printf "%.12f" "$VE")
    #echo "NE: $NE, EE: $EE, VE: $VE"
    U=`echo '('$VE')-('$VM')' | bc -l`
    #echo "$U"
    E=`echo '('$EE')-('$EM')' | bc -l`
    #echo "$E"
    N=`echo '('$NE')-('$NM')' | bc -l`
    #echo "$N"
    #echo "U: $U, E: $E, N: $N"
    daz_tide=$(echo "($E * s($heading * 3.14159/180) + $N * c($heading * 3.14159/180))*1000" | bc -l)
     
    # Append processed data to the output file
    #echo "$epochdate, $E, $N, $U, $daz_tide"
    echo "$epochdate,$lat, $lon, $E, $N, $U, $daz_tide" >> "$output_dir/$out_SET"
done
echo "SET calculcated..."