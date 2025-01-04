#!/bin/bash

# Navigate to the directory containing production files
cd david_PROD_CLA/
# cd PROD_SOD/
# cd PROD_POT/

# Alter for each ion 
if [ ! -d david_WHAM_CLA_unclamp ]; then
  mkdir david_WHAM_CLA_unclamp
  echo "Created directory david_WHAM_CLA_unclamp"
else
  echo "Directory david_WHAM_CLA_unclamp already exists"
fi

# Compile files_pullx.dat and files_tpr.dat
# Assumes files are named PROD_{index}.xvg and PROD_{index}.tpr
ls PROD_*_pullx.xvg > files_pullx.dat
ls PROD_*.tpr > files_tpr.dat

# Parameters for WHAM
bins=100
min=-4.8
max=4.8

# Loop through time intervals and perform WHAM analysis
for i in {5..21..2}
do
  echo "Performing WHAM for ${i} ns..."
  j=$((1000*i))  # Convert ns to ps
  gmx wham -ix files_pullx.dat -it files_tpr.dat -bins $bins -min $min -max $max -b 0 -e $j -o david_WHAM_CLA_unclamp/profile_${i}ns.xvg -unit kCal -cycl
done

# Perform additional WHAM analysis for specific intervals
echo "Performing WHAM for 0-7 ns..."
gmx wham -ix files_pullx.dat -it files_tpr.dat -bins $bins -min $min -max $max -b 0 -e 7000 -o david_WHAM_CLA_unclamp/profile_0-7ns.xvg -unit kCal -cycl

echo "Performing WHAM for 7-14 ns..."
gmx wham -ix files_pullx.dat -it files_tpr.dat -bins $bins -min $min -max $max -b 7000 -e 14000 -o david_WHAM_CLA_unclamp/profile_7-14ns.xvg -unit kCal -cycl

echo "Performing WHAM for 14-21 ns..."
gmx wham -ix files_pullx.dat -it files_tpr.dat -bins $bins -min $min -max $max -b 14000 -e 21000 -o david_WHAM_CLA_unclamp/profile_14-21ns.xvg -unit kCal -cycl

echo "WHAM analysis complete. Results saved in david_WHAM_CLA_unclamp/"
