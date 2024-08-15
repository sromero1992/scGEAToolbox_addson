#!/bin/bash

mkdir -p h5_files
cd h5_files

# Initialize metadata file
fname="metadata.csv"
echo "sample_id,grpid" > "$fname"
# Mapped data directory
path0="/mnt/SCDC/Bumblebee/2024_Sato/Project1_2024_SatoAnestacia_CircadianRhythm_MouseBreastCancerTumors_Data/mapped_data"

# Custom pools
pools=("Sato-Pool10" "Sato-Pool11" "Sato-Pool12" "Sato-Pool13" )

echo "Main path : ${path0}"
echo "Pools : ${pools[@]}"

# Expected subdirectory for pulling sample counts
subdir0="outs/per_sample_outs"

for name in "${pools[@]}" ; do 
	echo "Current path $name"
	# Builing whole path
	pathi="${path0}/${name}/${subdir0}"
	# Getting subdirectories in the current variable
	str=`ls ${pathi}`
	subdirs=($str) 
	for namesub in "${subdirs[@]}"; do
		# Writting metadata
		echo "${namesub},${namesub}" >> "$fname"
		echo "      ${namesub}"
		# Copying sample_filtered_features_bc_matrix/h5 to target
		#path1="${pathi}/${namesub}/count/sample_filtered_feature_bc_matrix"
		path1="${pathi}/${namesub}/count/sample_filtered_feature_bc_matrix.h5"
		# Target directory
		target="./${namesub}"
		mkdir -p ${target}
		cp -rv ${path1} ${target}/.
	done
done
cd ../
