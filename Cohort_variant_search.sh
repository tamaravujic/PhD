#!/bin/bash

#Path to the master file containing all the file paths 
master_file="file_paths.txt"

#variant or region you are looking for 
Variant= ""

#loop through each file path in the masterfile 
While IFS= read -r file path; do
	echo "checking $filepath"
     #use zgrep to search within the compressed VCF file 
       if zgrep -q "$variant" "$filepath"; then 
	echo "Variant found in $filepath"
      else 
	echo "Variant not found in $filepath"
      fi
Done < "$master_file"
