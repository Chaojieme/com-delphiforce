#!/bin/bash

# Script Overview:
# This script automates the process of splitting a PDB file into antigen and antibody chains, 
# converting them into PQR format, running DelPhi calculations, and extracting forces and binding energies.
# It operates in multiple steps and calculates the total force and binding energy for different steps, 
# saving the results into corresponding output files.

# Usage:
# ./your_script.sh <start_step> <end_step> <step_increment>
# 
# - start_step: The initial step for the DelPhi calculations (integer).
# - end_step: The final step for the DelPhi calculations (integer).
# - step_increment: The step increment to move between calculation steps (integer).
#
# Example usage:
# ./your_script.sh 1 10 1
#
# This will run the calculations from step 1 to step 10, incrementing by 1.

# Dependencies:
# The following files and scripts must be available in your working directory:
# 1. masscentpdb_pqr.sh: A script to calculate the center of mass of a PQR file.
# 2. pdb2pqr.py: A Python-based script to convert PDB files to PQR format.
# 3. delphicpp: The executable file for running DelPhi calculations.
#
# Author: Cheng Haojie
#
# Description:
# This script processes PDB files containing antigen and antibody data, 
# and runs DelPhi calculations over a series of steps. The results for 
# total forces and binding energies are computed and saved into output files.
#
# WARNING:
# Make sure to have all the dependencies and the required PDB files in place before running this script.
#
# You can modify the working directory or input PDB file by adjusting the relevant variables in the script.


# Set the working path and parameters
WORKDIR=$(realpath $(pwd))
PDB2PQR=$WORKDIR/pdb2pqr/pdb2pqr.py
MASSCENT=$WORKDIR/masscentpdb_pqr.sh
DELPHICPP=$WORKDIR/delphicpp
input_pdb="HIT2_complex.pdb"
echo "Working directory: $WORKDIR"
start_step=$1
end_step=$2
step_increment=$3



# Traverse all folders to find input_pdb
find . -type f -name "$input_pdb" | while read -r pdb_file; do
    target_dir=$(dirname "$(realpath "$pdb_file")")   # Get the path of pdb_file
    # Prepare output file paths
    result_set="$target_dir/result_set.txt"
    # Prepare paths for antigen and antibody
    antigen_pdb="$target_dir/antigen.pdb"
    antibody_pdb="$target_dir/antibody.pdb"
    antigen_pqr="$target_dir/antigen.pqr"
    antibody_pqr="$target_dir/antibody.pqr"
    
    # 1. Split pdbfile for chainID (assuming chainIDs are PROA for antigen and PROB for antibody)
    awk '{if($12=="PROA")print$0}' "$pdb_file" > "$antigen_pdb"
    awk '{if($12=="PROB")print$0}' "$pdb_file" > "$antibody_pdb"
    
    # 2. Convert pdb file to pqr file
    "$PDB2PQR" --ff=charmm --chain "$antigen_pdb" "$antigen_pqr"
    "$PDB2PQR" --ff=charmm --chain "$antibody_pdb" "$antibody_pqr"

    temp_step=$start_step
    
    cd $target_dir
    # Loop over steps and perform DelPhi calculations
    while [ "$temp_step" -le "$end_step" ]; do
        
        # Create output directory 
        output_dir="Delphiforce_step${temp_step}"
        mkdir -p "$output_dir"
        
        cd "$output_dir"
       
        delphi_parameter="step_${temp_step}_param.txt"
        residue_filename="step_${temp_step}.residue"
        atom_filename="step_${temp_step}.atom"
        antibody_stepN="sep_${temp_step}_antibody.pqr"

        echo "Processing step: $temp_step"
        
        # 3. Execute script separate

        echo "Executing separate script."
        #####finding the center of mass of both structures

        x1=$("$MASSCENT" "$antigen_pqr" |awk '{print $1}')
        y1=$($MASSCENT "$antigen_pqr" |awk '{print $2}')
        z1=$("$MASSCENT" "$antigen_pqr" |awk '{print $3}')

        x2=$("$MASSCENT" "$antibody_pqr" |awk '{print $1}')
        y2=$("$MASSCENT" "$antibody_pqr" |awk '{print $2}')
        z2=$("$MASSCENT" "$antibody_pqr" |awk '{print $3}')

    #####calculating the normalization factor
        n=$(echo $x1 $y1 $z1 $x2 $y2 $z2 | awk '{print (($1-$4)^2+($2-$5)^2+($3-$6)^2)^0.5}')

        dnx=$(echo $x1 $x2 $n|awk '{print ($2-$1)/$3 }' )
        dny=$(echo $y1 $y2 $n|awk '{print ($2-$1)/$3 }' )
        dnz=$(echo $z1 $z2 $n|awk '{print ($2-$1)/$3 }' )

    #####finding final displacement

        dx=$(echo ${temp_step} $dnx |awk '{print $1*$2}')
        dy=$(echo ${temp_step} $dny |awk '{print $1*$2}')
        dz=$(echo ${temp_step} $dnz |awk '{print $1*$2}')
    ##### applying final displacement to the intial position of file2

    awk -v xx=$dx -v yy=$dy -v zz=$dz '{if(substr($0,1,4) == "ATOM" || substr($0,1,6) == "HETATM")  printf("%s%8.3f%8.3f%8.3f%s\n",substr($0,1,30), substr($0,31,8)+xx,substr($0,39,8)+yy,substr($0,47,8)+zz,substr($0,55,100)); else print $0 }' "$antibody_pqr" > "$antibody_stepN"

        # 4. Generate parameter file for DelPhi
        echo "Generating DelPhi parameter file."
        echo -e "perfil=70.0\nscale=2\nin(modpdb4,file=temp_${temp_step}_complex,format=pqr)\nindi=2.0\nexdi=80.0\nprbrad=1.4\nsalt=0.15\nbndcon=2\nmaxc=0.01\nlinit=800\nin(frc,file=${antibody_stepN})\nout(frc,file=frc.out)\nsite(a,p,f)\nenergy(s,c)" > "$delphi_parameter"


        # 5. Neutralize antibody and generate complex
        echo "Neutralizing antibody molecule"
        awk '{if(substr($0,1,6)=="ATOM  " || substr($0,1,6)=="HETATM")printf("%s%8.4f%7.4f\n",substr($0,1,54),0,substr($0,63,69)); else print $0}' "$antibody_stepN" > "temp_${temp_step}_f2"

        echo "Generating complex containing antigen and neutralized antibody."
        cat "$antigen_pqr" "temp_${temp_step}_f2" > "temp_${temp_step}_complex"


        # 6. Run DelPhi to get frc.out
        echo "Running DelPhi:"
        $DELPHICPP "$delphi_parameter" > "step_${temp_step}_delphi.log"
        echo "DelPhi run finished."

        # 7. Process frc.out to generate final result
        echo "Creating output files."
        awk 'BEGIN{printf("\n\n\n\n\n\n\n\n\n\n\n\n")} {if(substr($0,1,6)=="ATOM  " || substr($0,1,6)=="HETATM") printf("%10.4f\n",substr($0,55,62))}' "$antibody_stepN" > "temp_${temp_step}_q"
        paste frc.out "temp_${temp_step}_q" | awk '{printf("%s%10.4f\n",substr($0,1,60),$NF)}' > "frc2.out"

        length=$(wc -l "frc2.out" | awk '{print $1}')

        awk -v l=$length -v fx=0 -v fy=0 -v fz=0 -v g=0 '{ if(NR>12 && NR < l ) {printf("%s%12.4f%10.4f%10.4f%10.4f%10.4f\n",substr($0,1,60),substr($0,61,10),substr($0,21,10)*substr($0,61,10),substr($0,31,10)*substr($0,61,10),substr($0,41,10)*substr($0,61,10),substr($0,51,10)*substr($0,61,10)); g=g+substr($0,21,10)*substr($0,61,10);fx=fx+substr($0,31,10)*substr($0,61,10);fy=fy+substr($0,41,10)*substr($0,61,10);fz=fz+substr($0,51,10)*substr($0,61,10)} else if(NR==12) print "ATOM DESCRIPTOR       GRID PT.    GRID FIELDS: (Ex, Ey, Ez)       CHARGE    ENERGY        FORCES:(Fx, Fy, Fz)"; else if(NR!=l) print $0 }END{printf("Total force: %15.4f%15.4f%15.4f\n",fx,fy,fz);printf("Binding energy:%15.4f\n",g)}' "frc2.out" > "$atom_filename"

        awk -v str="" -v str2="" -v l=$length -v q=0 -v g=0 -v fx=0 -v fy=0 -v fz=0 -v flag=0 -v flag2=0 'BEGIN{print "RESDUE ID      NET CHARGE         G        Fx        Fy        Fz"} {if(NR>12 && NR<=l) {flag=substr($0,9,12);str=substr($0,6,15);if(flag2!=flag && flag2!="0") {printf("%s%10.4f%10.4f%10.4f%10.4f%10.4f\n",str2,q,g,fx,fy,fz); q=0;g=0;fx=0;fy=0;fz=0};flag2=flag;str2=str;q=q+substr($0,61,12);g=g+substr($0,73,10);fx=fx+substr($0,83,10);fy=fy+substr($0,93,10);fz=fz+substr($0,103,10)} ; if (NR>=l){print $0} }' "$atom_filename" > "$residue_filename"

        # 8. Calculate G vector (center of mass difference)
        read ax ay az <<< $($MASSCENT "$antigen_pqr")
        read bx by bz <<< $($MASSCENT "$antibody_stepN")
        Gx=$(echo "$ax - $bx" | bc -l)
        Gy=$(echo "$ay - $by" | bc -l)
        Gz=$(echo "$az - $bz" | bc -l)

        G_norm=$(echo "scale=10; sqrt($Gx^2 + $Gy^2 + $Gz^2)" | bc -l)
        echo "G_vector: $Gx $Gy $Gz"

        # 9. Extract Total Force and Binding Energy
        read Fx Fy Fz <<< $(awk '/Total force:/ {print $(NF-2), $(NF-1), $NF}' "$residue_filename")
        BE=$(grep -i "Binding energy" "$residue_filename" | grep -oP '[-+]?[0-9]*\.?[0-9]+' | head -n 1)

        # Check if force or binding energy is missing
        if [[ -z $Fx || -z $Fy || -z $Fz || -z $BE ]]; then
            echo "❌ ERROR: Can't get force or binding energy."
            exit 1
        fi

        # Calculate force projection along G vector
        dot_product=$(echo "$Fx * $Gx + $Fy * $Gy + $Fz * $Gz" | bc -l)
        projection=$(echo "$dot_product / $G_norm" | bc -l)

        # Write results to the output file
        echo "Step ${temp_step}" > "result_step${temp_step}.txt"
        echo "G: $Gx $Gy $Gz" >> "result_step${temp_step}.txt"
        echo "Total Force: $Fx $Fy $Fz" >> "result_step${temp_step}.txt"
        echo "Binding Energy: $BE" >> "result_step${temp_step}.txt"
        printf "dot (Force, G): %.10f\n" "$projection" >> "result_step${temp_step}.txt"
        echo "\n" >> "result_step${temp_step}.txt"
        
        cat "result_step${temp_step}.txt" >> "$result_set"

        # Move to the next step
        temp_step=$((temp_step + step_increment))
        cd ../
        

    done
    cd "$WORKDIR"
    
done

echo "All processing complete."
