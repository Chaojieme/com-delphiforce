#Molecule Interaction Analycal Tool

# Lin Li:
# This tool is used to calculate electric field, forces, energy for a 2-atom system. The usage is:
# Usage: ./DelPhiForce.sh file1 file2 output_file
# file1 is a reference molecule, file2 is the prob molecule
# After the run, you will get an output_file.residue 

echo -e "###################################################################################################"
echo -e "#### DelPhiForce is used to calculate electric field, forces, energy between two molecules.    ####"
echo -e "#### Please cite this paper:                                                                   ####"
echo -e "#### Lin Li, Arghya Chakravorty, and Emil Alexov. \"DelPhiForce, a tool for electrostatic       ####" 
echo -e "#### force calculations: Applications to macromolecular binding.\"                              ####"
echo -e "#### Journal of Computational Chemistry 38.9 (2017): 584-593.                                  ####"
echo -e "###################################################################################################"

if [ $# -lt 3 ]
then
        echo "can't run, it need parameter... "
	echo "Usage: ./DelPhiForce.sh file1 file2 output_file"
	echo "file1 and file2 are in pqr format. file1 is a reference molecule, file2 is the prob molecule"
        exit
fi

echo -en "DelPhiForce started at "
date

### test if the pdb has chain ID ###
chainid=$(awk -v num=0 '{if((substr($0,1,6)=="ATOM  " ||substr($0,1,6)=="HETATM") && substr($0,22,1)==" ") num++}END{print num}' $2)
if [ $chainid -gt 0 ]
then
	echo -e "\nThe chain ID of the probe molecule is missing. \nPlease add the chain ID for the probe molecule.\n\n"
	exit
fi

### 0. generate $3_param.txt
echo  "generating delphi parameter file."
echo -e "perfil=70.0\nscale=1.0\nin(modpdb4,file="temp_$3_complex",format=pqr)\nindi=2.0\nexdi=80.0\nprbrad=1.4\nsalt=0.15\nbndcon=2\nmaxc=0.01\nlinit=800\nin(frc,file="temp_$3_frc")\nout(frc,file="frc.out")\nsite(a,p,f)\nenergy(s,c)" > $3_param.txt

#echo -e "gsize=165\nscale=2.0\nacenter(0.0,0.0,0.0)\nin(modpdb4,file="temp_$3_complex",format=pqr)\nindi=2.0\nexdi=80.0\nprbrad=1.4\nsalt=0.00\nbndcon=2\nmaxc=0.0001\nlinit=800\nin(frc,file="temp_$3_frc")\nout(frc,file="frc.out")\nsite(a,p,f)\nenergy(s,c)" > $3_param.txt

### 1. Neutralize f2 -> temp_$3_f2
echo "Neutralizing prob molecule"
awk '{if(substr($0,1,6)=="ATOM  " || substr($0,1,6)=="HETATM")printf("%s%8.4f%7.4f\n",substr($0,1,54),0,substr($0,63,69)); else print $0}' $2 > temp_$3_f2

### 2. Generate complex containing f1 and Neutralized f2 (temp_$3_f2) -> temp_$3_complex
echo "Generating complex containing reference molecule and neutralized prob molecule."
cat $1 temp_$3_f2 > temp_$3_complex

### 3. Generate input frc file temp_$3_frc
cp $2 temp_$3_frc

### 4. Run DelPhi to get the output frc file $3
echo "Running DelPhi:"
echo "delphicpp $3_param.txt > $3_delphi.log"
#~/soft/delphicpp_v73 $3_param.txt > $3_delphi.log
./delphicpp $3_param.txt > $3_delphi.log
echo "DelPhi run finished."
### 5. create the outputfile
echo "Creating the outputfile."
awk 'BEGIN{printf("\n\n\n\n\n\n\n\n\n\n\n\n")} {if(substr($0,1,6)=="ATOM  "|| substr($0,1,6)=="HETATM") printf("%10.4f\n",substr($0,55,62))}' $2 > temp_$3_q
paste frc.out temp_$3_q |awk '{printf("%s%10.4f\n",substr($0,1,60),$NF)}' > frc2.out

length=$(wc -l frc2.out |awk '{print $1}')


awk -vl=$length -vfx=0 -vfy=0 -vfz=0 -vg=0 '{ if(NR>12 && NR < l ) {printf("%s%12.4f%10.4f%10.4f%10.4f%10.4f\n",substr($0,1,60),substr($0,61,10),substr($0,21,10)*substr($0,61,10),substr($0,31,10)*substr($0,61,10),substr($0,41,10)*substr($0,61,10),substr($0,51,10)*substr($0,61,10)); g=g+substr($0,21,10)*substr($0,61,10);fx=fx+substr($0,31,10)*substr($0,61,10);fy=fy+substr($0,41,10)*substr($0,61,10);fz=fz+substr($0,51,10)*substr($0,61,10)} else if(NR==12) print "ATOM DESCRIPTOR       GRID PT.    GRID FIELDS: (Ex, Ey, Ez)       CHARGE    ENERGY        FORCES:(Fx, Fy, Fz)"; else if(NR!=l) print $0 }END{printf("Total force: %15.4f%15.4f%15.4f\n",fx,fy,fz);printf("Binding energy:%15.4f\n",g)}' frc2.out > $3.atom



awk -vstr="" -vstr2="" -vl=$length -vq=0 -vg=0 -vfx=0 -vfy=0 -vfz=0 -vflag=0 -vflag2=0 'BEGIN{print "RESDUE ID      NET CHARGE         G        Fx        Fy        Fz"} {if(NR>12 && NR<=l) {flag=substr($0,9,12);str=substr($0,6,15);if(flag2!=flag && flag2!="0") {printf("%s%10.4f%10.4f%10.4f%10.4f%10.4f\n",str2,q,g,fx,fy,fz); q=0;g=0;fx=0;fy=0;fz=0};flag2=flag;str2=str;q=q+substr($0,61,12);g=g+substr($0,73,10);fx=fx+substr($0,83,10);fy=fy+substr($0,93,10);fz=fz+substr($0,103,10)} ; if (NR>=l){print $0} }' $3.atom > $3.residue


#awk -vstr="" -vstr2="" -vl=$length -vq=0 -vg=0 -vfx=0 -vfy=0 -vfz=0 -vflag=0 -vflag2=0 'BEGIN{print "RESDUE ID      NET CHARGE         G        Fx        Fy        Fz"} {if(NR>12 && NR<l) {flag=substr($0,9,12);str=substr($0,6,15);if(flag2!=flag && flag2!="0") {printf("%s%10.4f%10.4f%10.4f%10.4f%10.4f\n",str2,q,g,fx,fy,fz); q=0;g=0;fx=0;fy=0;fz=0};flag2=flag;str2=str;q=q+substr($0,61,12);g=g+substr($0,73,10);fx=fx+substr($0,83,10);fy=fy+substr($0,93,10);fz=fz+substr($0,103,10)} else if (NR>=l){print $0} }' $3.atom > $3.residue


############# generate forece file in tcl format #############
echo "Generate forece file in tcl format."
#~/LinLi/soft/ll/mybash/center.sh $2|awk '{print "Prob_Center: ", $0}'>> $3.atom
awk -vsumx=0 -vsumy=0 -vsumz=0 -vnum=0 '{if($1=="ATOM" || $1=="HETATM") {sumx=sumx+substr($0,31,8);sumy=sumy+substr($0,39,8);sumz=sumz+substr($0,47,8);num=num+1 } }END{print sumx/num,sumy/num,sumz/num}' $2 |awk '{print "Prob_Center: ", $0}'>> $3.atom

x1=$(grep "Center" $3.atom |awk '{print $2}')
y1=$(grep "Center" $3.atom |awk '{print $3}')
z1=$(grep "Center" $3.atom |awk '{print $4}')

fx=$(grep "Total" $3.atom |awk '{print $3}')
fy=$(grep "Total" $3.atom |awk '{print $4}')
fz=$(grep "Total" $3.atom |awk '{print $5}')


x2=$(echo "$x1 + $fx"|bc -l)
y2=$(echo "$y1 + $fy"|bc -l)
z2=$(echo "$z1 + $fz"|bc -l)

#echo "{ $x1 $y1 $z1 } {$fx $fy $fz } {$x2 $y2 $z2}" 
echo "vmd_draw_arrow 0 {$x1 $y1 $z1} {$x2 $y2 $z2} red 30.00" > $3.tcl 


### 6. clean temporary files

rm frc2.out temp_$3_complex  temp_$3_f2  temp_$3_frc  temp_$3_q 

### 7. generate force on each residue
echo "Generating tcl file for each residue."
#for pn in $(awk '{if(substr($0,7,4)+0 >0) print substr($0,5,6)}' $3.residue)
#do 
#  echo $pn > list_$pn
#  ./ForceGen.sh list_$pn $2 $3.residue $3_$pn
#  cat $3_$pn.tcl >> $3_residue.tcl

  #rm list_$pn $3_$pn.tcl
#done

awk '{if(substr($0,7,4)+0 >0) print substr($0,5,6)}' $3.residue > list_temp

while IFS= read -r line
do
  chain=$(echo "$line" | awk '{print substr($0,1,1)}' )
  res=$(echo "$line" | awk '{print substr($0,2,length($0))+0 }' )

  pn=${chain}_${res}
  echo $line > list_${pn}

  ./ForceGen.sh list_$pn $2 $3.residue $3_$pn
  cat $3_$pn.tcl >> $3_residue.tcl

  rm list_$pn $3_$pn.tcl

done < list_temp

rm list_temp

echo -en "Calculation finished. DelPhiForce exit at "
date

