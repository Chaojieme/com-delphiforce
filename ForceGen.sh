#Molecule Interaction Analycal Tool

# Lin Li: 
# This tool is used to calculate electric force from a *.residue force file.
# Usage: ./force_generator.sh list_flie pqr_file force_residue_file output_file

if [ $# -lt 4 ]
then
        echo "can't run, it needs parameters... "
	echo "Usage: ./force_generator.sh list_file pqr_file force_residue_file output_file"
        exit
fi

############# generate force file and pqr file  #############

#for pn in $(cat $1); do awk -vn=$pn '{if(substr($0,23,4)+0==n) print $0}' $2; done > $4_temp_pqr
#for pn in $(cat $1); do awk -vn=$pn '{if(substr($0,7,4)+0==n) print $0}' $3; done > $4_temp_residue

### pqr file ###
while IFS= read -r line
do
  chain=$(echo "$line" | awk '{print substr($0,1,1)}' )
  res=$(echo "$line" | awk '{print substr($0,2,length($0))+0 }' )
  #echo "test:" $line "  chain: " $chain "   res: " $res
  awk -vchain=$chain -vres=$res '{if(substr($0,23,4)+0 == res && substr($0,22,1) == chain) print $0}' $2
done < $1 > $4_temp_pqr

### residue file ###
while IFS= read -r line
do
  chain=$(echo "$line" | awk '{print substr($0,1,1)}' )
  res=$(echo "$line" | awk '{print substr($0,2,length($0))+0 }' )
  #echo "test:" $line "  chain: " $chain "   res: " $res
  awk -vchain=$chain -vres=$res '{if(substr($0,7,4)+0 == res && substr($0,5,1) == chain) print $0}' $3
done < $1 > $4_temp_residue

awk -vx=0 -vy=0 -vz=0 '{x+=substr($0,36,10);y+=substr($0,46,10);z+=substr($0,56,10)}END{print "Total force: ",x,y,z}' $4_temp_residue >> $4_temp_residue


############# generate forece file in tcl format #############
awk -vsumx=0 -vsumy=0 -vsumz=0 -vnum=0 '{if($1=="ATOM" || $1=="HETATM") {sumx=sumx+substr($0,31,8);sumy=sumy+substr($0,39,8);sumz=sumz+substr($0,47,8);num=num+1 } }END{print sumx/num,sumy/num,sumz/num}' $4_temp_pqr|awk '{print "Prob_Center: ", $0}' > $4.temp
#~/LinLi/soft/ll/mybash/center.sh $4_temp_pqr|awk '{print "Prob_Center: ", $0}' > $4.temp

x1=$(grep "Center" $4.temp |awk '{print $2}')
y1=$(grep "Center" $4.temp |awk '{print $3}')
z1=$(grep "Center" $4.temp |awk '{print $4}')


fx=$(grep "Total" $4_temp_residue |awk '{print $3}')
fy=$(grep "Total" $4_temp_residue |awk '{print $4}')
fz=$(grep "Total" $4_temp_residue |awk '{print $5}')


x2=$(echo "$x1 + $fx"|bc -l)
y2=$(echo "$y1 + $fy"|bc -l)
z2=$(echo "$z1 + $fz"|bc -l)

# clean temp files
rm $4.temp $4_temp_pqr $4_temp_residue

#echo "{ $x1 $y1 $z1 } {$fx $fy $fz } {$x2 $y2 $z2}" 
echo "vmd_draw_arrow 0 {$x1 $y1 $z1} {$x2 $y2 $z2} red 30.00" > $4.tcl 

