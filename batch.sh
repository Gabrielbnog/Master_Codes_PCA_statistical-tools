#rm -f ~/Desktop/fort.dat
#rm -f ~/Desktop/grid_mod.cgns

#user=$(whoami)
#if [ ! -d /tmp/$user/xPOD_* ]; then
#  ninst=0
#else
#  ninst=$(ls -d /tmp/$user/xPOD_* | wc -l)
#  echo $ninst
#fi
#ninst=$(($ninst+1))

#folder="xPOD_${ninst}"
#work_path="/tmp/${user}/${folder}"
#mkdir -p $work_path

#cp -v parameters.in xPOD3D $work_path
#cd $work_path


#for var in "MomentumX" "MomentumY"
#do
#  for zon in $(seq 1 1)
#  do
#    echo $var > working_var.in
#    echo $zon > working_zone.in
#    ./xPOD3D
#  done
#done

for zon in $(seq 1 5)
do
  echo $zon > working_zone.in
  ./xPOD3D
done



