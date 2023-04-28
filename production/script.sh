#for era in "C" "Dv1" "Dv2" "E" "F" "G"
for era in "G"
do
  for i in 0 1 2 3 4 5 6 7
  do
    crab status -d crab_MuMuTrigger_2022${era}_DMLM${i}_19feb23
 #  crab resubmit -d crab_Bc_2022F_DMLM${i}_13feb23 --maxjobruntime=2000
    echo
    echo
    echo
  done
done
