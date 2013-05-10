
#!/bin/bash

flag=0;
while [  $flag -eq 0 ]; do
  a=`bjobs -J forWithBalance_VariableReadLength_c200_full_21_40_12 | grep "is not found" -`;
  if [ -f fileName ]
 then
  else
    
  fi
 if [[ "$a"==*"is not found"* ]]; then  
     echo "not submitteed";
     sleep 20;
  else 
    echo "was submitted";
    flag=1;
  fi
done


/seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/BAC/dataForSim/VariableReadLength_c200_full_21_40/forWithBalance_VariableReadLength_c200_full_21_40_12



% check if submitted
sleep 360
flag=0;
while [  $flag -eq 0 ]; do
      b=`bjobs -J forWithBalance_VariableReadLength_c200_full_21_40_34 | grep priority -`;
      a=`bjobs -J forWithBalance_VariableReadLength_c200_full_21_40_34 | grep "is not found" -`;
      if [[ "$a"==*"is not found"* ]]; then 
        echo "not submitteed";
        sleep 20;
      else
        echo "was submitted";
        flag=1;
      fi
   
    if [[ "$b"==*"priority"* ]]; then 
      echo "was submitted";
      flag=1;
    fi
done

% now - sure it was submitted

% if finished
flag=0;
while [ $flag -eq 0 ]; do
  if [ -f /seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/BAC/dataForSim/VariableReadLength_c200_full_21_40/VariableReadLength_c200_full_21_40_bac_dist_flag_0_Nbac_in_mixture_200_npower_10_readlen_75_numIter_14_Nreads_1000000_Noise_1_Correction_1_forBalance.o ]
  then
   echo "already finished"
   flag=1;
  else
    echo "priority still runnning"
    sleep 60;
  fi
done
