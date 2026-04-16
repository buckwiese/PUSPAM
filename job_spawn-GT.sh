#!/bin/bash
#SBATCH -p qcb
#SBATCH -A nmherrer_110

here=`pwd`
jobslist=${here}'/replicates.txt' # can change this to title of your joblist file
there='GT_runs'

#------Create run folder and move unique files into the run folder----------------------------------#
if [[ -d ${there} ]]
then
rm -r ${there}
mkdir ${there}

else
mkdir ${there}

fi

cp ${jobslist} ${there}
cp Global_cell_params.jl ${there}
cp GT_laminarin_params.jl ${there}
cp GT_laminarin_constr.jl ${there}

#-------------- Determine number of jobs you need to run  ----------------#
num_jobs=`wc -l ${jobslist} | awk '{print $1 }'`

echo 'Number of runs: ' ${num_jobs} '...'

#-------------- Loop over all polys (start after headed line)  ----------------------------#
iter=1

while [ ${iter} -lt ${num_jobs} ]
do
   let iter=${iter}+1
   curr_var=`head -${iter} ${jobslist} | tail -1 | cut -f 1`
   runname="run-rep-${curr_var}-GT"

   replicate=`head -${iter} ${jobslist} | tail -1 | cut -f 1`
      
   ####
   
   echo 'Run #'${iter}': creating files for '${runname}'...'
   
   #Create destination file names
   genJL=${there}"/run-rep-${curr_var}-GT"'.jl'
   genSLURM=${there}'/submit-'rep-${curr_var}'.slurm'
   
   #---- copy template and data files to runs folder 
cp -f ${here}'/GT_laminarin_optmsr.jl' ${genJL}
cp -f ${here}'/GT_batch_submit.slurm' ${genSLURM}
   # # # NOTE: This is currently making all files in the /Runs/ directory.
   # # # You can split the files to have them save in separate folders per job as well.
   
# ------- Editing the julia code file ------- #
   joboutname=${runname}'_out'
   joberrname=${runname}'_err'
   
   sed -i -e s@sub_replicate@${replicate}@g    ${genJL}
  
            
# ------- Editing the Slurm files ------- #
joboutname=${there}'/'${runname}'_out'
joberrname=${there}'/'${runname}'_err'

sed -i -e s@sub_out@${joboutname}@g  \
       -e s@sub_err@${joberrname}@g  \
       -e s@sub_julia@${genJL}@g    ${genSLURM}

echo 'Run files successfully modified'
  
done
