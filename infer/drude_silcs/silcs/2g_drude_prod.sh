#!/bin/bash
# Copyright SilcsBio LLC, 2016
#
# This file setup SILCS simulations for independent simulation boxes,
# and submit jobs to the queue.
#

load_core_library()
{
  if [[ "$SILCSBIO_LIB_DIR" == "" ]]; then
    SILCSBIO_LIB_DIR="$SILCSBIODIR/lib"
  fi
  if [[ "$SILCSBIO_TEMPLATE_DIR" == "" ]]; then
    SILCSBIO_TEMPLATE_DIR="$SILCSBIODIR/templates"
  fi

  core="$SILCSBIO_LIB_DIR/core.sh"
  enc_core="$SILCSBIO_LIB_DIR/core.enc"
  CORE_LIB=0
  if [[ ! -f "$core" ]]; then
    if [[ -f "$enc_core" ]]; then
      source <(unzip -P "4DNE=ihz;rKZ" -p "${enc_core}" 2>/dev/null)
    else
      echo "Error: library file not found"
      exit 1
    fi
  else
    source "${core}"
  fi
  if [[ "$CORE_LIB" -ne 1 ]]; then
    echo "Error: library is not correctly read"
    exit 1
  fi
}

load_core_library


show_usage()
{
  cat <<@EOF
This script will submit GCMC/MD jobs to the queue system.

Usage: \$SILCSBIODIR/silcs-rna/2a_run_gcmd prot=<RNA PDB file>

Optional Parameters:
  temp=<simulation temperature; default=$temp>
  halogen=<use halogen probes true/false; default=$halogen>

Additional parameters:
  nproc=<\# of cores per simulation; default nproc=$nproc>
  batch=<only generate inputs and not submit jobs true/false; default=$batch>
@EOF
}
prot=""
sysname=""
toppardir=""
fcsilcs=""
workdir=`echo $PWD`
setupdir="${workdir}/2a_run_gcmd"
nproc=""
qname=""
account=""
email=""
customres=""
lig=""
numsys=10
end=100
batch=false
halogen=false
memb=false
standard=true
gpu=false
cleanup=true
protmap=true
spacing=1.0
drudedir="${workdir}/2e_run_drude"
temp=298



if [[ $# -lt 1 ]]; then
    show_usage
    if ["`ls ${setupdir}/* -d| wc -l`" -gt 0 ]; then
        numsys=`ls ${setupdir}/* -d| wc -l`
        printf "\n *** DETECTED : numsys=$numsys ***\n"
    fi
    printf "\n"
    exit 0
fi

parse_args "$@"

# required arguments
for arg in prot
do
  if [[ "${!arg}" == "" ]]; then
    echo "argument \"${arg}\" is required."
    exit 1
  fi
done

for arg in sysname
do
  if [[ "${!arg}" == "" ]]; then
    echo "argument \"${arg}\" is required."
    exit 1
  fi
done

for arg in fcsilcs
do
  if [[ "${!arg}" == "" ]]; then
    fcsilcs=50.208
    exit 1
  fi
done

#echo "\"${customres}\""

PROT_PDB=$(basename $prot)
PROT_PDB="${PROT_PDB%.*}"

echo "########################################################"
echo "prot         = " ${PROT_PDB}
echo "sysname      = " $sysname
echo "numsys       = " $numsys
echo "customres    = " $customres
echo "probe        = " $probe
echo "temp         = " $temp
echo "########################################################"
echo -e "\n"

totruns=$numsys

run_drude () 
{
  local setupdir=$1
  local drudedir=$2
  local probe=$3

  # define jobname
  p=`echo $probe | awk '{print toupper($0)}'`
  jobname=DS${p:0:1}-${PROT_PDB:0:4}

  if [[ ! -d ${drudedir} ]]; then mkdir ${drudedir}; fi

# first loop for drude preparation
  i=1
  while [[ "$i" -le "$totruns" ]];do
    for j in {10..100..10}; do
     if [ ! -z ${lig} ]; then
        nligatoms=`grep -e ATOM ${setupdir}/${LIGAND}_gmx.pdb | wc -l`
        if [[ ! -f ${drudedir}/${i}/${j}/equil/${sysname}.${i}.drude.silcs.${j}.npt.pdb ]]; then
          echo "NO DRUDE EQUIL PDB FOUND FOR THIS RUN"
          exit 1
        fi
      else
        nligatoms=0
        if [[ ! -f ${drudedir}/${i}/${j}/equil/${sysname}.${i}.drude.silcs.${j}.npt.pdb ]]; then
        echo "${i} ${j} NO DRUDE EQUIL PDB FOUND FOR THIS RUN"
        continue 1
        fi
      fi
      #  check for drude directory for appropriate setup and interval
      if [[ ! -d ${drudedir}/${i}/${j}/prod ]]; then mkdir -p ${drudedir}/${i}/${j}/prod ; fi

      equildir=${drudedir}/${i}/${j}/equil
      proddir=${drudedir}/${i}/${j}/prod

      cp ${equildir}/${sysname}.${i}.drude.silcs.${j}.npt.pdb ${proddir}
      cp ${equildir}/${sysname}.drude.silcs.xplor.psf ${proddir}
        
    done
  i=$((i+1))

  done
  
  cd ${drudedir}

  # activate conda env to use MDAnalysis
  #module reset
  #module load Anaconda3/2020.11

  module load Anaconda3/2020.11

  source activate /projects/lemkul_lab/software/infer/openmm/7.7.0/conda_envs/openmm7.7.0 
#  which python

  # create restraint file
  if [[ ${customres} == "" ]]; then
    python ${SILCSBIODIR}/drude_silcs/scripts/make_restraint_silcs.py -pdb ${drudedir}/1/10/equil/${sysname}.1.drude.silcs.10.npt.pdb 
  else 
    echo "Running SILCS with custom restraints on "\"${customres}\"""
    echo "python make_restraint_silcs.py -pdb ${sysname}.drude.silcs.npt.pdb -use_custom "\"${customres}\"""
    python ${SILCSBIODIR}/drude_silcs/scripts/make_restraint_silcs.py -pdb ${drudedir}/1/10/equil/${sysname}.1.drude.silcs.10.npt.pdb -use_custom "\"${customres}\""
  fi

  # create production submission files
  for i in {1..10}; do
    sed -e "s/<sysname>/${sysname}/g" \
        -e "s~<drudedir>~${drudedir}~g" \
        -e "s/<run>/${i}/g" \
        -e "s/<fcsilcs>/${fcsilcs}/g" \
        -e "s/<jobname>/${jobname}/g" \
        -e "s/<nproc>/${nproc}/g" \
        -e "s~<SILCSBIODIR>~${SILCSBIODIR}~g" \
        -e "s/<batch>/${batch}/g" \
        -e "s/<gpu>/${gpu}/g" \
        -e "s/<account>/${account}/g" \
 	      -e "s/<email>/${email}/g" \
        -e "s/<qname>/${qname}/g" \
        ${SILCSBIODIR}/drude_silcs/templates/job_drude_silcs_prod.tmpl > ${drudedir}/sub_drude_silcs_prod_${i}.sh
    nbadlinks=`find . -xtype l | wc -l`

    if [[ "$nbadlinks" -eq 0 ]]; then
      if [[ "$batch" == false ]]; then
        output1=`sbatch sub_drude_silcs_prod_${i}.sh`
        echo ${output1}
      taskid=`echo $output | awk '{print $3}' | tr -d '<>'`
      echo "$taskid" > .progress
      echo "0 $end" >> .progress
      fi
    else
      echo "${drudedir}/${i} has bad symlinks; not submitting the job;"
    fi
  done 
}
   
if [[ "${standard}" == "true" ]]; then
  echo -e "Running Drude SILCS - Prod \n-----------------------\n"
  run_drude "${setupdir}" "${drudedir}" "standard"
  echo -e ""
fi
   
    







