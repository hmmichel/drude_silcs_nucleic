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
workdir=`echo $PWD`
setupdir="${workdir}/2a_run_gcmd"
ions=""
nproc=""
qname=""
account=""
email=""
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

#echo ${workdir}
#echo ${drudedir}

if [[ $# -lt 1 ]]; then
  show_usage
  if ["`ls ${setupdir}_neutral/* -d| wc -l`" -gt 0 ]; then
    numsys=`ls ${setupdir}/* -d| wc -l`
    printf "\n *** DETECTED : numsys=$numsys ***\n"
  fi
  printf "\n"
  exit 0
fi

parse_args $@

#echo $@
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

PROT_PDB=$(basename $prot)
PROT_PDB="${PROT_PDB%.*}"

echo "########################################################"
echo "prot         = " ${PROT_PDB}
echo "sysname     = " $sysname
echo "numsys       = " $numsys
echo "nproc        = " $nproc
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
        if [[ ! -f ${drudedir}/${i}/${j}/build/${sysname}.silcs.min.drude.pdb ]]; then
          echo "NO DRUDE MIN PDB FOUND FOR THIS RUN"
          exit 1
        fi
      else
        nligatoms=0
        if [[ ! -f ${drudedir}/${i}/${j}/build/${sysname}.silcs.min.drude.pdb ]]; then
        echo "NO DRUDE MIN PDB FOUND FOR THIS RUN"
        exit 1
        fi
      fi
# check for drude directory for appropriate setup and interval
      if [[ ! -d ${drudedir}/${i}/${j}/equil ]]; then mkdir -p ${drudedir}/${i}/${j}/equil; fi

      builddir=${drudedir}/${i}/${j}/build
      equildir=${drudedir}/${i}/${j}/equil

      cp ${builddir}/${sysname}.silcs.min.drude.crd ${equildir}
      cp ${builddir}/${sysname}.drude.silcs.xplor.psf ${equildir}
    
      cd ${equildir}
    
      # create boxsize file using system_properties.str in prep/
      python ${SILCSBIODIR}/drude_silcs/scripts/create_box_rst.py ${builddir}/system_properties.str ${sysname}
      
      # create restraint file
      if [[ ions != "" ]]; then
        python ${SILCSBIODIR}/drude_silcs/scripts/make_restraint.py -crd ${sysname}.silcs.min.drude.crd -addsegments HETA
      else
        python ${SILCSBIODIR}/drude_silcs/scripts/make_restraint.py -crd ${sysname}.silcs.min.drude.crd
      fi

    done
    
    cd ${drudedir}/${i}
    
    sed -e "s~<toppardir>~${toppardir}~g" \
        ${SILCSBIODIR}/drude_silcs/silcs-rna/toppar_drude_silcs_openmm.tmpl > ${drudedir}/toppar_drude_silcs_openmm.str
          
    # create equilibration submission file
    sed -e "s/<sysname>/${sysname}/g" \
        -e "s~<drudedir>~${drudedir}~g" \
        -e "s/<run>/${i}/g" \
        -e "s/<int>/${j}/g" \
        -e "s/<temp>/${temp}/g" \
        -e "s/<jobname>/${jobname}/g" \
        -e "s/<nproc>/${nproc}/g" \
        -e "s~<SILCSBIODIR>~${SILCSBIODIR}~g" \
        -e "s/<batch>/${batch}/g" \
        -e "s/<gpu>/${gpu}/g" \
	-e "s/<account>/${account}/g" \
  	-e "s/<email>/${email}/g" \
	-e "s/<qname>/${qname}/g" \
        ${SILCSBIODIR}/drude_silcs/templates/job_drude_silcs_equil.tmpl > ${drudedir}/${i}/sub_drude_silcs_equil.sh

    
    nbadlinks=`find . -xtype l | wc -l`

    if [[ "$nbadlinks" -eq 0 ]]; then
      if [[ "$batch" == false ]]; then
        output=`sbatch sub_drude_silcs_equil.sh`
        echo $output
        taskid=`echo $output | awk '{print $3}' | tr -d '<>'`
        echo "$taskid" > .progress
        echo "0 $end" >> .progress
      else
        if [[ ${i} == ${totruns} ]]; then
          echo "> Your job submission scripts are ready!"
        fi
      fi
    else
      echo "${drudedir}/${i} has bad symlinks; not submitting the job;"
    fi

    i=$((i+1))
  done
}
   
if [[ "${standard}" == "true" ]]; then
  echo -e "Set 1 - Neutral Probes\n-----------------------\n"
  run_drude "2a_run_gcmd_neutral" "${drudedir}_neutral" "neutral"
  echo -e "\nSet 2 - Charged Probes\n-----------------------\n"
  run_drude "2a_run_gcmd_charged" "${drudedir}_charged" "charged"
  echo -e ""
fi 
