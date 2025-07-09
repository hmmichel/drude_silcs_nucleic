#!/bin/bash
# Copyright SilcsBio LLC, 2016
#
# This file generates FragMaps.
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
This script will submit jobs for generating FragMaps.

Usage: \$SILCSBIODIR/silcs-rna/2b_gen_maps prot=<RNA PDB file>

OPTIONAL parameters:
  margin=<size of margin in the map; default=$margin>
  ref=<reference PDB; default=same PDB file given in prot>
  numsys=<\# of simulations; default=$numsys>
  protmap=<true/false; generate protein side chain maps; default=false>
  halogen=<true/false; default=false>
  spacing=<map spacing; default=1>
  begin=<cycle number to begin for map generation; default=$begin>
  end=<cycle number to end for map generation; default=$end>
@EOF
}

prot=""
sysname=""
workdir=`echo $PWD`
setupdir="1_setup_neutral"
drudedir="${workdir}/2e_run_drude"
qname=""
account=""
email=""
if [ -e ${setupdir}/ligand.txt ]; then lig=`cat ${setupdir}/ligand.txt`; else lig=""; fi
if [ -e ${setupdir}/numsys.txt ]; then numsys=`cat ${setupdir}/numsys.txt`; else numsys=10; fi
margin=10
protmap=false
batch=false
ref=false
halogen=false
standard=true
spacing=1.0
begin=1
end=100


if [[ $# -lt 1 ]]; then
  show_usage
  if ["`ls ${setupdir}/* -d| wc -l`" -gt 0 ]; then
    numsys=`ls ${setupdir}/* -d| wc -l`
    printf "\n *** DETECTED : numsys=$numsys ***\n"
  fi
  printf "\n"
  exit 0
fi

parse_args $@

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

if [[ "$workdir" != "false" ]]; then
  workdir=`realpath $workdir`
  if [[ ! -d "${workdir}" ]]; then
    echo "custom workdir is given, but the directory is not found"
    exit 0
  fi
fi

if [[ "$ref" != "false" ]]; then
  REF_PDB="${ref%.*}"
  if [[ ! -f "${REF_PDB}.pdb" ]]; then
    echo "reference PDB is not found: ${REF_PDB}.pdb"
    exit 0
  fi
else
  REF_PDB=${PROT_PDB}
fi

if [[ "$standard" == "true" ]]; then
  halogen=false
  setupdir="1_setup"
  mapdir="2h_drude_maps"
  probe=standard
fi

#if [[ "$standard" == "neutral" ]]; then
#  halogen=false
#  setupdir="1_setup"
#  drudedir="2e_run_drude"
#  mapdir="2h_drude_maps"
#  probe=standard
#fi

if [[ ! -d "${drudedir}_neutral" ]] || [[ ! -d "${drudedir}_charged" ]]; then
  echo "cannot find Drude dir ${drudedir}_neutral/charged; exit"
  exit 1
fi

echo "########################################################"
echo "prot         = " ${PROT_PDB}
if [ ! -z ${lig} ]; then
echo "lig          = " ${LIGAND}
fi
echo "numsys       = " $numsys
echo "probe        = " $probe
echo "########################################################"
echo -e "\n"

totruns=$numsys

gen_maps () {
  local setupdir=$1
  local drudedir=$2
  local mapdir=$3
  local probe=$4
  # define jobname
  i=$((${#sysname}-4))
  j=`echo $probe | awk '{print toupper($0)}'`
  jobname=DM${p:0:1}-${sysname:$i}
  
  gcx=`grep "NELEMENTS" ${workdir}/silcs_fragmaps_${PROT_PDB}/maps/${PROT_PDB}.benc.gfe.map | awk '{print substr($2,1,2)}' | awk '{print $1/2}'`
  gcy=`grep "NELEMENTS" ${workdir}/silcs_fragmaps_${PROT_PDB}/maps/${PROT_PDB}.benc.gfe.map | awk '{print substr($3,1,2)}' | awk '{print $1/2}'`
  gcz=`grep "NELEMENTS" ${workdir}/silcs_fragmaps_${PROT_PDB}/maps/${PROT_PDB}.benc.gfe.map | awk '{print substr($4,1,2)}' | awk '{print $1/2}'`
  
  #  while true;
  #  do
      sed -e "s/<gcx>/${gcx}/g" \
          -e "s/<gcy>/${gcy}/g" \
          -e "s/<gcz>/${gcz}/g" \
          ${SILCSBIODIR}/drude_silcs/templates/map_prm.tmpl > ${workdir}/map.${sysname}.prm
  #  done
  
  if [[ ! -d ${mapdir} ]]; then mkdir ${mapdir}; fi

  i=1

  while [[ "$i" -le $totruns ]];
  do
    sed -e "s/<run>/${i}/g" \
        -e "s/<jobname>/${jobname}/g" \
        -e "s/<sysname>/${sysname}/g" \
        -e "s~<drudedir>~${drudedir}~g" \
        -e "s/<prot>/${PROT_PDB}/g" \
        -e "s/<lig>/${LIGAND}/g" \
        -e "s/<protmap>/${protmap}/g" \
        -e "s~<SILCSBIODIR>~${SILCSBIODIR}~g" \
        -e "s~<mapdir>~${mapdir}~g" \
        -e "s~<workdir>~${workdir}~g" \
        -e "s/<probe>/${probe}/g" \
	-e "s/<batch>/${batch}/g" \
        -e "s/<ref>/${REF_PDB}/g" \
        -e "s/<margin>/${margin}/g" \
        -e "s/<probe>/${probe}/g" \
        -e "s/<spacing>/${spacing}/g" \
        -e "s/<begin>/${begin}/g" \
        -e "s/<end>/${end}/g" \
        -e "s/<account>/${account}/g" \
        -e "s/<email>/${email}/g" \
        -e "s/<qname>/${qname}/g" \
	${SILCSBIODIR}/drude_silcs/silcs-rna/job_gen_maps_drude.tmpl > ${mapdir}/sub_maps_drude.${i}.sh

  cd ${mapdir}

  if [[ "$batch" == false ]]; then
    output=`sbatch sub_maps_drude.${i}.sh`
    echo $output
    sleep 0.3
  fi

  cd ../

  i=$((i+1))
  done
}

if [[ "${standard}" == "true" ]]; then
  echo -e "Set 1 - Neutral Probes\n-----------------------\n"
  gen_maps "${setupdir}_neutral" "${drudedir}_neutral" "${mapdir}_neutral" "neutral"
  echo -e "\nSet 2 - Charged Probes\n-----------------------\n"
  gen_maps "${setupdir}_charged" "${drudedir}_charged" "${mapdir}_charged" "charged"
fi



if [[ "${halogen}" == "true" ]]; then
  echo "Halogen probes are not supported in SILCS-RNA yet!!!"
  exit 1
  ##gen_maps $setupdir $gcmddir $mapdir "halogen"
fi

