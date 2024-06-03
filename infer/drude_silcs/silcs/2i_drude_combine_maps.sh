#!/bin/bash
# Copyright SilcsBio LLC, 2016-2018
#
# This file setup collates SILCS FragMaps.
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
This script will collate occupancy maps (computed using 2b_gen_maps) to create GFE FragMap.

Usage: \$SILCSBIODIR/silcs-rna/2c_fragmap prot=<RNA PDB file>

OPTIONAL parameters:
  numsys=<\# of simulations; default=$numsys>
  protmap=<true/false; generate protein side chain maps; default=false>
  skipoc=<true/false; skip overlap coefficient calculation; default=false>
  halogen=<true/false; generate halogen fragmaps; default=false>
  norm=<standard/water/volume; default=water>
  cns=<true/false; generate maps in CNS format; default=false>
  begin=<cycle number to begin for map generation; default=$begin>
  end=<cycle number to end for map generation; default=$end>
  temp=<simulation temperature; default=$temp>
@EOF
}

# default parameters
prot=""
sysname=""
setupdir="${workdir}/1_setup_neutral"
if [ -e ${setupdir}/ligand.txt ]; then lig=`cat ${setupdir}/ligand.txt`; else lig=""; fi
if [ -e ${setupdir}/numsys.txt ]; then numsys=`cat ${setupdir}/numsys.txt`; else numsys=10; fi
workdir=`echo $PWD`
drudedir="2e_run_drude"
qname=""
account=""
email=""
protmap=false
skipoc=false
halogen=false
standard=true
norm="water"
cns="false"
required_args="prot"

parse_args $@

if [[ "$halogen" == "true" ]]; then
  echo "Halogen probes are not supported in SILCS-RNA yet!!!"
  exit 1
  #standard=false
  #fragments="brbx clbx clex fetx flbc flbx tfec"
  #gfefragments="brbx clbx clex fetx flbc flbx tfec"
  #setupdir="1_setup_x"
  #gcmddir="2a_run_gcmd_x"
  #mapdir="2b_gen_maps_x"
  #ocfilename="overlap_coeff_x.dat"
  #probe=halogen
fi
if [[ "$standard" == "true" ]]; then
  halogen=false
  fragments_neutral="meoo imin iminh foro forn forc dmeo benc prpc swm4o gehc"
  fragments_charged="mamn acec aceo swm4o"
  gfefragments_neutral="meoo imin iminh foro forn forc dmeo benc prpc swm4o gehc apolar hbdon hbacc"
  gfefragments_charged="mamn aceo acec swm4o"
  setupdir="1_setup"
  drudeddir="${workdir}/2e_run_drude"
  mapdir="2h_drude_maps"
  ocfilename="oc_drude.dat"
  probe=standard
fi

PROT_PDB=$(basename $prot)
PROT_PDB="${PROT_PDB%.*}"

if [ ! -z "$lig" ];
then
  ligfilename=$lig
  ext="${ligfilename##*.}"

  ligfile=$(basename "${ligfilename}")
  ext="${ligfile##*.}"
  ligfilename="${ligfile%.*}"

  lig="${ligfilename%.*}"
  LIGAND=$ligfilename
  echo $LIGAND $ligfilename $ext
fi

if [ ! -z ${lig} ];
then
  pdb=${sysname}_${LIGAND}
else
  pdb="${sysname%.*}"
fi

#if [[ ! -d "${mapdir}" ]]; then                                                           
#    echo "cannot find mapdir "${mapdir}"; exit"                                           
#    exit 1                                                                                
#  fi 

#echo ${mapdir}
#if [[ ! -d "${mapdir}_neutral" ]] || [[ ! -d "${mapdir}_charged" ]]; then
#  echo "cannot find mapdir ${mapdir}_neutral/charged ; exit"
#  exit 1
#fi

bundle_maps () {
  # combine multiple maps and simply add them together to produce a single map
  # if third argument is passed, the number will be used to divide the sum
  local mapfiles=$1
  local outfile=$2
  local denorm=${3:-1}
#  echo `pwd`
  #echo $mapfiles >$(tty)
  head -6 `echo $mapfiles | cut -d " " -f1` > $outfile
  paste $mapfiles | awk -v denorm=${denorm} '{ if (NR > 6) { sum=0; for (i=1; i<=NF; i++) sum+=$i; printf "%10.3f\n", (sum/denorm) }}' >> $outfile
}

gen_gfe_map () {
  # generate GFE map by normalizing the occupancy maps with normalization constant
  local infile=$1
  local outfile=$2
  local natoms=$3
#  echo "in=" $infile >$(tty)
#  echo "out=" $outfile >$(tty)
#  echo "natoms=" $natoms >$(tty)
  spacing=$(head -6 $infile | grep SPACING | awk '{print $2*$2*$2}')
  #echo "eps is:" $eps >$(tty)
#  echo "bulk="${bulk1} >$(tty)
  awk -v kt=$kt -v bulk1=${bulk1} -v natoms=$natoms -v eps=$eps -v spacing=$spacing 'BEGIN{bulk=bulk1*natoms*spacing} { if (NR<7) { print $0 } else { printf "%10.3f\n", -kt*log(($1+eps)/(bulk+eps)) }}' $infile > $outfile
}

bundle_gfe_maps () {
  # combine multiple maps and generate GFE maps at the same time
  local mapfiles=$1
  local outfile=$2
  local natoms=$3
  #echo "maps=" $mapfiles >$(tty)
  #echo "out=" $outfile >$(tty)
  #echo "natoms=" $natoms >$(tty)
  #echo "out=" $outfile >$(tty)
  head -6 `echo $mapfiles | cut -d " " -f1` > $outfile
  spacing=$(head -6 $outfile | grep SPACING | awk '{print $2*$2*$2}')
  paste $mapfiles | awk -v kt=$kt -v bulk1=${bulk1} -v natoms=$natoms -v eps=$eps -v spacing=$spacing 'BEGIN{bulk=bulk1*natoms*spacing} { if (NR > 6) { sum=0; for (i=1; i<=NF; i++) sum+=$i; printf "%10.3f\n", -kt*log((sum+eps)/(bulk+eps)) }}' >> $outfile
}

bundle_frag_maps () {
  # helper function for bundling up fragment occ maps
  local frags=$1
  local outfile=$2
  local natoms=${3:-1}
  
  mapfiles=""
  
  for frag in $(echo $frags); do
    #echo "natoms=" $natoms >$(tty)
    #echo "frag=" $frag >$(tty)
    mapfiles="${mapfiles} ${mapdir}/${pdb}.${frag}.${startrun}-${endrun}.map"
  done
  #echo "mapfiles=" $mapfiles >$(tty)
  bundle_maps "${mapfiles}" "$outfile" $natoms
}

bundle_frag_gfe_maps () {
  # helper function for bundling up occ maps and create gfe map
  local frags=$1
  local outfile=$2
  local natoms=$3
  mapfiles=""
  for frag in $(echo $frags); do
    mapfiles="${mapfiles} ${mapdir}/${pdb}.${frag}.${startrun}-${endrun}.map"
  done
#  echo "out=" $outfile >$(tty)
  bundle_gfe_maps "$mapfiles" "$outfile" $natoms
}

avg_fragments_md () {

  # count the number of fragments during the Drude simulation
  local frag=$1
  local key=$2
  local nruns=$3
  #echo $frag $key >$(tty)
  nsys=$((nruns*10))
  #echo $endrun $nsys >$(tty)
  numfrag=""
  i="$startrun"
  while [[ "$i" -le "$endrun" ]]; do
    tot=`ls ${drudedir}/${i}/*/ -d| wc -l | awk '{print $1+1}'`
    for j in {10..100..10}; do 
      topfile="${drudedir}/${i}/${j}/build/system_properties.str"
      if [[ -f "$topfile" ]]; then
        if [[ "$frag" == "SWM4" ]]; then
          nf=`grep "nsolv" ${topfile} | awk '{print $3}'`
          numfrag=$((numfrag+nf))
        else
          nf=`grep "n${frag,,}" ${topfile} | awk '{print $3}'`
        #numfrag="$numfrag $nf"
        #echo $nf >$(tty)
          numfrag=$((numfrag+nf))
        fi
      fi
    done
    #echo "function test" $numfrag>$(tty)
    i=$((i+1))
  done
  avg_numfrag=$((numfrag/$nsys))
  echo $avg_numfrag
}

#avg_fragments_md () {
#  # count the number of fragments during the full MD simulation
#  local frag=$1
#  local key=$2
#  numfrag=""
#  i="$startrun"
#  while [[ "$i" -le "$endrun" ]]; do
#      topfile="${drudedir}/${i}/step4_equilibration_drude.crd"
#      if [[ -f "$topfile" ]]; then
#        nf=`grep "$frag" "$topfile" | grep "$key" | wc -l`
#        numfrag="$numfrag $nf"
#      fi
#    i=$((i+1))
#  done
#  ##echo $frag $key $numfrag
#  echo $(avg $numfrag)
#}

avg_volume () {
  # compute average volume
  local nrun=$3
  nsys=$((nrun*10))
  volume=""
  i="$startrun"
  while [[ "$i" -le "$endrun" ]]; do
    tot=`ls ${drudedir}/${i}/*/ -d| wc -l | awk '{print $1+1}'`
    pdbfile="${drudedir}/${i}/${j}/${i}_${j}_or.pdb"
    for j in {10..100..10}; do
      nf=`head "$pdbfile" | grep CRYST1 | awk '{print $2*$3*$4}'`
      #volume="$volume $nf"
      volume=$((volume+nf))
    done
    i=$((i+1))
  done
  echo $((volume/$numsys))
}

map_frag_mol () {
  # return the molecule name for a given fragment name
  local frag=$(tolower $1)
  FRAGARRAY=( "meoo:MEOH"
              "imin:IMID"
              "iminh:IMID"
              "foro:FORM"
              "forn:FORM"
              "forc:FORM"
              "mamn:MAMY"
              "aceo:ACEY"
              "benc:BENX"
              "prpc:PRPX"
              "aalo:AALD"
              "aalc:AALD"
              "gehc:IMID"
              "acec:ACEY"
              "brbc:BRBX"
              "brbx:BRBX"
              "clbc:CLBX"
              "clbx:CLBX"
              "clec:CLEX"
              "clex:CLEX"
              "dmec:DMEE"
              "dmeo:DMEE"
              "fetc:FETX"
              "fetx:FETX"
              "flbc:FLBX"
              "flbx:FLBX"
              "tfec:TFEX"
              "tfex:TFEX"
              "swm4o:SWM4" )

  for fragment in "${FRAGARRAY[@]}" ; do
    FRAGNAME="${fragment%%:*}"
    MOLNAME="${fragment##*:}"
    if [[ "$FRAGNAME" == "$frag" ]]; then
      echo $(toupper "$MOLNAME")
      return
    fi
  done

  echo "couldn't find molname for fragment: $frag" 1>&2
  return
}

map_frag_key () {
  # return the key atom name for a given fragment name
  local frag=$(toupper $1)
  FRAGATOMARRAY=( "MAMY:LPA"
              "ACEY:LPA"
              "BENX:LPA"
              "PRPX:LPA"
              "MEOH:LP1A"
              "DMEE:LP2A"
              "FORM:LPA"
              "IMID:LP1"
              "SWM4:OM" )

  for fragment in "${FRAGATOMARRAY[@]}" ; do
    FRAGNAME="${fragment%%:*}"
    KEYNAME="${fragment##*:}"
    if [[ "$FRAGNAME" == "$frag" ]]; then
      echo $(toupper "$KEYNAME")
      return
    fi
  done

  echo "couldn't find keyname for fragment: $frag" 1>&2
  return
}

density_prefac () {
  density=$1
  nsnapshots=$2
  avogadro=0.0006023
  #echo $density $nsnapshots >$(tty)
  bulk1=$(echo 1 | awk -v frag_density=${density} -v nsnapshots=${nsnapshots} -v avogadro=${avogadro} '{prefac=frag_density*nsnapshots*avogadro ; printf "%f\n",prefac}')
  # this is a hidden return. Do not comment out the next line
  echo $bulk1
}

collate_maps () {
  local drudedir=$1
  local mapdir=$2
  local startrun=$3
  local endrun=$4
  local protmap=$5
  local probe=$6
  #echo "Build maps combining data from : " $startrun"-"$endrun

  #nframes=101
  totcycles=0
  i=$startrun
  while [[ "$i" -le "$endrun" ]];
  do
    #tot=`find ./${drudedir}/${i}/ -name ${pdb}_silcs.${i}.prod.*.rec.pdb | wc -l`
    tot=2000
    totcycles=$((totcycles+tot))
    #echo $i $tot $totcycles >$(tty)
    i=$((i+1))
    
  done

  #nsnapshots=`echo $totcycles $nframes | awk '{printf "%d\n", $1*$2}'`
  nsnapshots=$totcycles

  if [[ "$nsnapshots" -eq 0 ]]; then
    echo "couldn't find GCMC trajectory files; exit"
    exit 1
  fi

  ## Collate Maps

  allfragments="$fragments"
  if [[ "$protmap" == true ]]; then
    allfragments="$fragments ala ser thr cys val leu ile met pro phe tyrc tyro trpc trpn asp glu asnn asno glnn glno hsdn hsdn2 hsen hsen2 hsp arg lys"
  fi
  for frag in $(echo $allfragments) ; do
    #echo ${frag} >$(tty)
    i=${startrun}
    mapfiles=""
    while [[ "$i" -le "${endrun}" ]];
    do
      mapfiles="${mapfiles} ${mapdir}/maps_frm_runs/${pdb}.${frag}.${i}.map"
      i=$((i+1))
    done
    bundle_maps "${mapfiles}" "${mapdir}/${pdb}.${frag}.${startrun}-${endrun}.map"
  done

  #GFE Maps
  frag_density=0.001138
  nfrags=8
  fragbulk=`echo 1 | awk -v frag_density=${frag_density} -v nsnapshots=${nsnapshots} -v nfrags=${nfrags} '{prefac=frag_density*nsnapshots/nfrags ; printf "%f\n",prefac}'`
  #echo "bulk1" $bulk1

  #Wat_density value is again obtained from the frag-only simulation. Its different from TierI value (which was ~0.028, causing a minor ~ 0.05 kcal/mol diff )
  wat_density=0.030539
  watbulk=`echo 1 | awk -v frag_density=${wat_density} -v nsnapshots=${nsnapshots} '{prefac=frag_density*nsnapshots ; printf "%f\n",prefac}'`
  #echo "watbulk" $watbulk

  kt=0.592
  eps=0.1

  if [[ "$norm" == "water" ]]; then
    numwat=$(avg_fragments_md "SWM4" "OM" "$(((endrun-startrun)+1))" | awk '{print $1/55}')
    #echo "numberofwater" $numwat
    frag_density=""
  fi
  if [[ "$norm" == "volume" ]]; then
    volume=$(avg_volume "${endrun}")
    frag_density=""
  fi

  for frag in $(echo $fragments); do
    natoms=0
    bulk1=${fragbulk}
    natoms=`head -2 ${SILCSBIODIR}/data/drude_selection/${frag}.txt | tail -1`
    if [[ "$frag" == swm4o ]]; then
      bulk1=${watbulk}
    fi
    if [[ "$natoms" -eq 0 ]]; then continue; fi
#    echo "before gfe" >$(tty)
    gfemap="${mapdir}/${pdb}.${frag}.${startrun}-${endrun}.gfe.map"
    if [[ "$norm" == "water" ]]; then
      molname=$(map_frag_mol "$frag")
      #echo $molname
      keyname=$(map_frag_key "$molname")
      #echo $keyname
      numfrag=$(avg_fragments_md "$molname" "$keyname" "$(((endrun-startrun)+1))")
      #echo $molname $numfrag
      #echo "eachfrag"  $molname $numfrag
      density=$(echo $numfrag $numwat | awk '{print $1/$2}')
      #echo $frag $density $nsnapshots >$(tty)
      bulk1=$(density_prefac $density $nsnapshots)
      frag_density="${frag_density} ${frag}:${density}"
    fi
    if [[ "$norm" == "volume" ]]; then
      molname=$(map_frag_mol "$frag")
      numfrag=$(avg_fragments_md "$molname")
      density=$(density "$numfrag" "$volume")
      bulk1=$(density_prefac $density $nsnapshots)
      frag_density="${frag_density} ${frag}:${density}"
    fi
    #echo $frag $frag_density $gfemap $natoms >$(tty)
 #   echo $bulk1 >$(tty) 
    gen_gfe_map "${mapdir}/${pdb}.${frag}.${startrun}-${endrun}.map" $gfemap $natoms 
  done
  #echo "after gfe"
  #
  # building generic maps
  #
  bulk1=${fragbulk}

  avg_bulk () {
    ave_density=""
    n=0
    for frag in $(echo $frags); do
      for fragment in $(echo $frag_density) ; do
        FRAGNAME="${fragment%%:*}"
        DENSITY="${fragment##*:}"
        #echo $frags >$(tty)
        #echo $FRAGNAME $DENSITY >$(tty)
        if [[ "$frag" == "$FRAGNAME" ]]; then
          #echo $frag >$(tty)
          tot_density=`echo $tot_density $DENSITY | awk '{print $1 + $2}'`
          #echo $fragment ${n} ${tot_density} >$(tty)
          n=$((n+1))
          break
        fi
      done
    done
    #echo ${n} >$(tty)
    ave_density=`echo ${tot_density} ${n} | awk '{print $1 / $2}'`
    #echo ${ave_density} >$(tty)
    bulk1=$(density_prefac $ave_density $nsnapshots)
    echo $bulk1 
  }
  #echo $probe >$(tty)
  if [[ "$probe" == "neutral" ]]; then
#	  echo "neutral" >$(tty)
    if [[ "$standard" == "true" ]]; then
      frags="foro dmeo imin"
      gfemap="${mapdir}/${pdb}.hbacc.${startrun}-${endrun}.gfe.map"
      if [[ "$norm" == "water" ]]; then
        bulk1=$(avg_bulk)
      fi
      if [[ "$norm" == "volume" ]]; then
        bulk1=$(avg_bulk)
      fi
      bundle_frag_gfe_maps "$frags" "$gfemap" 3

      frags="forn iminh"
      gfemap="${mapdir}/${pdb}.hbdon.${startrun}-${endrun}.gfe.map"
      if [[ "$norm" == "water" ]]; then
        bulk1=$(avg_bulk)
      fi
      if [[ "$norm" == "volume" ]]; then
        bulk1=$(avg_bulk)
      fi
      bundle_frag_gfe_maps "$frags" "$gfemap" 2

      frags="benc prpc gehc"
      gfemap="${mapdir}/${pdb}.apolar.${startrun}-${endrun}.gfe.map"
      if [[ "$norm" == "water" ]]; then
        bulk1=$(avg_bulk)
      fi
      if [[ "$norm" == "volume" ]]; then
        bulk1=$(avg_bulk)
      fi
    #echo $gfemap
      bundle_frag_gfe_maps "$frags" "$gfemap" 12 
    fi
  fi

  #
  # excl & zero maps
  #
  bundle_frag_maps "$fragments" "allprobes.map"

  cat allprobes.map | awk '{if(NR<=6) {print $0} else {s=0 ;if($1<0.9){s=1000}else{s=0} ; printf "%10.3f\n",s}}' > ${mapdir}/${pdb}.excl.${startrun}-${endrun}.map

  #rm -f allprobes.map

  cat ${mapdir}/${pdb}.${frag}.${startrun}-${endrun}.map | awk '{if(NR<=6)print $0 ; else{printf "%10.3f\n",0}}' > ${mapdir}/${pdb}.zero.map

  #
  # create softlinks
  #
  for frag in $(echo $gfefragments); do
  #echo $gfefragments >$(tty)
    if [[ ! -f ${pdb}.${frag}.map ]]; then
      if [[ -f "${pdb}.${frag}.${startrun}-${endrun}.map" ]]; then
        rm -f ${mapdir}/${pdb}.${frag}.map
        ln -s ${pdb}.${frag}.${startrun}-${endrun}.map ${mapdir}/${pdb}.${frag}.map;
      fi
    fi
    if [[ ! -f ${pdb}.${frag}.gfe.map ]]; then
      rm -f ${mapdir}/${pdb}.${frag}.gfe.map
      ln -s ${pdb}.${frag}.${startrun}-${endrun}.gfe.map ${mapdir}/${pdb}.${frag}.gfe.map;
    fi
  done

  if [[ ! -f ${pdb}.excl.map ]]; then
    rm -f ${mapdir}/${pdb}.excl.map
    ln -s ${pdb}.excl.${startrun}-${endrun}.map ${mapdir}/${pdb}.excl.map;
  fi

  #
  # collate protein side-chain map
  #
  if [[ "$protmap" == true ]]; then

    frags="ala val leu ile met pro phe tyrc trpc"
    bundle_frag_maps "$frags" ${mapdir}/${pdb}.pro.apolar.${startrun}-${endrun}.map

    frags="ser thr cys tyro"
    bundle_frag_maps "$frags" ${mapdir}/${pdb}.pro.meoo.${startrun}-${endrun}.map

    frags="trpn asnn glnn hsdn2 hsen2"
    bundle_frag_maps "$frags" ${mapdir}/${pdb}.pro.hbdon.${startrun}-${endrun}.map

    frags="asno glno hsdn hsen"
    bundle_frag_maps "$frags" ${mapdir}/${pdb}.pro.hbacc.${startrun}-${endrun}.map

    frags="asp glu"
    bundle_frag_maps "$frags" ${mapdir}/${pdb}.pro.aceo.${startrun}-${endrun}.map

    frags="hsp arg lys"
    bundle_frag_maps "$frags" ${mapdir}/${pdb}.pro.mamn.${startrun}-${endrun}.map

    frags="apolar hbdon hbacc meoo aceo mamn"
    for frag in $(echo $frags); do
      if [[ ! -f ${pdb}.${frag}.pro.map ]]; then
        rm -f ${mapdir}/${pdb}.${frag}.pro.map
        ln -s ${pdb}.pro.${frag}.${startrun}-${endrun}.map ${mapdir}/${pdb}.pro.${frag}.map;
      fi
    done
  fi
}

compute_overlap_coefficient() {
  IFS=' ' read -r -a fragorder <<< "${fragments}"
  local mapdir=$1
  set1=$2
  set2=$3
  orderindex=0
  totfrag=${#fragorder[@]}
  while [[ $orderindex -lt $totfrag ]];do
  #echo `pwd` >$(tty)
  #echo ${mapdir}/${pdb}.${fragorder[$orderindex]}.${set1}.map >$(tty) 
    oc=`paste ${mapdir}/${pdb}.${fragorder[$orderindex]}.${set1}.map \
              ${mapdir}/${pdb}.${fragorder[$orderindex]}.${set2}.map \
        | awk 'BEGIN {oc=0}
               NR>6 {sum1=sum1+$1; sum2=sum2+$2; set1[NR]=$1; set2[NR]=$2}
               END {for (i=7; i<=NR; i++) {if (set1[i]<set2[i]) { oc += set1[i]/sum1 } else { oc += set2[i]/sum2 }}; printf "%3.6f", oc}'`
    #echo $sum1 $sum2 >$(tty)
    if [[ $orderindex -eq 0 ]]; then
      echo -e ${fragorder[$orderindex]} "\t" $oc > $ocfilename
    else
      echo -e ${fragorder[$orderindex]} "\t" $oc >> $ocfilename
    fi
    #echo ${fragments} >$(tty)
    echo $overlapcoeff ${fragorder[$orderindex]}| awk '{if($1 < 0.5){printf "For %s fragment, OC = %3.6f, please check your simulations\n", $2,$1}}'

    orderindex=$((orderindex + 1))
  done

  #rm -f a.nom.map b.nom.map
}

echo "########################################################"
echo "prot         = " ${prot}
echo "system       = " ${pdb}
echo "probe        = " ${probe}
echo "temp         = " ${temp}
echo "norm         = " ${norm}
echo "protmap      = " ${protmap}
echo "numsys       = " ${numsys}
echo "cns          = " ${cns}
echo "Will collate maps from individual runs to generate a visualization directory : silcs_fragmaps_${pdb} "
echo "For subsequent steps (SILCS-MC, etc), use maps from silcs_fragmaps_${pdb} directory"
echo "########################################################"
echo -e "\n"

if [[ "$numsys" -gt 1 && "$skipoc" == "false" ]]; then
  startrun=1
  midrun1=$((numsys/2 + numsys%2))
  echo "round 1 : collate maps from runs : ${startrun}-${midrun1}"
  for probe in neutral charged; do
      echo " --> $probe"
      if [ $probe == "neutral" ]; then fragments="${fragments_neutral}"; elif [ $probe == "charged" ]; then fragments="${fragments_charged}"; fi
      collate_maps ${drudedir}_${probe} ${mapdir}_${probe} $startrun $midrun1 false $begin $end $probe
  done

  midrun2=$((numsys/2 + 1))
  endrun=$((numsys))
  echo "round 2 : collate maps from runs : ${midrun2}-${endrun}"
  for probe in neutral charged; do
      echo " --> $probe"
      if [ $probe == "neutral" ]; then fragments="${fragments_neutral}"; elif [ $probe == "charged" ]; then fragments="${fragments_charged}"; fi
      collate_maps ${drudedir}_${probe} ${mapdir}_${probe} $midrun2 $endrun false $begin $end $probe
  done

  echo "Running QC to ensure simulations converged; will calculate overlap coefficients (OC) for each of the probe maps."

  set1="1-${midrun1}"
  set2="${midrun2}-${numsys}"
  for probe in neutral charged; do
      echo " --> $probe"
      if [ $probe == "neutral" ]; then fragments="${fragments_neutral}"; elif [ $probe == "charged" ]; then fragments="${fragments_charged}"; fi
      compute_overlap_coefficient ${mapdir}_${probe} $set1 $set2
      mv oc_drude.dat oc_drude_${probe}.dat
  done

  echo "Finished QC check; creating directory silcs_fragmaps_${pdb} for visualization"
  echo ""

else
  if [[ "$skipoc" == "true" ]]; then
    echo "skipping overlap coefficient calculation"
  else
    echo "too little numsys; skipping overlap coefficient calculation"
  fi
fi

#
# collate maps
#
startrun=1
endrun=$numsys
echo "final: collate maps from runs : ${startrun}-${endrun}"
echo ""
for probe in neutral charged; do
    echo " --> $probe"
    if [ $probe == "neutral" ]; then fragments="${fragments_neutral}"; elif [ $probe == "charged" ]; then fragments="${fragments_charged}"; fi
    collate_maps ${drudedir}_${probe} ${mapdir}_${probe} $startrun $endrun false $begin $end $probe
done

#
# collate maps
#
fragmapdir="${pdb}_drude_fragmaps"

# Create a new visualization directory

if [[ "$standard" == "true" ]]; then
  for probe in neutral charged; do
  if [ ! -d ${fragmapdir}_${probe}/maps ]; then mkdir -p ${fragmapdir}_${probe}/maps; fi
  if [ $probe == "neutral" ]; then gfefragments="${gfefragments_neutral}"; elif [ $probe == "charged" ]; then gfefragments="${gfefragments_charged}"; fi
  for frag in $(echo "${gfefragments}") ; do
    if [[ -f ${mapdir}_${probe}/${pdb}.${frag}.${startrun}-${endrun}.gfe.map ]]; then
      rm -f ${fragmapdir}_${probe}/maps/${pdb}.${frag}.gfe.map
      cp ${mapdir}_${probe}/${pdb}.${frag}.${startrun}-${endrun}.gfe.map ${fragmapdir}_${probe}/maps/${pdb}.${frag}.gfe.map
    fi
  done

  if [[ -f ${mapdir}_${probe}/${pdb}.excl.map ]]; then
    rm -f ${fragmapdir}_${probe}/maps/${pdb}.excl.map
    cp ${mapdir}_${probe}/${pdb}.excl.${startrun}-${endrun}.map ${fragmapdir}_${probe}/maps/${pdb}.excl.map
  fi

  if [[ -f ${mapdir}_${probe}/${pdb}.zero.map ]]; then
    rm -f ${fragmapdir}_${probe}/maps/${pdb}.zero.map
    cp ${mapdir}_${probe}/${pdb}.excl.${startrun}-${endrun}.map ${fragmapdir}_${probe}/maps/${pdb}.zero.map
  fi
  done
  if [[ "$protmap" == true ]]; then
    fragments="apolar hbdon hbacc meoo"
    probe="neutral"
    for frag in $(echo $fragments); do
      if [[ -f ${mapdir}_${probe}/${pdb}.pro.${frag}.${startrun}-${endrun}.map ]]; then
        rm -f ${fragmapdir}_${probe}/maps/${pdb}.pro.${frag}.map
        cp ${mapdir}_${probe}/${pdb}.pro.${frag}.${startrun}-${endrun}.map ${fragmapdir}_${probe}/maps/${pdb}.pro.${frag}.map
      fi
    done
    fragments="mamn aceo acec"
    probe="charged"
    for frag in $(echo $fragments); do
      if [[ -f ${mapdir}_${probe}/${pdb}.pro.${frag}.${startrun}-${endrun}.map ]]; then
        rm -f ${fragmapdir}_${probe}/maps/${pdb}.pro.${frag}.map
        cp ${mapdir}_${probe}/${pdb}.pro.${frag}.${startrun}-${endrun}.map ${fragmapdir}_${probe}/maps/${pdb}.pro.${frag}.map
      fi
    done
  fi
fi
if [[ "$halogen" == "true" ]]; then
  echo "Halogen probes are not supported in SILCS-RNA yet!!!"
  exit 1
  ##if [ ! -d ${fragmapdir}/maps ]; then mkdir -p ${fragmapdir}/maps; fi
  ##
  ##for frag in $(echo $gfefragments) ; do
  ##  if [[ -f ${mapdir}/${pdb}.${frag}.${startrun}-${endrun}.gfe.map ]]; then
  ##    rm -f ${fragmapdir}/maps/${pdb}.${frag}.gfe.map
  ##    cp ${mapdir}/${pdb}.${frag}.${startrun}-${endrun}.gfe.map ${fragmapdir}/maps/${pdb}.${frag}.gfe.map
  ##  fi
  ##done
fi

if [[ -f ${pdb}.pdb ]]; then
  for probe in neutral charged; do
    rm -f ${fragmapdir}_${probe}/${pdb}.pdb
    cp ${pdb}.pdb ${fragmapdir}_${probe}/${pdb}.pdb
  done
fi

# Update plugins
for probe in neutral charged; do
  sed -e "s/<prot>/${pdb}/g" ${SILCSBIO_TEMPLATE_DIR}/silcs-rna/vmd_fragmap.tmpl > ${fragmapdir}_${probe}/view_maps.vmd
  sed -e "s/<prot>/${pdb}/g" ${SILCSBIO_TEMPLATE_DIR}/silcs-rna/pymol_fragmap.tmpl > ${fragmapdir}_${probe}/view_maps.pml
  rm -rf ${fragmapdir}_${probe}/plugins
  cp -r ${SILCSBIODIR}/utils/plugins ${fragmapdir}_${probe}/
done

# Merge neutral and charged maps
cp -rf ${fragmapdir}_neutral ${fragmapdir}
cp ${fragmapdir}_charged/maps/${pdb}.mamn.gfe.map ${fragmapdir}/maps/${pdb}.mamn.gfe.map
cp ${fragmapdir}_charged/maps/${pdb}.acec.gfe.map ${fragmapdir}/maps/${pdb}.acec.gfe.map
cp ${fragmapdir}_charged/maps/${pdb}.aceo.gfe.map ${fragmapdir}/maps/${pdb}.aceo.gfe.map

# CNS maps
if [[ "$cns" == "true" ]]; then
  echo "converting maps to CNS format"
  ${SILCSBIODIR}/utils/map_to_cns.sh prot=${prot} mapsdir=${fragmapdir}/maps > /dev/null
  sed -e "s/<prot>/${pdb}/g" ${SILCSBIO_TEMPLATE_DIR}/silcs-rna/moe_fragmap.tmpl > ${fragmapdir}/view_maps.svl
fi

# Copy the version number of SILCSBIO to the fragmaps map folder
cp ${SILCSBIODIR}/VERSION ${fragmapdir}/maps/.

cat <<@EOF
########################################################
Copy the folder '${pdb}_drude_fragmaps' to your local workstation to visualize the FragMap.
To visualize the FragMaps :

cd ${pdb}_drude_fragmaps
vmd -e view_maps.vmd

or

pymol view_maps.pml

For subsequent steps (SILCS-MC, SILCS-PHARM, etc), use maps from drude_fragmaps_${pdb} directory
########################################################
@EOF

