Drude SILCS Scripts
==================

This repository contains scripts for running **Drude polarizable simulations**
in **OpenMM** as part of the Drude SILCS-Nucleic workflow. 

The workflow is designed for nucleic-acid systems (DNA and RNA), with support for ions and multiple chains. **Drude SILCS-Nucleic** was developed as an extension of the additive SILCS process using **CHARMM36 (C36)**. As such, Drude polarizable simulation must be run only after PDB files have been generated from the additive SILCS workflow. 


The force field files provided are developmental but correspond to the files used in the Drude SILCS-Nucleic paper. Official force field files should be obtained from http://mackerell.umaryland.edu for all other purposes.

## Overview
**SILCS (Site Identification by Ligand Competitive Saturation)** is a physics-based method for mapping functional group interactions around biomolecular targets. This repository extends the SILCS methodology to nucleic acids using the Drude polarizable FF, for improved treatment of electrostatic and polarization effects. 

The scripts here automate key steps required to:
- Prepare nucleic acid systems for Drude simulations
- Run OpenMM-based Drude SILCS simulations
- Generate FragMaps based on simulation trajectories

---

## Repository Structure

## Dependencies 

The Drude SILCS-Nucleic workflow requirese the following software and libraries

### Required Software
- SILCSBio v2023.5
- OpenMM v
- Python
- LOOS
- Anaconda3
- CHARMM

### Python Packages

The following Python packages are required (may be installed via `pip` or `conda`:

-`numpy`
-`scipy`
-`mdanalysis`
-`openmm`

### Hardware Notes

- GPU acceleration is strongly recommended for Drude simulations
  
## Usage

### Step 2e: Build Drude Systems
```
$ bash $SILCSBIODIR/drude_silcs/silcs-rna/2e_drude_build.sh prot=system.pdb sysname=pdb_code toppardir=/expanse/lustre/projects/vpt121/hmichel/toppar_drude_silcs_rna qname=gpu-shared account=vpt121 email=user@example.com
````
### Step 2f: Drude Equilibration
````
$ bash $SILCSBIODIR/drude_silcs/silcs-rna/2f_drude_equil.sh prot=system.pdb sysname=pdb_code toppardir=/path/to/ff_files qname=partition_name account=your_account email=user@example.com
````
### Step 2g: Drude Production
````
$ bash $SILCSBIODIR/drude_silcs/silcs-rna/2g_drude_prod.sh prot=system.pdb sysname=pdb_code fcsilcs=50.208 customres="segid DNAA and (resid 3 or resid 4 or resid 14) and name P" qname=partition_name account=your_account email=user@example.com
````

### Step 2h: Generate Drude FragMaps
````
$ bash ${SILCSBIODIR}/drude_silcs/silcs-rna/2h_gen_maps_drude.sh prot=system.pdb sysname=pdb_code qname=partition_name account=your_account email=user@example.com
````

### Step 2i: Combine Drude FragMaps
````
$ bash ${SILCSBIODIR}/drude_silcs/silcs-rna/2i_drude_combine_maps.sh prot=system.pdb sysname=pdb_code qname=partition_name account=your_account email=user@example.com
````

Principal Authors and Contributors
----------------------------------

Haley M. Michel<br>
Virginia Tech Department of Biochemistry<br>
hmichel@vt.edu

Justin A. Lemkul<br>
Virginia Tech Department of Biochemistry<br>
jalemkul@vt.edu

Alexander D. MacKerell, Jr.<br>
University of Maryland, Baltimore Department of Pharmaceutical Sciences<br>
alex@outerbanks.umaryland.edu
