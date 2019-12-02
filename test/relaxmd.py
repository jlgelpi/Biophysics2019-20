#!/usr/bin/env python
# coding: utf-8

# Based on the official Gromacs tutorial: http://www.mdtutorials.com/gmx/lysozyme/01_pdb2gmx.html

import sys

#Import module
from biobb_io.api.pdb import Pdb

downloaded_pdb_path = sys.argv[1]

base_path = sys.argv[2]

# ## Modeling the missing heavy atoms in the structure side chains
# This BB will reconstruct missing atoms from residue side chains and detect all kinds of clashes
#Check & Fix PDB
#Import module
from biobb_model.model.fix_side_chain import FixSideChain
# Create prop dict and inputs/outputs
fixed_pdb_path = base_path + '/fixed.pdb'
#Create and launch bb
FixSideChain(input_pdb_path=downloaded_pdb_path, 
             output_pdb_path=fixed_pdb_path).launch()

# ## Generate the topology
# Generate the topology using [Pdb2gmx](http://manual.gromacs.org/documentation/2019/onlinehelp/gmx-pdb2gmx.html) module.
# The default forcefield is [amber99sb-ildn](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2970904/) and the default water_type is [spce](http://www.sklogwiki.org/SklogWiki/index.php/SPC/E_model_of_water).
# This BB will generate 2 main files:
# -  A compressed ZIP file containing:
# > -  The [Gromacs topology file](http://manual.gromacs.org/documentation/2019/reference-manual/file-formats.html#top) (.top)
# > -  The [Gromacs position restraint file](http://manual.gromacs.org/documentation/2019/reference-manual/file-formats.html#itp) (.itp)
# -  A post-processed [Gromacs structure file](http://manual.gromacs.org/documentation/2019/reference-manual/file-formats.html#gro) (.gro)

#Create system topology
# Import module
from biobb_md.gromacs.pdb2gmx import Pdb2gmx
# Create inputs/outputs
output_pdb2gmx_gro_path = base_path + '/pdb2gmx.gro'
output_pdb2gmx_top_zip_path = base_path + '/pdb2gmx_top.zip'
#Create and launch bb
Pdb2gmx(input_pdb_path=fixed_pdb_path, 
        output_gro_path=output_pdb2gmx_gro_path, 
        output_top_zip_path=output_pdb2gmx_top_zip_path).launch()

# Note that hydrogen atoms had been added to the structure.
# ***

# ## Create the solvent box
# Create the solvent box using the [Editconf](http://manual.gromacs.org/documentation/2019/onlinehelp/gmx-editconf.html) module. The box will be centered, the distance between the solute and the box will be 1.0nm and the box shape will be cubic by default. The main output of this BB will be an updated post-procesed Gromacs structure file (.gro). 

# Editconf: Create solvent box
# Import module
from biobb_md.gromacs.editconf import Editconf
# Create prop dict and inputs/outputs
output_editconf_gro_path = base_path + '/editconf.gro'
#Create and launch bb
Editconf(input_gro_path=output_pdb2gmx_gro_path, 
         output_gro_path=output_editconf_gro_path).launch()


# ## Fill the solvent box with water molecules
# Fill the solvent box using the Gromacs Solvate module. The "spc216.gro" will be default solvent model. The main output of this BB will be an updated post-procesed Gromacs structure file (.gro) and a zip file containing the updated topology file (.top) and the restriction files (.itp).

# Solvate: Fill the box with water molecules
from biobb_md.gromacs.solvate import Solvate
# Create prop dict and inputs/outputs
output_solvate_gro_path = base_path + '/solvate.gro'
output_solvate_top_zip_path = base_path + '/solvate_top.zip'
#Create and launch bb
Solvate(input_solute_gro_path=output_editconf_gro_path, 
        output_gro_path=output_solvate_gro_path, 
        input_top_zip_path=output_pdb2gmx_top_zip_path, 
        output_top_zip_path=output_solvate_top_zip_path).launch()

# ## Preprocess ion generation
# Create the [portable binary run file (.tpr)](http://manual.gromacs.org/documentation/2019/reference-manual/file-formats.html#tpr) for ion generation, using the Grompp module. The main default parameters for this execution are:
# -  integrator  = steep         ; Algorithm (steep = steepest descent minimization)
# -  emtol       = 1000.0        ; Stop minimization when the maximum force < 1000.0 kJ/mol/nm
# -  emstep      = 0.01          ; Minimization step size
# -  nsteps      = 50000         ; Maximum number of (minimization) steps to perform
# 
# The main output of this BB will be the portable binary run file (.tpr).

# Grompp: Creating portable binary run file for ion generation
from biobb_md.gromacs.grompp import Grompp
# Create prop dict and inputs/outputs
prop = {'mdp':{'type': 'minimization', 'nsteps':'5000'}}
output_gppion_tpr_path =  base_path + '/gppion.tpr'
#Create and launch bb
Grompp(input_gro_path=output_solvate_gro_path, 
       input_top_zip_path=output_solvate_top_zip_path, 
       output_tpr_path=output_gppion_tpr_path,  
       properties=prop).launch()

# ## Ion generation
# Replace solvent molecules to neutralice the system and then reach a 0.05 mol/litre concentration by default. Using the Genion module.
# The main output of this BB will be an updated post-procesed Gromacs structure file (.gro) and a zip file containing the updated topology file (.top) and the restriction files (.itp).

# Genion: Adding ions to reach a 0.05 nm concentration
from biobb_md.gromacs.genion import Genion
# Create prop dict and inputs/outputs
prop={'neutral':True, 'concentration':0.05}
output_genion_gro_path = base_path + '/genion.gro'
output_genion_top_zip_path = base_path + '/genion_top.zip'
#Create and launch bb
Genion(input_tpr_path=output_gppion_tpr_path, 
       output_gro_path=output_genion_gro_path, 
       input_top_zip_path=output_solvate_top_zip_path, 
       output_top_zip_path=output_genion_top_zip_path, 
       properties=prop).launch()

# ## Preprocess energy minimization
# Create the [portable binary run file (.tpr)](http://manual.gromacs.org/documentation/2019/reference-manual/file-formats.html#tpr) for energy minimization, using the Grompp module. The main default parameters for this execution are:
# -  integrator  = steep         ; Algorithm (steep = steepest descent minimization)
# -  emtol       = 500.0         ; Stop minimization when the maximum force < 1000.0 kJ/mol/nm
# -  emstep      = 0.01          ; Minimization step size
# -  nsteps      = 50000         ; Maximum number of (minimization) steps to perform
# 
# The main output of this BB will be the portable binary run file (.tpr).

# Grompp: Creating portable binary run file for mdrun
from biobb_md.gromacs.grompp import Grompp
# Create prop dict and inputs/outputs
prop = {'mdp':{'type': 'minimization', 'nsteps':'5000', 'emtol':'500'}}
output_gppmin_tpr_path = base_path + '/gppmin.tpr'
#Create and launch bb
Grompp(input_gro_path=output_genion_gro_path, 
       input_top_zip_path=output_genion_top_zip_path, 
       output_tpr_path=output_gppmin_tpr_path,  
       properties=prop).launch()

# ## Execute system equilibration
# Execute the energy minimization using the MDrun module and the input the portable binary run file (.tpr) as the main input.
# The main output of this BB will be updated structure file (.gro) and the [Gromacs trajectory file](http://manual.gromacs.org/documentation/2019/reference-manual/file-formats.html#trr) (.trr).

# Mdrun: Running minimization
from biobb_md.gromacs.mdrun import Mdrun
# Create prop dict and inputs/outputs
output_min_trr_path =  base_path + '/min.trr'
output_min_gro_path =  base_path + '/min.gro'
output_min_edr_path =  base_path + '/min.edr'
output_min_log_path =  base_path + '/min.log'
#Create and launch bb
Mdrun(input_tpr_path=output_gppmin_tpr_path, 
      output_trr_path=output_min_trr_path, 
      output_gro_path=output_min_gro_path, 
      output_edr_path=output_min_edr_path, 
      output_log_path=output_min_log_path).launch()

# ## Preprocess system temperature equilibration
# Equilibrate the solvent and ions restraining the protein heavy atoms
# Create the [portable binary run file (.tpr)](http://manual.gromacs.org/documentation/2019/reference-manual/file-formats.html#tpr) for system equilibration, using the Grompp module. The main default parameters for this execution are:
# -  Define                   = -DPOSRES
# -  integrator               = md
# -  dt                       = 0.002
# -  nsteps                   = 5000
# -  pcoupl                   = no
# -  gen_vel                  = yes
# -  gen_temp                 = 300
# -  gen_seed                 = -1
# 
# The main output of this BB will be the portable binary run file (.tpr).

## Grompp: Creating portable binary run file for Equilibration
#from biobb_md.gromacs.grompp import Grompp
## Create prop dict and inputs/outputs
#prop = {'mdp':{'type': 'nvt', 'nsteps':'5000'}}
#output_gppnvt_tpr_path = '1aki_gppnvt.tpr'
##Create and launch bb
#Grompp(input_gro_path=output_min_gro_path, 
#       input_top_zip_path=output_genion_top_zip_path, 
#       output_tpr_path=output_gppnvt_tpr_path,  
#       properties=prop).launch()
#
#
## ## Execute system temperature equilibration
## Execute the system equilibration using the MDrun module and the input the portable binary run file (.tpr) as the main input.
## The main output of this BB will be updated structure file (.gro) and the [Gromacs trajectory file](http://manual.gromacs.org/documentation/2019/reference-manual/file-formats.html#trr) (.trr).
#
## In[ ]:
#
#
## Mdrun: Running Equilibration NVT
#from biobb_md.gromacs.mdrun import Mdrun
## Create prop dict and inputs/outputs
#output_nvt_trr_path = '1aki_nvt.trr'
#output_nvt_gro_path = '1aki_nvt.gro'
#output_nvt_edr_path = '1aki_nvt.edr'
#output_nvt_log_path = '1aki_nvt.log'
#output_nvt_cpt_path = '1aki_nvt.cpt'
##Create and launch bb
#Mdrun(input_tpr_path=output_gppnvt_tpr_path, 
#      output_trr_path=output_nvt_trr_path, 
#      output_gro_path=output_nvt_gro_path, 
#      output_edr_path=output_nvt_edr_path, 
#      output_log_path=output_nvt_log_path, 
#      output_cpt_path=output_nvt_cpt_path).launch()
#
#
## ## Preprocess system pressure equilibration
## Equilibrate the solvent and ions restraining the protein heavy atoms.
## Create the [portable binary run file (.tpr)](http://manual.gromacs.org/documentation/2019/reference-manual/file-formats.html#tpr) for system equilibration, using the Grompp module. The main default parameters for this execution are:
## -  Define                   = -DPOSRES
## -  integrator               = md
## -  dt                       = 0.002
## -  nsteps                   = 5000
## -  pcoupl = Parrinello-Rahman
## -  pcoupltype = isotropic
## -  tau_p = 1.0
## -  ref_p = 1.0
## -  compressibility = 4.5e-5
## -  refcoord_scaling = com
## -  gen_vel = no
## 
## The main output of this BB will be the portable binary run file (.tpr).
#
## In[ ]:
#
#
## Grompp: Creating portable binary run file for mdrun
#from biobb_md.gromacs.grompp import Grompp
## Create prop dict and inputs/outputs
#prop = {'mdp':{'type': 'npt', 'nsteps':'5000'}}
#output_gppnpt_tpr_path = '1aki_gppnpt.tpr'
##Create and launch bb
#Grompp(input_gro_path=output_nvt_gro_path, 
#       input_top_zip_path=output_genion_top_zip_path, 
#       output_tpr_path=output_gppnpt_tpr_path, 
#       input_cpt_path=output_nvt_cpt_path,  
#       properties=prop).launch()
#
#
## ## Execute system pressure equilibration
## Execute the system equilibration using the MDrun module and the input the portable binary run file (.tpr) as the main input.
## The main output of this BB will be updated structure file (.gro) and the [Gromacs trajectory file](http://manual.gromacs.org/documentation/2019/reference-manual/file-formats.html#trr) (.trr).
#
## In[ ]:
#
#
## Mdrun: Running minimization NPT
#from biobb_md.gromacs.mdrun import Mdrun
## Create prop dict and inputs/outputs
#output_npt_trr_path = '1aki_npt.trr'
#output_npt_gro_path = '1aki_npt.gro'
#output_npt_edr_path = '1aki_npt.edr'
#output_npt_log_path = '1aki_npt.log'
#output_npt_cpt_path = '1aki_npt.cpt'
##Create and launch bb
#Mdrun(input_tpr_path=output_gppnpt_tpr_path, 
#      output_trr_path=output_npt_trr_path, 
#      output_gro_path=output_npt_gro_path, 
#      output_edr_path=output_npt_edr_path, 
#      output_log_path=output_npt_log_path, 
#      output_cpt_path=output_npt_cpt_path).launch()
#
#
## ## Preprocess free dynamics
## Free molecular dynamics simulation.
## Create the [portable binary run file (.tpr)](http://manual.gromacs.org/documentation/2019/reference-manual/file-formats.html#tpr) for free molecular dynamics simulation, using the Grompp module. The main default parameters for this execution are:
## -  integrator               = md
## -  dt                       = 0.002
## -  nsteps                   = 5000
## 
## The main output of this BB will be the portable binary run file (.tpr).
#
## In[ ]:
#
#
## Grompp: Creating portable binary run file for mdrun
#from biobb_md.gromacs.grompp import Grompp
## Create prop dict and inputs/outputs
#prop = {'mdp':{'type': 'free', 'nsteps':'15000'}}
#output_gppmd_tpr_path = '1aki_gppmd.tpr'
##Create and launch bb
#Grompp(input_gro_path=output_npt_gro_path, 
#       input_top_zip_path=output_genion_top_zip_path, 
#       output_tpr_path=output_gppmd_tpr_path, 
#       input_cpt_path=output_npt_cpt_path, 
#       properties=prop).launch()
#
#
## ## Execute free molecular dynamics simulation
## Execute the free molecular dynamics simulation  using the MDrun module and the input the portable binary run file (.tpr) as the main input.
## The main output of this BB will be updated structure file (.gro) and the [Gromacs trajectory file](http://manual.gromacs.org/documentation/2019/reference-manual/file-formats.html#trr) (.trr).
#
## In[ ]:
#
#
## Mdrun: Running free dynamics
#from biobb_md.gromacs.mdrun import Mdrun
## Create prop dict and inputs/outputs
#output_md_trr_path = '1aki_md.trr'
#output_md_gro_path = '1aki_md.gro'
#output_md_edr_path = '1aki_md.edr'
#output_md_log_path = '1aki_md.log'
#output_md_cpt_path = '1aki_md.cpt'
##Create and launch bb
#Mdrun(input_tpr_path=output_gppmd_tpr_path, 
#      output_trr_path=output_md_trr_path, 
#      output_gro_path=output_md_gro_path, 
#      output_edr_path=output_md_edr_path, 
#      output_log_path=output_md_log_path, 
#      output_cpt_path=output_md_cpt_path).launch()
#
#
## In[ ]:
#
#
##Show trajectory
#nglview.show_simpletraj(nglview.SimpletrajTrajectory(output_md_trr_path, output_md_gro_path), gui=True)
#
