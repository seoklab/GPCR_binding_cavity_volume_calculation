# GPCR binding cavity volume calculation

## Introduction
 &nbsp;&nbsp;&nbsp;This program calculates the volume of the binding cavity of a given GPCR structure with a grid-based algorithm.

## Pre-requisites
 &nbsp;&nbsp;&nbsp;Make sure you have installed the following prerequisites on your
environment. All prerequisites is publicly available. Mismatch of version might not be a problem, but I have not tested about it.
 * Python (version 3.8.10)
 * Numpy (version 1.21.0)
 * Torch (version 1.10+cu102)

## Input preparation
  &nbsp;&nbsp;&nbsp;A query GPCR structure must be aligned to the z-axis properly. The direction of the extracellular side should be aligned with the negative direction of the z-axis.
 This kind of alignment can be done by ppm3 webserver from OPM database(https://opm.phar.umich.edu/ppm_server3).
 You can also simply superpose a query structure to demo input PDB files already aligned.
 The query structure has to be an all-hydrogen topology structure. 
 Many protein viewer programs support superposition and adding hydrogen, ex UCSF chimera.



## Installation
  &nbsp;&nbsp;&nbsp;No specific installation step is required.

## Usage
  ```
	usage: python run.py [-h] -p PDB_FN -toggle TOGGLE_RESNO -tip_s TIP_RESNO_S
	                     [-trim_Nterm TRIM_NTERM] [-trim_Cterm TRIM_CTERM]
	                     [-exclude_TM1_side_truncation EXCLUDE_TM1]
	
	optional arguments:
	  -p:           
			         A query pdb file. ex) -p demo_5zbq.pdb
	  -toggle:      
			         Residue number of toggle switch residue. ex) -toggle 276
	  -tip_s:       
			         Tip residues of TM1~7.ex) -tip_s 39,103,109,177,205,289,295
	  -trim_Nterm:
			         Residues, which have lower residue number than this variable, 
				are excluded for calculating center of protein. (default:35)
	  -trim_Cterm:
			         Residues, which have lower residue number than this variable, 
				are excluded for calculating center of protein. (default:9999)
	  -exclude_TM1_side_truncation:
	                         Exclude TM1 for side truncation.(default:True)
  ```

  * &nbsp;&nbsp;&nbsp;Additional expalnation for option

 &nbsp;&nbsp;&nbsp;'trim_Nterm' and 'trim_Cterm' option has no significant effect on the results. This option has been added for a fair comparison of results according to the presence or absence of the C-terminal region or N-terminal region.  
  &nbsp;&nbsp;&nbsp;'Exclude TM1' option can be used when there is a volume area between TM1, TM2, and TM7 that is separated from the main binding cavity, and it is reasonable to exclude it from the binding cavity volume. Note that TM1 is excluded only for side truncation. For the current NPY1R target analysis, this option has no significant effect on the results.

## Instruction for reproduce data on paper
 &nbsp;&nbsp;&nbsp;Run the following command 

 For 5zbq volume calculation : `python run.py demo_input/demo_5zbq_trim.pdb -toggle 276 -tip_s 39,103,109,177,205,289,295`

 For 5zbh volume calculation : `python run.py demo_input/demo_5zbh_trim.pdb -toggle 276 -tip_s 39,103,109,177,205,289,295`

 For 7vgx volume calculation : `python run.py demo_input/demo_7vgx.pdb -toggle 276 -tip_s 39,103,109,177,205,289,295`

## Expected run time
 &nbsp;&nbsp;&nbsp;For normal desktop envrionment, it takes about 80 secs.

## Result files description 
 * `*_rec.pdb`: Repositioned protein receptor pdb file.
 * `*_truncated.pdb`: Visualized side and roof truncated grid space.
 * `*_volume.pdb`: Visualized binding cavity grid points. 

## Acknowledgement 
&nbsp;&nbsp;&nbsp;Thanks for Jinsol Yang

## Contact
  &nbsp;&nbsp;&nbsp;Hyeonuk Woo : dngusdnr1@gmail.com
