Welcome to NBODY6MP, a modified version of Sverre Aarseth's NBODY6 and Florent Renaud's NBODY6TT. The purpose of NBODY6MP is to be able to dynamically model the evolution of a star cluster consisting of two stellar sub-population's that have different helium abundances. Moving forward, I am assuming you (the user) have experience using and installing NBODY6. In necessary, Florent Renaud's README.md (also including here) shows how to install this version.

Modifications have been made throughout the code, marked by comments of "BEGIN NBODY6MP" and "END NBODY6MP". The modifcations yield the following changes to the standard version of NBODY6:
- Option to initialze NBODY6MP mode through KZ(50), where the actual value corresponds to the number of sub-populations the user wishes to model. When this is >1, the number of stars in each sub-population and their respective values of Y and Z need to be included in the input file (see /Ncode/define.f for format)
- NBODY6MP mode also requires the initial conditions are being read-in, since it needs to be clear which sub-population each star belongs to. A specific option has been set for KZ(22), specifically KZ(22)=5, which reads in an additional column from fort.10. The new column is for the PTYPE parameter.
- The new variable PTYPE (population type), tags stars based on their population number. Hence modifications can be found throughout the code that alters the PTYPE array whenever the order of stars within the array is altered.
- When the values of Y and Z are read-in for each population, different sets of stellar evolution parameters (ZPARS) are created for each PTYPE. Hence whenever ZPARS is used by a subroutine, an if statement has been added to load the correct ZPARS for the star in question. 
-ZPARS(12) still maintains the He abundance a star would have given the user defined Z. This value is needed for calculating all stellar parameters EXCEPT main sequence lifetime. The user input value of Y is stored in ZPARS(15). When ZPARS(15) is not zero, the main sequence lifetime of the star is scaled via the relation found in Fare, Webb, & Sill 2018 (submitted). 
-Output files have been customized the output ptype as well, so the sub-populations can be easily separated in post-processing. 



Edited the following files in /Ncode

params.h - Add NPOPSMAX which sets the maximum number of sub-populations

coal.f - Coalesence of Roche/CE binary. Give PYTPE of most massive star

common.h - Add arrays for population number (NPOPS), helium content (YPOP), metalicity (ZPOP), ZPARS for a given populations ZPARSP and population type (PTYPE)

cmbody.f - Formation of c.m. body by collision. Give body PTYPE of most massive star 

comenv.f - In case of coalesence use PTYPE of more massive star

data.f - If KZ(50) .GT. 1 then read in number of stars in each population as well as the Y and Z of each sub population
       - If KZ(22).EQ.5 then fort.10 has an extra column to read in PTYPE 

events.f - ??????

expel.f - edit CALL comenv to include star IDs

hrplot.f - Output PTYPE to fort.82 and fort.83

input.f - Set default for KZ(50) and output number of populations if not equal to 0

ksreg.f - Need temporary save variables for PTYPE

ksterm.f - Need temporary save variables for PTYPE

mix.f - Give merger remnant PTYPE of higher mass star

mydump.f - New common block for NBODY6MP treatment

remove.f - Need to adjust PTYPE array when star is removed

roche.f - When necessary, assumed ptype of higher mass star

tail0.f - Copy escaper info

In GPU2:

kspreg.f - save PTYPE
swap.f - keep track of PTYPE

Need to update ZPARS before calling star or hrdiag in the following:
brake.f, expel.f, comenv.f, instar.f, mdot.f, newtev.f, synch.f, roche.f

Make same changes to similar files in ARint, ARchain, Chain, Nchain, and Block
