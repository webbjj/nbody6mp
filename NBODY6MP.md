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
