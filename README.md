Welcome to NBODY6MP, a modified version of Sverre Aarseth's NBODY6 and Florent Renaud's NBODY6TT. The purpose of NBODY6MP is to be able to dynamically model the evolution of a star cluster consisting of two stellar sub-population's that have different helium abundances. Moving forward, I am assuming you (the user) have experience using and installing NBODY6. In necessary, Florent Renaud's README.md (also including here) shows how to install this version.

Modifications have been made throughout the code, marked by comments of "BEGIN NBODY6MP" and "END NBODY6MP". The modifcations yield the following changes to the standard version of NBODY6:
- Option to initialze NBODY6MP mode through KZ(50), where the actual value corresponds to the number of sub-populations the user wishes to model. When this is >1, the number of stars in each sub-population and their respective values of Y and Z need to be included in the input file (see /Ncode/define.f for format)
- NBODY6MP mode also requires the initial conditions are being read-in via fort.10, since it needs to be clear which sub-population each star belongs to. A specific option has been set for when KZ(22)>=2 and KZ(50)>1, where an additional column is read from fort.10. The new column is for the PTYPE parameter.
- The new variable PTYPE (population type), tags stars based on their population number. Hence modifications can be found throughout the code that alters the PTYPE array whenever the order of stars within the array is altered.
- When the values of Y and Z are read-in for each population, different sets of stellar evolution parameters (ZPARS) are created for each PTYPE. Hence whenever ZPARS is used by a subroutine, an if statement has been added to load the correct ZPARS for the star in question. In some cases, subroutines had to be edited to accept ZPARS and sometimes even location of star in the array since ZPARS was used within.
-ZPARS(12) still maintains the He abundance a star would have given the user defined Z. This value is needed for calculating all stellar parameters EXCEPT main sequence lifetime. The user input value of Y is stored in ZPARS(15). When ZPARS(15) is not zero, the main sequence lifetime of the star is scaled via the relation found in Fare, Webb, & Sill 2018 (submitted). 
-Output files have been customized the output ptype as well, so the sub-populations can be easily separated in post-processing. 

To make life a little easier, a sample input file to model two sub-populations with different Helium abundances would be:

1 1000000.0 <br />
40000 1 100 27774 400 1 <br />
0.02 0.02 0.15 2.0 10.0 12000.0 1.0 6.34505264 0.60884957 <br />
2 0 0 0 1 0 1 0 0 0 <br />
0 3 0 3 2 1 0 1 3 0 <br />
1 4 2 0 0 2 2 1 0 1 <br />
1 0 2 2 1 0 0 2 0 3 <br />
0 0 0 0 0 0 0 0 0 2 <br />
2.0E-06 0.0001 0.2 1.0 1.0E-06 0.001 <br />
2.3 50.0 0.1 0 0 0.0001 0.0 5.0 <br />
20000 0.0 0.0001 <br />
20000 0.3 0.0001 <br />
0.5 0.0 0.0 50.0 1.0 <br />
1.5D+10 5.0D+10 4.0 0.5 220.0 8.5 0.0 0.0 0.0 <br />
20.0 0.0 0.0 0.0 231.8 0.0 <br />
0.0 1.0 1.0 1.0 <br />
1.0 <br />

The new rows required by NBODY6MP are rows 11 and 12, each of which contain the parameters NPOP(I), YPOP(I), and ZPOP(I). NPOP(I) is the number of stars in population the i'th population, YPOP(I) is the He abundance and ZPOP(I) the metallicity. Note when YPOP(I)=0 the default value for Y in NBODY6 is used for the given value of ZPOP(I). The corresponding fort.10 file would containt M,X,Y,Z,VX,VY,VZ,PTYPE for 40,000 stars.

That's about everything. Please contact me with any issues or questions you may have (webb.jjw@gmail.com is likely my most permanent email address at the moment).

***Note on Binaries and Mergers***
I should warn you that the code has primarily been testing using single stars only and no primordial binaries. It "should" be setup that primordial binaries have the same population type, and any newly formed binaries and any merger products retain the population type of the most massive star. Again, this has not been tested thoroughly so include binaries at your own risk. Any bugs please bring to my attention.

********
Below is an "incomplete" list of files that were edited in the making of NBODY6MP. I am working to get this up to date in case anyone wants to make the NBODY6MP changes to their own working version of NBODY6. Again, while this is incomplete, all changes to the code can be found via grep "BEGIN NBODY6MP"

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
