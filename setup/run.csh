#!/bin/csh -f

# $Id: run.csh 4 2012-12-05 04:19:17Z wangsl2001@gmail.com $ 

# Note: Please use Amber10, 
# ptraj will generate PDB file with continuous atom index 

rm -f *.pdb *.zpdb *.xyz *.exe *.input

if($AMBERHOME:t != "amber10" && $AMBERHOME:t != "amber11" ) then
    echo " Please setenv AMBERHOME and use amber10"
    exit
endif

# Step 1: 
# conver Amber restart file restrt-5000 to PDB file

cat <<EOF | $AMBERHOME/exe/ptraj SET79.prmtop
trajin restrt-5000
trajout SET79 pdb
image center
go
EOF

if(! -e SET79.1) then
    echo SET79.1 file does not exit
    exit
else
    mv SET79.1 SET79-Amber.pdb
endif


# Step 2: 
# Cut sphere: remove water molecules
# Sphere center coordinates is CE (SAM)
# Box dimension: from Amber restart file 

cat <<EOF >! cut-sphere.input
PDB-Input SET79-Amber.pdb
PDB-Output sphere-30.pdb
Center-Coordinate 29.311  29.267  31.941
Sphere-Radius 30.0
Box-Dimension 85.2212166  62.9545826  60.1530603
EOF

perl cut-sphere.pl cut-sphere.input

# Step 3: 
# Convert sphere pdb file to zpdb file

cat <<EOF >! pdb2zpdb.input
PDB-Input sphere-30.pdb
ZPDB-Output SET79-30.zpdb
N-Terminal-Resid 1 251
C-Terminal-Resid 250 260
Amber-Tinker-Table amber2tinker_hmt.table
EOF

perl pdb2zpdb.pl pdb2zpdb.input


# Setp 4:
# Convert zpdb file to Tinker XYZ file

set bond_search_src = bond-search.C
set bond_search_exe = bond-search.exe
if(! -e $bond_search_exe) g++ -DMAXNBONDS=10 -O3 -o $bond_search_exe $bond_search_src

cat <<EOF >! zpdb2xyz.input
Tinker-Parameter amber99_m.prm
ZPDB-Input SET79-30.zpdb
ZPDB-Output SET79-qmmm-34.zpdb
Tinker-XYZ SET79-qmmm-34.xyz
Bond-Search-Exe ./bond-search.exe

QM-Atoms 3863 3864 3865 3866 3867 3868 3869 3870 3871
QM-Atoms 3983 3984 3985 3986 3987 3988 3989 3990 
QM-Atoms 3991 3992 3993 3994 3995 3996 3997 3998 3999
QM-Atoms 4000 4001 4002 4003 4004 4005
Pseudo-Atoms 3860 4006
Zero-Charge-Atoms 3857 3858 3859 3861 3862 4007 4008 4009
EOF

perl zpdb2xyz.pl zpdb2xyz.input

cat <<EOF >! zpdb2xyz.input
Tinker-Parameter amber99_m.prm
ZPDB-Input SET79-30.zpdb
ZPDB-Output SET79-qmmm-66.zpdb
Tinker-XYZ SET79-qmmm-66.xyz
Bond-Search-Exe ./bond-search.exe

QM-Atoms 3983 3984 3985 3986 3987 3988 3989
QM-Atoms 3990 3991 3992 3993 3994 3995 3996 3997 3998 3999 
QM-Atoms 4000 4001 4002 4003 4004 4005 4006 4007 4008 4009 
QM-Atoms 4010 4011 4012 4013 4014 4015 4016 4017 4018 4019 
QM-Atoms 4020 4021 4022 4023 4024 4025 4026 4027 4028 4029 
QM-Atoms 4030 4031 4032
QM-Atoms 3857 3858 3859 3860 3861 3862 3863 3864 3865 3866 
QM-Atoms 3867 3868 3869 3870 3871
Pseudo-Atoms 3855
Zero-Charge-Atoms 3853 3854 3856 3872 3873
EOF

perl zpdb2xyz.pl zpdb2xyz.input
