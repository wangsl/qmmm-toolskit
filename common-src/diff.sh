#!/bin/sh

src=$*

echo $src

vmd_plugin=/home/sw77/QChem4.0-Tinker4.2/qmmm-toolskit/vmd-plugin
tnc2_traj=/home/sw77/QChem4.0-Tinker4.2/qmmm-toolskit/tnc2-traj
qchem=/home/sw77/QChem4.0-Tinker4.2/qchem4.0.1-tinker4.2/tinkerman

vmd_src=$vmd_plugin/$*
tnc2_src=$tnc2_traj/$*

OK=1
if [ ! -e $vmd_src ]; then
    echo "$vmd_src does not exit"
    OK=0
fi

if [ ! -e $tnc2_src ]; then
    echo "$tnc2_src does not exit"
    OK=0
fi

if [ $OK -eq 1 ]; then
    diff $vmd_src $tnc2_src
else
    echo -n "diff "
    if [ -e $vmd_src ]; then
	echo -n "$vmd_src"
    fi

    if [ -e $tnc2_src ]; then
	echo -n "$tnc2_src"
    fi
fi


echo " $qchem/$*"


