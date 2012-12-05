#!/bin/csh -f

set pwd = `pwd`

cd Save

set src = *.[hC]

cd $pwd

echo $src

set Tinker_dir = $QCSOURCE/tinkerman

echo $Tinker_dir

foreach src ($src)
    set target_file = $Tinker_dir/$src
    if(-e $target_file) then
	ln -s $target_file .
    else
	echo $src does not exist
    endif
end
    
