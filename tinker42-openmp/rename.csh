#!/bin/csh -f

foreach f (*.f) 
    set new_name = $f:r.F
    echo $f
    mv $f $new_name
end
