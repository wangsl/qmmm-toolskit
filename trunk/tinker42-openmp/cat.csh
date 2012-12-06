#!/bin/csh -f

foreach f (*.[Fi]) 
    echo $f
    set tmp_name = $f-tmp
    if(-e $tmp_name) rm $tmp_name
    echo " " > $tmp_name
    echo "C     "\$Id\$ >> $tmp_name
    cat $f >> $tmp_name
    mv $tmp_name $f
end

foreach f (*.[hC]) 
    echo $f
    set tmp_name = $f-tmp
    if(-e $tmp_name) rm $tmp_name
    echo " " > $tmp_name
    echo "/* "\$Id\$ "*/" >> $tmp_name
    cat $f >> $tmp_name
    mv $tmp_name $f
end
