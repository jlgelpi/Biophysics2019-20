#!/bin/csh
foreach f (pdbs/*pdb) 
echo $f
perl fixchains.pl his.dat $f > ${f:r}_f.pdb
python binding.py ${f:r}_f.pdb  > ${f:r}.bind.log
rm ${f:r}_f.pdb
end
