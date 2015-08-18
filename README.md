# pystella

## Flux from ph-files to the bands

#### Install:
```bash
>>> git clone https://github.com/baklanovp/pystella.git

>>> cd pystella

>>> ./ubv.py [params]
```

#### Usage:

  ubv.py [params]

-  -b <bands>: string, default: U-B-V-R-I, for example U-B-V-R-I-u-g-i-r-z-UVW1-UVW2.
     Available: U-B-V-R-I-u-g-i-r-z-UVM2-UVW1-UVW2-U_UVOT-B_UVOT-V_UVOT
-  -i <model name>.  Example: cat_R450_M15_Ni007_E7
-  -d <model directory>, default: ./
-  -e <model extension> is used to define model name, default: tt 
-  -s  silence mode: no info, no plot
-  -w  write magnitudes to file, default 'False'
-  -h  print usage


```bash
>>>./ubv.py  -i cat_R1000_M15_Ni007_E15  -d ~/Sn/Release/seb_git/res/tt  -b U-B-V
```

Run script for all *.ph-files in the DIR:
```bash
>>> find DIR  -name "*.ph" | sed -r 's/\.ph$//' | while read fn; do ./ubv.py -i $(basename  $fn)  -d $(dirname $fn) -s; done
```
also the same could be done with key '-d' without model's name. 


Acknowledgments:
    The  bands data  was taken from  SNPY, see http://csp.obs.carnegiescience.edu/data/snpy/