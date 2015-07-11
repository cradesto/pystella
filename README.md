# pystella

## Flux from ph-files to the bands

#### Install:
```bash
>>> git clone https://github.com/baklanovp/pystella.git

>>> cd pystella

>>> ./sn.py [params]
```

#### Usage:

  sn.py [params]
  
-  -b <bands>: string like U-B-V-R-I-g-r-i-UVM2-UVW1-UVW2, default: U-B-V-R-I
-  -i <model name>.  Ex: cat_R1000_M15_Ni007_E15
-  -d <model directory>, default: ./
-  -s silence mode: no info, no plot
-  -h: print usage

```bash
>>>./sn.py  -i cat_R1000_M15_Ni007_E15  -d ~/Sn/Release/seb_git/res/tt  -b U-B-V
```

Run script for all *.ph-files in the DIR:
```bash
>>> find DIR  -name "*.ph" | sed -r 's/\.ph$//' | while read fn; do ./sn.py -i $(basename  $fn)  -d $(dirname $fn) -s; done
```



Acknowledgments:
    The  bands data  was taken from  SNPY, see http://csp.obs.carnegiescience.edu/data/snpy/