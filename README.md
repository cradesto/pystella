# pystella

## Flux from ph-files to the bands

#### Usage:

  ubv.py [params]

-  -b "bands": string, default: U-B-V-R-I, for example U-B-V-R-I-u-g-i-r-z-UVW1-UVW2.
    shift: to move lc along y-axe (minus is '_', for example -b R:2-V-I:_4-B:5  
     Available bands: B-F105W-F125W-F140W-F160W-F435W-F606W-F814W-H-HSCg-HSCi-HSCr-HSCy-HSCz-I-J-K-
                      -Kepler-PS1g-PS1i-PS1r-PS1z-R-U-UVM2-UVOTB-UVOTU-UVOTV-UVW1-UVW2-
                      -V-bol-g-i-r-u-w-y-z
-  -i "model name".  Example: cat_R450_M15_Ni007_E7
-  -p "model directory", default: ./
-  -e "extinction, E(B-V)" is used to define A_nu, default: 0
-  -c "callback function" [plot_tolstov, popov[:R:M:E[FOE]:Mni]]. You can add parameters in format func:params
-  -d "distance" [pc].  Default: 10 pc
-  -q  turn off quiet mode: print info and additional plots
-  -s  "file-name" without extension. Save plot to pdf-file. Default: ubv_<file-name>.pdf
-  -t  plot time points
-  -x  "xbeg:xend" - xlim, ex: 0:12. Default: None, used all days.
-  -y  "ybeg:yend" - ylim, ex: 26:21. Default: None, used top-magnitude+-5. 
-  -w  write magnitudes to file, default 'False'
-  -z "redshift".  Default: 0
-  -l  write plot label
-  -h  print usage


```bash
>>>./ubv.py  -i cat_R1000_M15_Ni007_E15  -d stella/res/tt  -b U:2-B-V-R:_1 -z 2 -d 16e9
```

Run script for all *.ph-files in the DIR:
```bash
>>> find DIR  -name "*.ph" | sed -r 's/\.ph$//' | while read fn; do ./ubv.py -i $(basename  $fn)  -d $(dirname $fn) -s; done
```
also the same could be done with key '-p' without model's name. 


#### Install:
```bash

>>> git clone https://github.com/baklanovp/pystella.git

>>> cd pystella

>>> ./ubv.py [params]
```

### Libraries

To fit observations using Gaussian process 
```bash
>>> sudo apt-get install -y python3-pip
>>> pip3 install numpy
>>> pip install cython
>>> pip install gptools
>>> pip install corner
```

Acknowledgments:
    Some of the passbands was taken from  SNPY (see http://csp.obs.carnegiescience.edu/data/snpy/) and
    SnCosmost (see https://github.com/srodney/sncosmost)