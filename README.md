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

-  -b "bands": string, default: U-B-V-R-I, for example U-B-V-R-I-u-g-i-r-z-UVW1-UVW2.
     Available: B-H-HSCg-HSCi-HSCr-HSCy-HSCz-I-J-K-PS1g-PS1i-PS1r-PS1z-R-U-UVM2-UVOTB-UVOTU-UVOTV-UVW1-UVW2-V-g-i-r-u-w-y-z
-  -i "model name".  Example: cat_R450_M15_Ni007_E7
-  -p "model directory", default: ./
-  -e "model extension" is used to define model name, default: tt 
-  -c "callback function" [plot_tolstov].
-  -d "distance" [pc].  Default: 10 pc
-  -z "redshift".  Default: 0
-  -s  silence mode: no info, no plot
-  -t  plot time points
-  -w  write magnitudes to file, default 'False'
-  -h  print usage


```bash
>>>./ubv.py  -i cat_R1000_M15_Ni007_E15  -d stella/res/tt  -b U-B-V -z 2 -d 16e9
```

Run script for all *.ph-files in the DIR:
```bash
>>> find DIR  -name "*.ph" | sed -r 's/\.ph$//' | while read fn; do ./ubv.py -i $(basename  $fn)  -d $(dirname $fn) -s; done
```
also the same could be done with key '-d' without model's name. 


Acknowledgments:
    Some of the data of bands was taken from  SNPY (see http://csp.obs.carnegiescience.edu/data/snpy/) and
    SnCosmost (see https://github.com/srodney/sncosmost)