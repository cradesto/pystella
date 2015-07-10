# pystella

## Flux from ph-files to the bands

Usage:
  sn.py [params]
  -b <bands>: string like U-B-V-R-I-g-r-i-UVM2-UVW1-UVW2, default: U-B-V-R-I
  -i <model name>.  Ex: cat_R1000_M15_Ni007_E15
  -d <model directory>, default: ./
  -h: print usage

```bash
>>>./sn.py  -i cat_R1000_M15_Ni007_E15  -d ~/Sn/Release/seb_git/res/tt  -b U-B-V
```



Acknowledgments:
    The  bands data  was taken from  SNPY, see http://csp.obs.carnegiescience.edu/data/snpy/