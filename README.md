# PySTELLA

IMPORTANT: the package is under heavy development. 
Keep this in mind if you want to use it on your research.

## Command shell

  Run the script "PATH-to-PySTELLA/pystella.py" in the directory with models.
  See help if needed.


## Load the flux from ph-files  and convert them to the bands

#### Usage:   ubv.py [params]

-  -b <bands:shift>: string, default: U-B-V-R-I, for example U-B-V-R-I-u-g-i-r-z-UVW1-UVW2.
     shift: to move lc along y-axe (minus is '_', for example -b R:2-V-I:_4-B:5 
-  -i <model-name OR pattern, like '*R450*'>.  Example: cat_R450_M15_Ni007_E7
-  -p <model directory>, default: ./
-  -e <extinction, E(B-V)> is used to define A_nu, default: 0 
-  -c <callback> [lcobs:fname:marker:dt:dm, velobs:fname:marker:dt:vnorm(1e8), popov[:R:M:E[FOE]:Mni], lcobssm as lcobs, but for sm-format data-files]. You can add parameters in format func:params
-  -d <distance> [pc].  Default: 10 pc
-  -g <single, grid, gridm, gridl> Select plot view.  single [default] = all models in one figure, grid = for each band separate figure.
-  -o <is_axes_right>.  Default: empty string
-  -m <magnification>.  Default: None, used for grav lens
-  -q  turn off quiet mode: print info and additional plots
-  -t  plot time points
-  -s  <file-name> without extension. Save plot to pdf-file. Default: ubv_<file-name>.pdf
-  -x  <xbeg:xend> - xlim, ex: 0:12. Default: None, used all days.
-  -y  <ybeg:yend> - ylim, ex: 26:21. Default: None, used top-magnitude+-5.
-  -v  <swd OR ttres[ttresold]> - plot model velocities computed from swd OR tt-res files[ttresold for old res format].
-  -w  write magnitudes to out-file. Use '1' for the default name of out-file
-  -z <redshift>.  Default: 0
-  --dt=<t_diff>  time difference between two spectra
-  --curve-old  - use old procedure
-  --curve-tt  - take curves from tt-file: UBVRI+bol
-  -l  write plot label
-  -h  print usage
    
#### Available bands:
  - AtlasC   AtlasO  
  - BesB     BesI     BesR     BesU     BesV    
  - F105W    F125W    F140W    F160W    F435W    F606W    F814W   
  - GaiaG   
  - HSCY     HSCg     HSCi     HSCr     HSCz    
  - KaitB    KaitI    KaitR    KaitU    KaitV    Kepler  
  - LcoH     LcoJ     LcoK     Lum     
  - PS1g     PS1i     PS1r     PS1w     PS1y     PS1z    
  - Sdssg    Sdssi    Sdssr    Sdssu    Sdssz    SwiftB   SwiftU   SwiftV  
  - USNO40g  USNO40i  USNO40r  USNO40u  USNO40z  UVM2     UVW1     UVW2    
  - bol      bolq    
  - massH    massJ    massK   
  - ubvri   

#### Available aliases of bands: 
  - U => BesU       B => BesB       V => BesV       R => BesR       I => BesI    
  - g => Sdssg      r => Sdssr      i => Sdssi      u => Sdssu      z => Sdssz   
  - J => massJ      H => massH      K => massK  

### Run

```shell
>>>./ubv.py  -i cat_R1000_M15_Ni007_E15  -p stella/res/tt  -b U:2-B-V-R:_1 -z 2 -d 16e9
```

Run script for all *.ph-files in the DIR:
```shell
>>> find DIR  -name "*.ph" | sed -r 's/\.ph$//' | while read fn; do ./ubv.py -i $(basename  $fn)  -d $(dirname $fn) -s; done
```
also the same could be done with key '-p' without model's name. 


#### Install:
```bash
>>> git clone https://github.com/baklanovp/pystella.git
>>> cd pystella
```
and work

```shell
>>> ./pystella.py
```
OR
```shell
>>> ./ubv.py [params]
```
OR

```shell
>>> ipython
ipython> import sys
ipython> sys.path.append('path-to-root-pystella')
ipython> import pystella as ps
```

### Requirements

```shell
>>> pip3 install -r requirements.txt
```

OR manually

To plot light curves [ubv] and shock wave details [swd]
```shell
>>> apt-get install python3
>>> apt-get install python3-numpy python3-scipy
>>> apt-get install python3-matplotlib
??? apt-get install python3-tk
```

To fit observations  
```shell
>>> apt-get install -y python3-pip
>>> pip3 install emcee
>>> pip3 install cython
>>> pip3 install gptools
>>> pip3 install corner
```

## Tests

```shell
>>> python3 -m unittest discover ./tests/
```

Testing can take a long time.

Some tests plot the comparison charts, they must be closed to continue the work.
It is not necessary that all tests show OK. Some of them are non-working yet.

### Citing PySTELLA

Please add some form of acknowledgement that you used this code.

### Acknowledgments:
    Some of the passbands was taken from  SNPY (see http://csp.obs.carnegiescience.edu/data/snpy/),
    SnCosmost (see https://github.com/srodney/sncosmost), and http://svo2.cab.inta-csic.es/svo/theory/fps3/index.php.
  