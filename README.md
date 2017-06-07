# pymc4c
A python based approach to processing MC4C data

### Base calling
Refer to [wiki](https://github.com/UMCUGenetics/pymc4c/wiki/Converting-raw-signals-(i.e.-Squiggle)-to-FAST5) to convert the squiggle signals to FAST5 reads.


### Running pyMC4C

#### 1. Create a config file
You need to have a config file that contains experiment's detail (e.g. primer sequence). You can find an example of such a file in "config_dir" folder.

#### 2. Make a fasta file out of primer sequences
Example:

```
	mc4c.py makeprimerfa ./config_dir/CTp-VpID.cfg ./workspace_dir/Cleave_Reads/PRM_ZZZ-TMP.fasta
```