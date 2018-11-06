# Structure originally present in Harvey's 1998 publication: C6H5+

"The singlet and triplet states of phenyl cation. A hybrid approach for locating minimum energy crossing points
between non-interacting potential energy surfaces"
1998, Theoretical Chemistry Accounts 99(2):95-99
DOI: 10.1007/s002140050309


To re-run the calculation, make sure `g09` is available in `$PATH` and all needed environment variables set. For example:

```
export g09root="$HOME/.local/g09/E6L-103X"
export GAUSS_EXEDIR="$HOME/.local/g09/E6L-103X/g09"
export LD_LIBRARY_PATH="$HOME/.local/g09/E6L-103X/g09:$LD_LIBRARY_PATH"
export TMPDIR="$HOME/tmp/scratch"
export PATH="$PATH:$HOME/.local/g09/E6L-103X/g09"
```

Since the input files are already prepared (`geom`, `Input_Header_A`, `Input_Header_B`, `footer`, `ProgFile`) and the `MECP.x` is properly compiled for 11 atoms (`C6H5+`), all you need to do is:

```
cd calc
rm -rf JOBS Gau* *.log *.gjf *.chk
csh -f Make_First_inputs
csh -f master_script
```
