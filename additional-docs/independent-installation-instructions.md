# Independent installation instructions

If you would like to perform package and dependency installation without the conda environments as described in the main `README.md` file [here](./README.md##snakemakeconda-installation), you can do so after confirming that you have R version 4.1 (the Intel version if you are on a Mac) installed, you will want to make sure that all of the R packages are installed as well.
First install the `renv` package by your preferred method.
Then, from within the `scpca-downstream-analyses` directory, run the following command to install all of the additional required packages:

```
Rscript -e "renv::restore()"
```

Note that pandoc must also be installed and in your path to successfully run the `Snakefile`.
You can install pandoc system-wide by following [pandoc's instructions](https://pandoc.org/installing.html), or you can add it to your conda environment with `mamba install pandoc`.

##### Apple Silicon installations

If you are on an Apple Silicon (M1/M2/Arm) Mac and are not using Snakemake and `conda` to handle dependencies, you will need to be sure that you have the Intel version of R, as Bioconductor packages do not currently support the Arm architecture.
Clicking [this link](https://cran.r-project.org/bin/macosx/base/R-4.2.1.pkg) will download the Intel version of R, version 4.2.1, and you can install R by following the installation instructions.
You will also need to install `gfortan`, a Fortran compiler, to facilitate building certain R packages.
Clicking [this link](https://mac.r-project.org/tools/gfortran-8.2-Mojave.dmg) will download the `gfortran` compiler, and again follow the installation instructions to install it.

If you experience library-related errors that indicate R can't find the Fortran compiler while setting up `renv`, you will want to create the file and/or add the following lines to the file `~/.R/Makevars`:

```
FC  = /usr/local/gfortran/bin/gfortran
F77 = /usr/local/gfortran/bin/gfortran
FLIBS = -L/usr/local/gfortran/lib/gcc
```

These lines will ensure that R can find the newly-installed `gfortran` compiler.
If you need to take these steps, you may need to restart R/terminal to proceed with your setup.

