# Independent installation instructions

If you would like to perform package and dependency installation without the conda environments as described in the main `README.md` file [here](./README.md##snakemakeconda-installation), you can do so after confirming that you have R version 4.2 installed.
Then follow the below instructions to ensure that you have all of the necessary R packages to run the workflow installed as well.
First install the `optparse` and `renv` packages by your preferred method.
Then, from within the `scpca-downstream-analyses` directory, run the following command to install all of the additional required packages:

```
Rscript -e "renv::restore()"
```

Note that pandoc must also be installed and in your path to successfully run the `Snakefile`.
You can install pandoc system-wide by following [pandoc's instructions](https://pandoc.org/installing.html), or you can add it to your conda environment with `mamba install pandoc`.

## Apple Silicon installations

As of Bioconductor 3.16, binary packages of Bioconductor packages are now available for Apple Silicon (M1/M2/Arm) Macs.
These binary packages should make installation using native code without Snakemake/conda possible, but we have not extensively tested this.

[This link](https://cran.rstudio.com/bin/macosx/big-sur-arm64/base/R-4.2.2-arm64.pkg) will download the Arm64 version of R, version 4.2.2, and you can install R by following the installation instructions.

You may also need to install `gfortan`, a Fortran compiler, to facilitate building certain R packages.
To do so, follow the instructions at https://mac.r-project.org/tools/ to install the Apple Silicon (arm64) version of `gfortran`.
Briefly, you will download the latest version of the `gfortran` package, currently at https://mac.r-project.org/tools/gfortran-12.0.1-20220312-is-darwin20-arm64.tar.xz
After download, unpack the file with:

```
sudo tar fxz gfortran-12.0.1-20220312-is-darwin20-arm64.tar.xz -C /
```

After initial installation, update the link to the macOS SDK with:

```
sudo /opt/R/arm64/gfortran/bin/gfortran-update-sdk`
```

You will also want to add `/opt/R/arm64/bin` to your `PATH`  by modifying `~/.zshrc` or `~/.bashrc` (depending on your shell) with:

```
export PATH=$PATH:/opt/R/arm64/gfortran/bin
```
