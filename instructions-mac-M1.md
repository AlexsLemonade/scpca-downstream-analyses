This workflow for downstream analysis of ScPCA data uses `renv` to manage R and R package versions.
When working on a _MacOS machine with an M1 (silicon) chip_ ("M1 computer"), additional steps have to be taken to ensure R and `renv` can be properly configured. 
The instructions walk through how to set up a local environment on an M1 computer.


1. First, ensure that command line tools are installed on your computer. 
You can install these by opening Terminal, and running `xcode-select --install`. 
The subsequent installation process will take 10-20 minutes.

2. Install the **intel version** of R, which still works perfectly fine on M1 computers using Rosetta under the hood, specifically R version `4.1.2`. 
This R version is necessary for certain Bioconductor packages to work.
[**Click here**](https://cran.r-project.org/bin/macosx/base/R-4.1.2.pkg) to automatically install the appropriate version of R, and follow all subsequent installation instructions as normal.

3. Install the Desktop version of **RStudio** [here](https://www.rstudio.com/products/rstudio/download/#download), and follow all subsequent installation instructions as normal.

4. Next, you will need to install `gfortran` to be able to build R packages (instructions adapted from [here](https://mac.r-project.org/tools/)):
    * Click [**this link**](https://mac.r-project.org/tools/gfortran-12.0.1-20220312-is-darwin20-arm64.tar.xz), which will download a compressed file `gfortran-12.0.1-20220312-is-darwin20-arm64.tar.xz` to your computer.
    * Open the Terminal, and navigate to your `Downloads/` directory (or to a different location where this file downloaded). Enter the following in Terminal to install `gfortran` into its proper location:
        ```
        # Note because of sudo you'll be prompted for your password
        sudo tar fxz gfortran-12.0.1-20220312-is-darwin20-arm64.tar.xz -C /
        ```
    * You will then need to modify your `PATH` environment variable to use this newly-installed `gfortran`. Add the following line to your profile (e.g. `.bash_profile`, `.bashrc`, `.zshrc`, etc.):
            ```
            export PATH=$PATH:/opt/R/arm64/gfortran/bin
            ```
5. You also need to install `openmp`, as follows, from the Terminal:
    ```
    curl -O https://mac.r-project.org/openmp/openmp-12.0.1-darwin20-Release.tar.gz
    # Again, you may be promptd for your password here:
    sudo tar fvxz openmp-12.0.1-darwin20-Release.tar.gz -C /
    ```

6. You now need to modify (or create if it doesn't exist yet!) the file `~/.R/Makevars` where you will provide information about how R can use these newly-downloaded softwares (`gfortran` and `openmp`) when building R packages.
Add these four lines to the `~/.R/Makevars` file:
    ```
    CPPFLAGS+=-I/usr/local/include -Xclang -fopenmp
    LDFLAGS+=-L/usr/local/lib -lomp

    FC=/opt/R/arm64/gfortran/bin/gfortran -mtune=native
    FLIBS=-L/opt/R/arm64/gfortran/lib/gcc/aarch64-apple-darwin20.6.0/12.0.1 -L/opt/R/arm64/gfortran/lib -lgfortran -lemutls_w -lm
    ```

7. From the local respository directory `scpca-downstream-analysis`, open the RStudio Project file `scpca-downstream-analysis.Rproj` by double-clicking to launch RStudio.

8. Install the `renv` package by executing the following in the R Console within RStudio: `install.packages("renv")`. 
The `renv` package will then install for this version (`4.1.2`) of R.

9. At this point, you should be ready to set up the rest of the project with `renv`. 
In the R Console within Rstudio, run: `renv::restore()` and follow prompts to install all package versions needed for the workflow.
**IMPORTANTLY,** you will have to run `renv::restore()` several times.
This is because your Mac will need to use these newly-downloaded softwares (`gfortran` and `openmp`) as well as certain compiler libraries installed along with them.
Since this will be the first time those softwares/libaries are "opened," you will receive standard Mac prompts, "XYZ can't be opened because Apple cannot check it for malicious software."
When you see these prompt, click "Open," and then also open `System Preferences -> Security and Privacy -> General`," and click "Allow" (you may need to click the unlock icon first) so that R can use it to build packages.
You will need to iterate through this process *several times*: 
    * Run `renv::restore()`.
    * When the system prompt appears, click "Open," and then click "Allow" in `System Preferences`.
    * The `renv::restore()` will then fail, but the _next time_ you run `renv::restore()` the particular software you allowed will work smoothly.
    * Run `renv::restore()` again!
After several rounds, there will be no further system prompts and all R packages should install.




