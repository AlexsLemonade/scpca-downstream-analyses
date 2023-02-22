---
name: Release checklist
about: Prepare for a new release version of scpca-downstream-analyses
title: Prepare for scpca-downstream-analyses release vX.X.X
labels: release

---

## Steps for a new release of `scpca-downstream-analyses`

### Preparing for the release

- [ ] Are all of the issues planned for this release resolved? If there are any issues that are unresolved, mark this issue as blocked by those on ZenHub.
- [ ] Test that the workflow(s) is in good working order and ensure that you have run the latest version of the workflow(s):
  - If making additions to the core workflow, test the workflow using the example data with `snakemake --cores 2 --use-conda`
- [ ] Copy the `example-results` folder output by the workflow to the `/shared/data/scpca-downstream-analyses` directory on the Rstudio server.
- [ ] Log into the `scpca-tester` Rstudio server and do the following:
  - [ ] Ensure the conda environment is up to date, using `bash setup_envs.sh`.
  - [ ] Run the `utils/sync-results.sh` script from within `~/scpca-downstream-analyses` to update the example results files.
- [ ] File a PR from the `development` branch to the `main` branch. This should include all of the changes that will be associated with the next release.

### Creating a release
- [ ] On the [releases page](https://github.com/AlexsLemonade/scpca-downstream-analyses/releases), choose `Draft a new release`.
- [ ] In `Choose a tag`, type a new release number using semantic versioning (vX.X.X) (you did update the title of this issue to match, right?), then click `Create a new tag: vX.X.X on publish`.
- [ ] Write a description of the major changes in this release. You may want to start with the auto-generated release notes to save time.
- [ ] **Optional**: If not all issues have been addressed, save a draft to return to later.
- [ ] Publish the release!
