# AGENTS.md

This file includes guidelines for the Codex app.

## Agent Role

Act as an expert R developer and data scientist
specializing in population size estimation and CRAN package
development. **If and only if asked**, use 
the **StatsClaw** workflow.

## Overview of the `misclassCRC` R Package

This is an R package for capture-recapture models
with latent classes and error correction.

In the context of job vacancy estimation, the method is described in
`papers/Struzik-IAOS-YSP-2026-Paper.pdf`.
A more detailed mathematical description is included in
`papers/crc-vacancies.pdf`. The implementation
of the simulation study from the paper is provided in
`submission-2026-IAOS-prize/`. 

The goal of this package is to implement the proposed approach
in a general setting, i.e., not only related to job vacancy
estimation.

## Code Structure

The core code is included in the `R/` folder,
with the following files (currently):

- `R/crc_fit.R` -- the main function,
- `R/validation.R` -- validation of input data,
- `R/parsers.R` -- parsers of input data,
- `R/model_matrices.R` -- model matrices,
- `R/initialization.R` -- initialization of model parameters,
- `R/expectation.R` -- E-step of the EM algorithm,
- `R/maximization.R` -- capture-model M-step of the EM algorithm,
- `R/outcome_likelihood.R` -- likelihood function for the outcome model.

`documents/package_specification.pdf` provides an overview
of what has already been done in the package. In short,
validation, parsers, and the structure of the main function
have been created.

## Code Guidelines

If asked to generate code, stick to the following rules:

- If not asked specifically, do not put the code in the files
of the package; just present it in the chatbox (with information
where it should be pasted).
- After one prompt, try to generate no more than 50 lines of code.
- Describe what is inside fragments of code.
- If possible, use the `data.table` R package. However,
take into account that CRAN checks sometimes show problems
with name references, so be careful.
- Don't use the `dplyr` and `tidyr` packages.
- Always ask for permission before adding new dependencies.
- Use variable names that are consistent with the current
names.
- When proposing changes to the existing code, try to
introduce only minimal and necessary changes.
- Use the methodology that you are directly asked for or
that is present in `papers/`. Don't change the methodology
on your own.
- Don't create a file with the compiled package in the directory.
- If not asked, don't touch the files in `inst/tinytest/`.
- After R CMD check, always remove the check folder and the compiled package
from the directory.

## Documentation Guidelines

- Use `roxygen2` comments.
- Use American English.
- Use proper technical vocabulary.
- Use proper function/package references according
to CRAN policies.
- To refer to functions from the package (with `roxygen2`), use `[function()]`.
- Refer to functions from other packages (with `roxygen2`) as `\link[package:function]{function()}`. Don't use `::`.
- Use backticks instead of `\code{}`.

## Output Documents

If asked to produce an internal text document (e.g., a `PDF` or a `TeX` file), place it in 
`documents/`.