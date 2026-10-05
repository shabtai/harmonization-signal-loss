# Is data that is not clean or normalized a problem? Or maybe it's the other way round

Oct 5, 2026 · Ran Tene

Pharma and clinical research still spend years building a common data layer: mapping codes, merging categories, standardizing names and resolving differences between sources. Then maintaining it as the sources and standards change.

The assumption is that models cannot analyse the data reliably until this work is complete.

My intensive benchmarks challenge that assumption. I tested current AI agents (Claude Opus 5.5 and Astra 6) on raw data with mixed codes, inconsistent names, unit changes, split identities, conflicting sources and planted defects. I tried to make them fail on those defects. I did not succeed.

Within the tasks I tested, a permanently harmonized copy was not a prerequisite for reliable analysis.

Harmonization and curation do more than consume years of work. They can erase distinctions that carry scientific signal. When analysis accepts only harmonized data, sources that have not been harmonized are excluded before the analysis even starts.

**The loss is not minor.** In my experience, most of the signal in the data never reaches the analysis. It is lost in three ways:

1. **Data I could not curate.** Sources that were too costly or too complex to map are left out.
2. **Data judged not ready.** Sources that are not yet harmonized are excluded before the analysis starts.
3. **Data lost in curation.** In the data that is curated, merges and mappings erase distinctions that carry signal.

Section 1 shows one published case of the third loss: a single merge, with no failed step and no warning.

The note has two parts:

1. A recent published example in which one harmonization step changed a main result.
2. Results from my tests of current models on many kinds of broken, unharmonized data.

The conclusion describes what this means for data work and pipelines.

## 1. Case study: one merge in a published survival analysis

**Study.** Calonaci et al. (Nature Genetics, 2026) classify mutations by mutant dosage and test the link with survival in MSK-MET (Nguyen et al., Cell, 2022). For KRAS in pancreatic adenocarcinoma they report a hazard ratio (HR) of 3.19 for high dosage against KRAS wild type (WT).

**The harmonization step.** MSK-MET records two pancreatic diagnoses as separate OncoTree codes: adenocarcinoma (PAAD) and neuroendocrine tumour (PANET). These are different diseases, with different driver genes and different survival. The data preparation script (`1.1.prepare_msk_met.R`, line 91) merges them:

```r
tumor_type %in% c('PAAD', 'PANET') ~ 'PAAD'
```

The preprint (medRxiv 2024.05.13.24307238) does not mention PANET, neuroendocrine tumours, or this merge. It calls the merged group "pancreatic adenocarcinoma". I have not checked the published Nature version.

**Effect on the comparison groups.** PANET tumours rarely carry KRAS mutations. In the analysed KRAS model (MSK-MET, 1,768 patients), the 144 PANET patients fall mostly into the reference group:

| Group | PAAD | PANET | PANET share |
| --- | --- | --- | --- |
| WT (reference) | 82 | 132 | 62% |
| Low dosage | 533 | 5 | 1% |
| Balanced dosage | 663 | 3 | <1% |
| High dosage | 346 | 4 | 1% |

**Method.** I loaded the authors' saved Cox model and its input. Refitting it reproduced the coefficients, covariance matrix, log-likelihood, N and number of deaths (tolerance 1e-10). I then refit the same model, with the same covariates (age, TMB, sample type, FGA, sex) and Efron ties, on patients whose original code is PAAD. Nothing else was changed.

**Result** (MSK-MET, KRAS, pancreatic):

| Term | Merged PAAD+PANET (published) | PAAD only |
| --- | --- | --- |
| N (deaths) | 1,768 (1,088) | 1,624 (1,035) |
| Low vs WT | 2.05 (1.57–2.67), p 1.5e-7 | 1.19 (0.88–1.62), p 0.26 |
| Balanced vs WT | 2.29 (1.77–2.97), p 4.5e-10 | 1.32 (0.98–1.79), p 0.069 |
| High vs WT | 3.19 (2.45–4.17), p 1.4e-17 | 1.73 (1.26–2.37), p 0.0006 |
| High vs balanced | 1.39, p 4.9e-5 | 1.31, p 0.0013 |
| FGA (covariate) | 0.91 (0.64–1.29), p 0.59 | 2.46 (1.57–3.85), p 8e-5 |

HR (95% CI), Cox model, Efron ties.

The log-HR against WT falls by 53% (high), 66% (balanced) and 76% (low). High dosage remains associated with worse survival, and the high-vs-balanced contrast holds. Low and balanced dosage are no longer significant at p < 0.05. The FGA covariate goes from no association to HR 2.46, so the merge also distorted the covariate estimates.

**Is it only the smaller sample?** No. I removed 144 random samples (the same number as PANET) from the merged data, 300 times. The high-dosage HR stayed at 3.20 (90% range 3.02–3.43), and every term stayed significant in all 300 runs. The drop to 1.73 comes from which samples were removed, not from how many.

**Other merges.** The script merges 32 OncoTree codes into 9 groups. I refit all 227 published MSK-MET 2021 models in merged groups; 220 reproduce exactly. Three more groups show the same pattern:

| Group (MSK-MET 2021) | Model | Published | Main code only |
| --- | --- | --- | --- |
| Ovarian, 4 histologies | TP53, high dosage | 1.22, p 0.29 | High-grade serous: 2.66, p 0.032 |
| Melanoma, 9 codes | NRAS, balanced | 1.30, p 0.073 | Skin melanoma: 1.85, p 0.00053 |
| Breast, ductal + lobular | CDH1, balanced | 1.37, p 0.042 | Not testable: 91% of mutants are lobular |

In the ovarian group, 81% of TP53 wild-type samples are not high-grade serous, so the published null result compares histologies. For ovarian and melanoma, removing the same number of random samples did not reproduce the change. The other merges (colorectal, bladder, lung, endometrial) change results little. Full details and code: [INCOMMON tumour-type merges: a detailed check](https://shabtai.github.io/incommon-merge-check/).

**Limits of this check.** It changes only the analysed cohort. The dosage priors in INCOMMON are tumour-type specific and were built on the merged group; PANET adds 12 of about 1,560 KRAS-mutant samples (<1%), and I did not rebuild the priors. P-values are per term, without the paper's FDR correction.

**What this shows.** The original code was still in the raw files. No step failed. After the merge, the distinction was gone from every later step. The error is large, and it goes both ways. In pancreatic KRAS the merge doubled the effect on the log scale. In ovarian TP53 it removed about 80% of the effect, and in melanoma TP53 88% (MSK-MET 2021, refits on the main OncoTree code). Removing the same number of random samples did not do this, so the loss comes from the merge, not from sample size.

## 2. Current models on broken data

Harmonize-first made sense when fixed code read the data, because a script cannot decide what an unfamiliar code means. Early language models did not remove this need. On tabular data they made many errors, and these errors were silent and hard to find.

This has changed. In the benchmarks described above, the few errors left were small, and a set of semantic preprocessing rules can reduce them further. The main point is that current models read messy data correctly without a harmonized copy.

So the main reason for harmonize-first no longer holds in the same way. The cost shown in Section 1 remains.

## Conclusion

A harmonization step is a decision about meaning. When it is stored in the data, every later analysis inherits it without seeing it. In the INCOMMON KRAS pancreatic analysis (MSK-MET), one such decision changed the high-dosage HR from 1.73 to 3.19. In other groups of the same study it removed up to 88% of real effects. In my tests, current models handled broken, unharmonized data without failing on the defects. The risk of losing signal through harmonization is real; the need for it is smaller than it was.

This points to a different future for data work and pipelines:

- **Raw data as the stored asset.** Original values are kept and treated as read-only. Curation effort moves from transforming data to documenting where it came from and what its values mean.
- **Merges at question time.** Groupings, mappings and cohort definitions are built by the agent for each analysis, from the scientific question, and kept as code with the result.
- **Decisions made visible.** Each result carries the merges it used, so a reviewer can see them, test them, and rerun the analysis with another choice.
- **Pipelines as short-lived steps.** Instead of one fixed pipeline that feeds every analysis, each question gets its own preparation, built when needed and discarded after.

In this model, the question is no longer how to clean the data once for everyone. It is how to keep the data complete, and make each decision about it explicit.

## References

1. Calonaci N. et al. Gene mutant dosage is associated with prognosis and metastatic tropism in 60,000 clinical cancer samples. Nature Genetics (2026). Preprint: medRxiv 2024.05.13.24307238.
2. Nguyen B. et al. Genomic characterization of metastatic patterns from prospective clinical sequencing of 25,000 patients. Cell (2022). Data: MSK-MET 2021, cBioPortal.
3. Chen et al. npj Precision Oncology 10:34 (2026). doi 10.1038/s41698-025-01227-7.
4. INCOMMON analysis code: `scripts/1.1.prepare_msk_met.R`, line 91; `scripts/4.survival_analysis.R`.
5. Reproduction code and outputs: [github.com/shabtai/incommon-merge-check](https://github.com/shabtai/incommon-merge-check).
