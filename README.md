# BIOSTAT 707: ICU Mortality Prediction

**Author:** Haiyi Xiong  
**Dataset:** PhysioNet/Computing in Cardiology Challenge 2012  
**Outcome:** In-hospital death  
**Predictor window:** The first 48 hours of an intensive care unit (ICU) stay

## Project overview

This project develops and evaluates models for in-hospital mortality using early ICU measurements. The work is organized into three checkpoints: characterize and prepare the data, compare models using training-set cross-validation, and evaluate the selected model on held-out data.

| Stage | Data used | Main purpose | Main source | Report or results |
|---|---|---|---|---|
| Checkpoint 1 | Set-a | Cohort description, exploratory data analysis, data quality and missingness | `checkpoint1.R` | Tables and figures in `output/` |
| Checkpoint 2 | Set-a | Predictor construction, preprocessing, cross-validation, model comparison and final model selection | `checkpoint2.qmd` | [Checkpoint 2 report](output/checkpoint2.html) |
| Checkpoint 3 | Set-b and set-c; set-a for benchmark mappings and explanation backgrounds | Evaluation of the saved final model without retraining or tuning on test data | `checkpoint3.qmd` | [Checkpoint 3 report](output/checkpoint3.html) |

Set-a contains 4,000 stays and is used for model development. Set-b and set-c each contain 4,000 stays and are used for evaluation. The outcome is joined by `RecordID`; it is not inferred from measurement values or row order. The evaluation sets share the historical data source and do not establish transportability to a different hospital or current clinical practice.

## Data and repository organization

Data source: [PhysioNet Challenge 2012](https://physionet.org/content/challenge-2012/1.0.0/).

Raw data are kept outside version control. Each set directory must contain its 4,000 extracted individual `.txt` records, rather than only a compressed archive.

```text
data/
  set-a/
  set-b/
  set-c/
  Outcomes-a.txt
  Outcomes-b.txt
  Outcomes-c.txt

checkpoint1.R
checkpoint2.qmd
checkpoint3.qmd
checkpoint3_features.R
checkpoint3_analysis.R
report.scss
pixi.toml
pixi.lock
scripts/
  quarto_windows.bat
output/
  checkpoint2.html
  checkpoint3.html
  ... analysis tables, figures and saved models
```

The outcome files include in-hospital death and the clinical severity scores used as benchmarks. Length of stay and survival are not model predictors. Each record is an ICU stay; distinct record identifiers do not necessarily establish distinct people if repeat admissions cannot be linked.

## Checkpoint 1: Cohort characterization and exploratory analysis

Checkpoint 1 describes set-a before modeling. It reads the individual patient records, extracts baseline descriptors and organizes repeated measurements into long and wide tables. It then summarizes the cohort and examines missingness and unusual measurements.

The analysis includes:

- Baseline descriptors, including age, recorded gender, height, weight and ICU type.
- A long table of repeated measurements, retaining measurement names and times.
- A wide table with first, last, minimum, maximum and mean measurement summaries.
- A descriptive cohort table and an outcome summary.
- A missingness map and a comparison of mortality between stays with and without recorded measurements.
- Exploratory identification of unusual measurements. Flagging an unusual value does not by itself establish that it is erroneous.

| Output | Purpose |
|---|---|
| `output/set-a_long.csv` | Long-format time-varying measurements |
| `output/set-a_wide.csv` | Stay-level measurement summaries for exploratory analysis |
| `output/table1.csv` | Descriptive cohort characteristics |
| `output/outcomes.csv` | Outcome summary |
| `output/missingness_map.png` | Visualization of measurement missingness |
| `output/missingness_vs_death.csv` | Measurement availability and observed mortality |

Missingness may reflect clinical concern, test-ordering practices and monitoring intensity. It should not automatically be treated as random. Checkpoint 1 is exploratory; it does not select the final prediction model or evaluate held-out performance.

## Checkpoint 2: Model development and comparison

Checkpoint 2 rebuilds the modeling table directly from raw set-a, so it does not depend on previously generated Checkpoint 1 CSV files. Each row represents one stay. The target is `In-hospital_death`, coded as death versus survival.

### Predictors and preprocessing

The candidate predictors include baseline age, gender, height, weight and ICU type. For each time-varying measurement, the code derives first, last, minimum, maximum, mean and valid-observation count within hours 0–48. This gives 221 candidate predictor columns before model-specific exclusions.

Compared with the exploratory Checkpoint 1 table, the modeling code treats the sentinel value `-1` as missing throughout measurements and excludes nonpositive baseline height and weight. Measurement counts are zero when no valid observation is available. Other extreme values are retained rather than removed using new, unverified clinical cutoffs.

Preprocessing is learned within the applicable training split. It includes model-specific missingness and zero-variance exclusions, median imputation for numerical predictors, mode imputation for categorical predictors, and encoding or scaling where required. The general sparse-predictor exclusion rule is missingness greater than 98%. Ridge predictors are standardized; the forest does not require standardization. The generalized additive model uses its own selected predictors and preprocessing, including `log1p` for count predictors.

`RecordID` is an identifier, not a predictor. Death, length of stay, survival and benchmark scores are excluded from the three candidate models' input features.

### Models and validation design

| Method | Role and implementation |
|---|---|
| Ridge logistic regression | L2-penalized logistic regression; handles correlated summaries by shrinking coefficients. The penalty is tuned using inner cross-validation. |
| Generalized additive model (GAM) | Logistic model with smooth terms for age, last blood urea nitrogen (BUN), last glucose, maximum heart rate and mean temperature, plus selected additional numerical predictors and categorical gender/ICU type. Additional predictor selection occurs within each outer training fold. |
| Random forest | A 300-tree ensemble capturing nonlinearities and interactions. Candidate predictors per split (`mtry`) and minimum node size (`min_n`) are tuned within training data. |
| Simplified Acute Physiology Score I (SAPS-I) | Clinical benchmark, mapped to mortality probability using a logistic model fitted within each training fold. |
| Sequential Organ Failure Assessment (SOFA) | Clinical benchmark evaluated with the same training-fold mapping approach. |

All methods use the same stratified five outer folds. Ridge and random forest use stratified three-fold inner cross-validation for hyperparameter tuning by log loss. Held-out outer-fold predictions are combined as out-of-fold (OOF) predictions. Performance is summarized using area under the receiver operating characteristic curve (AUROC) and the Brier score, the mean squared probability error.

The final candidate is selected by highest mean outer-fold AUROC among ridge, GAM and random forest, with mean Brier score breaking ties. The selected model is then fitted on all set-a and saved with its preprocessing. Choosing the winner using CV results can introduce selection optimism.

### Results and outputs

The saved selected model is the **random forest**. In the supplied Checkpoint 2 results, its mean five-fold AUROC is **0.8409**, the fold standard deviation is **0.0092**, and the mean Brier score is **0.0927**. Fold variability is descriptive; it is not a confidence interval for the CV mean because training sets overlap.

| Output | Purpose |
|---|---|
| `output/checkpoint2.html` | Full model-development report, with explanations, tables and figures |
| `output/cp2_modeling_table.csv` | Modeling table rebuilt from raw set-a |
| `output/cp2_benchmarks.csv` | Outcomes and clinical benchmark scores |
| `output/cp2_folds.rds` | Saved outer cross-validation splits |
| `output/cp2_oof_predictions.csv` | Aligned OOF predictions for all five methods |
| `output/cp2_all_fold_performance.csv` | Fold-level AUROC and Brier scores |
| `output/cp2_model_comparison.csv` | Model-comparison summary and fold AUROCs |
| `output/cp2_dropped_predictors.csv` | Audit of preprocessing exclusions |
| `output/cp2_final_model.rds` | Selected fitted model and preprocessing for Checkpoint 3 |
| `output/cp2_final_model_selection.csv` | Selected model name and selection criterion |
| `output/cp2_roc_all_models.png` | Receiver operating characteristic (ROC) comparison |
| `output/cp2_ridge_coefficients_top25.png` | Largest standardized ridge coefficients |
| `output/cp2_gam_or_Age.png` | GAM age effect, reused for comparison in Checkpoint 3 |
| `output/cp2_gam_or_HR_max.png`, `output/cp2_gam_or_Temp_mean.png` | Additional GAM effect plots |
| `output/cp2_rf_variable_importance.png` | Random-forest permutation importance |

Some intermediate files in this table are generated by running Checkpoint 2 rather than committed to the repository. The final model, model-selection file, comparison summary and GAM age figure are the existing artifacts needed by Checkpoint 3.

## Checkpoint 3: Evaluation of the final model

Checkpoint 3 loads the saved random forest and its fitted set-a recipe. It does not retrain, tune or recalibrate that model on set-b or set-c. `checkpoint3_features.R` contains the feature-building code extracted from the existing Checkpoint 2, parameterized by dataset name; the original Checkpoint 2 source can remain unchanged. Test predictors use the same definitions and are ordered according to the saved training recipe.

The evaluation covers:

- **Discrimination:** AUROC with 2,000-resample stratified bootstrap 95% confidence intervals, SAPS-I/SOFA comparisons, comparison with Checkpoint 2 CV mean and standard deviation, and an interval for the set-b versus set-c AUROC difference.
- **Precision and recall:** A set-b precision-recall curve with the observed death rate as its no-skill baseline.
- **Calibration:** Flexible loess calibration curves with 200-resample pointwise bootstrap bands and predicted-risk rugs; calibration-in-the-large, slope, Brier score and scaled Brier score. Calibration regressions are diagnostic and do not update predictions.
- **Clinical usefulness:** Decision curve analysis (DCA) comparing the forest, SOFA, review-all and review-none over 5%–50% risk thresholds. At a fixed 15% review threshold, report sensitivity, positive predictive value (PPV), flags per 100 stays and Challenge 2012 event 1, defined as `min(sensitivity, PPV)`.
- **Subgroups:** Set-b AUROC and calibration-in-the-large by age group, recorded sex and ICU type, with stay/death counts and uncertainty intervals.
- **Interpretation:** Saved permutation importance and an age partial dependence plot compared with the original GAM age-effect curve.
- **Confident errors:** The five highest-risk survivors and five lowest-risk deaths on set-b, with local donor-replacement sensitivity, missingness/range checks and case-specific judgments.
- **Conclusions:** Clinical usefulness, single-source and 48-hour-window limitations, leakage risks, and why predictors should not be interpreted causally.

SAPS-I and SOFA probability mappings are fitted only to set-a. Set-a also supplies the explanation background. These uses do not refit the final forest. Local donor-replacement sensitivities are not SHAP values and do not provide an additive or causal decomposition of risk.

### Evaluation summary

| Measure | Set-b | Set-c |
|---|---:|---:|
| ICU stays | 4,000 | 4,000 |
| Deaths | 568 | 585 |
| Mortality rate | 14.20% | 14.63% |
| Random-forest AUROC | 0.8305 | 0.8356 |
| Bootstrap 95% AUROC CI | 0.8117–0.8474 | 0.8180–0.8528 |
| Calibration intercept | −0.0767 | −0.0394 |
| Calibration slope | 0.8995 | 0.9967 |
| Brier score | 0.0950 | 0.0956 |

The model provides useful discrimination within this historical source, but these results do not justify immediate clinical deployment. The proposed 15% clinician-review threshold is illustrative and was not optimized on test outcomes. The report discusses alert burden, missed deaths and calibration. Prior exploratory exposure to the test sets is disclosed; this is not represented as a prospectively registered, first-ever test-set evaluation.

| Output | Purpose |
|---|---|
| `output/checkpoint3.html` | Complete evaluation report with embedded figures |
| `output/cp3_cohorts.csv`, `output/cp3_predictions.csv` | Cohort counts and frozen stay-level predictions |
| `output/cp3_discrimination.csv`, `output/cp3_cv_comparison.csv` | Held-out discrimination, uncertainty and CV comparison |
| `output/cp3_calibration.csv`, `output/cp3_calibration_curve.csv`, `output/cp3_risk_regions.csv` | Calibration metrics, smooth curve and descriptive risk regions |
| `output/cp3_decision_curve.csv`, `output/cp3_net_benefit_15.csv`, `output/cp3_threshold.csv` | Net benefit and threshold-specific operating characteristics |
| `output/cp3_subgroups.csv` | Subgroup sample sizes, discrimination and calibration |
| `output/cp3_importance.csv`, `output/cp3_age_pdp.csv` | Global importance and age partial dependence |
| `output/cp3_local_sensitivity.csv`, `output/cp3_confident_errors.csv` | Local explanations and ten error-case assessments |
| `output/cp3_*.png` | Evaluation figures |
| `output/cp3_runtime.csv`, `output/cp3_session_info.txt` | Measured runtime, model checksum and software versions |

## Running the project

The project uses Pixi to manage R, Quarto and analysis dependencies. `pixi.toml` defines the tasks and `pixi.lock` fixes the dependency environment. Run commands from the repository root after placing raw data in `data/`.

| Checkpoint | Command | Required inputs |
|---|---|---|
| 1 | `pixi run --locked checkpoint1` | Raw set-a and `Outcomes-a.txt` |
| 2 | `pixi run --locked checkpoint2` | Raw set-a and `Outcomes-a.txt` |
| 3 | `pixi run --locked checkpoint3` | All three raw sets/outcomes and the saved Checkpoint 2 artifacts listed below |

Checkpoint 1 and 2 are independently runnable because Checkpoint 2 rebuilds its own modeling table. Checkpoint 3 needs these existing files:

```text
output/cp2_final_model.rds
output/cp2_final_model_selection.csv
output/cp2_model_comparison.csv
output/cp2_gam_or_Age.png
```

When those artifacts are already present, there is no need to rerun Checkpoint 1 or 2 to evaluate Checkpoint 3. Rerunning Checkpoint 2 intentionally regenerates its training outputs; Checkpoint 3 subsequently evaluates the saved model that is present.

The Checkpoint 3 task includes `--no-clean` to preserve existing model artifacts in the shared `output/` directory. The original Checkpoint 1 and 2 task commands are retained in the incremental update. Their original Quarto command may clean its output directory when rerun, so the committed report remains the reference until a new run finishes. On Windows, `scripts/quarto_windows.bat` resolves the Quarto Deno DOM library from the active environment.

All raw inputs are read from `data/`; generated analysis artifacts are written to `output/`. Checkpoint 3 additionally reads the existing Checkpoint 2 artifacts from `output/`. No new dependency or `fastshap` installation is needed for the delivered Checkpoint 3 implementation. HTML reports are self-contained and can be downloaded and opened in a browser; GitHub source preview does not render them as a webpage.

### Observed runtime

| Checkpoint | Runtime | Measurement context |
|---|---:|---|
| 1 | Approximately 0.65 minutes | Verification run on the available Windows/R 4.5.3 environment |
| 2 | Approximately 6.99 minutes | Existing completed run recorded in `output/cp2_runtime.csv` |
| 3 | Approximately 2.3–2.4 minutes | Analysis time recorded in `output/cp3_runtime.csv`; HTML packaging is additional |

These timings exclude first-time environment installation and vary with disk speed and CPU. The development runs used compact tuning grids and completed in less than 30 minutes. A fresh offline environment installation could not finish because a package was absent from the local cache; full fresh-environment installation remains unverified. Report generation was verified using the existing R 4.5.3 Pixi environment.

## Reproducibility and interpretation safeguards

- Random seeds control folds, tuning-related resampling, bootstrap intervals and explanation sampling.
- Checkpoint 2 estimates preprocessing within training splits and uses held-out outer-fold predictions for evaluation.
- Checkpoint 3 checks stay counts, outcome completeness, identifier overlap, predictor availability and valid prediction ranges.
- A checksum verifies that the final model is unchanged during Checkpoint 3. The evaluated model matches GitHub blob SHA `8371d2ba4784865537a4b19823c921c74fe55b8a`.
- `.gitignore` excludes raw data and environment files while allowing the rendered reports and selected evaluation outputs to be committed.
- Performance and explanation results are predictive associations, not estimates of causal treatment effects. Monitoring intensity, sedation and treatment response can influence recorded predictors.
