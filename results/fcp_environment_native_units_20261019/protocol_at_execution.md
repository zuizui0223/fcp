# FCP environmental transport and common physical units: additional-analysis protocol

Fixed on 2026-10-19 before executing the new comparisons. The existing measured roster, colour classification, predictor choice and 1,574-species panel have already been examined. This is a post hoc internal validation of an existing research programme, not a newly preregistered or external data validation.

## Scope and inputs

Use the authentic 447,500-photo CSV (SHA256 `8e7f10f0427927208485b2da6fa9c46bb681d6f3611e53b8dfdc606cfc8dabe6`) and frozen WorldClim rasters. Preserve the seven-variable complete-case selection, at least 40 classified photos and at least five white and five chromatic photos. The expected panel is 1,574 species / 86,114 photos. Use the existing BIO5, BIO14 and SRADmean predictors. All existing outputs remain intact.

## A. Common physical units

Transform each existing primary joint slope and covariance to increments of 1 degree C BIO5, 10 mm BIO14 and 1,000 kJ m^-2 day^-1 SRADmean. If the existing species SD is s and the physical increment is d, the new slope is b*d/s, its variance v*(d/s)^2 and cross-covariance V_jk*d_j*d_k/(s_j*s_k). This is exact re-expression of fitted models, not new photo fits.

Report equal-species means and the same normal-REML / modified-HKSJ summaries. Compare GIFT self-fertilization present-minus-absent, annual-minus-perennial and woody-minus-non-woody using the same admitted categories and species, with 5,000 shared-genus bootstrap draws per trait. Keep all nine comparisons and their nominal and approximate family9 intervals. Do not equate different numerical results with contradiction: weighting and the estimand change with units. Record small environmental ranges and uncertain estimates; no response-driven trimming.

## B. Transfer to held-out genera

Assign entire genera (first species-name word) deterministically to 10 folds using SHA256 of a documented seed and genus string. A held-out species and every species of its genus must be absent from coefficient training.

Learn two fixed common-score models from the other genera: (1) the BIO5 component of the existing joint fits, and (2) all three joint components. Aggregate training-species slopes in common physical units by marginal REML; save all training fold coefficients. Do not use held-out species' fitted coefficients, labels, or colour prevalence to choose the transferred score. No tuning or selection of the winning model.

Evaluate the transferred score within each held-out species. Primary metric: mean pairwise conditional log-loss gain over log(2), using every white-versus-chromatic photo pair in that species and then equal species weight. For score difference delta = score(white)-score(chromatic), loss = log(1+exp(-delta)); positive gain is better. AUC (ties count one half) is secondary. This tests within-species discrimination, not absolute white-frequency calibration or selection coefficients. Pairs are not independent replicates; uncertainty resamples genus clusters, not pairs.

Use 2,000 genus-cluster bootstrap draws of the fixed out-of-fold species metrics, with nominal 95% and approximate family2 intervals for the two primary loss-gain comparisons. These intervals describe empirical stability conditional on the fitted cross-validation procedure; they do not refit the entire training process or remove all cross-fold dependence.

## C. Prediction at held-out spatial cells within known species

Assign global 1-degree latitude-longitude cells deterministically to five folds with a distinct documented SHA256 seed. All photos from a cell enter the same fold. Run two prespecified variants: no additional distance buffer and a 50 km exclusion of training photos around every held-out conspecific photo. Use great-circle distances, including across the antimeridian. A degree cell alone does not imply a minimum distance.

Training-fold eligibility is at least 20 retained photographs, based on training counts and feature identifiability only. Do not require either response state to be present in training: constant-response fits remain eligible if numerically estimable. Test folds must contain at least one photo. Fit four models using training data only: intercept; three climate predictors; linear latitude plus circular-unwrapped longitude; geography plus climate. Predictor centering and scaling and the longitude origin use training coordinates only. Apply the existing longitude identifiability criteria and numerical rank/convergence checks. Preserve all attempted fits, failures, exclusions and extrapolation diagnostics. No adaptive predictor dropping, response-informed exclusion or alternate distance threshold.

Pre-result clarification on 2026-10-19: before any C production fits or held-out scores were inspected, fix a common smoothed intercept recalibration for prediction. For each non-intercept Firth model, hold its slopes fixed and solve for a scalar intercept shift such that its mean fitted probability on the retained training rows equals (n_white + 0.5)/(n_train + 1), the intercept-only Jeffreys-smoothed estimate on those same rows. Paired models must use identical retained training rows. Save the raw intercept, shift, adjusted intercept, target and achieved training mean. This prevents model-dimension differences in Firth's penalty from inducing different training mean calibration when adding predictors; it remains finite for constant-response training. This is an explicitly smoothed recalibration, not standard FLIC and not a guarantee of calibration on held-out locations. Use the recalibrated predictions for all C comparisons; do not choose between raw and adjusted scores after seeing outcomes. B's pairwise score differences are unaffected.

Compare climate versus intercept, and geography-plus-climate versus geography, on exactly matched successful held-out predictions for each comparison. Primary metric is baseline minus augmented binary log loss, averaged over all eligible held-out photos within species and then equally over species. Brier-score gain is secondary. Save actual test photo IDs, fold assignment and prediction values in compressed form as appropriate; save species and fold-level summaries plus model/selection diagnostics.

The four loss-gain contrasts (two comparisons times two buffer variants) form this additional family. Use 2,000 shared-genus bootstrap draws for each paired comparison and report nominal 95% and approximate family4 intervals. Report changing species/photo coverage explicitly. The buffered comparisons target transport beyond nearby observations, often under environmental extrapolation; failure does not establish absence of a local association.

## D. Molecular database extension

Expand the earlier ENA title query to include plant flower/petal pigmentation, anthocyanin and carotenoid terms. Freeze the actual query, all returned small metadata rows and project dispositions. Keep the original 131-run screen and new query denominators distinct. Exact matching to the measured roster and environmental panel is required; absent hits remain unknown. Prioritize biologically replicated white-versus-chromatic contrasts among environment-panel species for primary-paper review. Do not treat run counts as biological replicates, cultivar contrasts as natural polymorphism, repeated projects as independent origins, or candidate expression as experimentally proven causal variants. No unbounded raw-read download is required by this extension.

## Verification and interpretation

Use meaningful analytic/synthetic checks for unit and covariance conversion, genus/cell exclusion, train-only scaling, antimeridian/buffer distances, constant-response handling, matched metric arithmetic, ties, and bootstrap cluster weights. Independently review the executed result and manuscript-facing claims. Keep negative and ambiguous findings. Do not extend claims to complete phylogenetic control, causal adaptation, measured pigment costs, or fitness tradeoffs. Structured evaluation follows the distinction between sampling and transfer targets discussed by Roberts et al. (2017), https://doi.org/10.1111/ecog.02881.

Firth prediction bias toward one half and intercept correction are discussed by Puhr et al. (2017), https://doi.org/10.1002/sim.7273. The C procedure above deliberately uses a smoothed training target to retain finite constant-response training folds; it is not that paper's unmodified FLIC procedure.
