# LGB Population Estimation from Survey + ACS (Block-Group Resolution)

This repository contains an R workflow that:

1. **Estimates LGB identification rates** by **state × sex × age bin** using Household Pulse Survey (HPS) microdata (pooled across weekly files).
2. **Projects LGB population counts to small areas** using **ACS 5-year age–sex counts (table B01001)** at the **block-group** level.
3. **Reweights those projections spatially** using the observed **concentration of same-sex couples** (ACS relationship tables) as a proxy for within-state geographic clustering.

The intended use is *descriptive*: producing **expected** LGB population counts/shares consistent with (i) survey-based LGB rates by demographics and (ii) ACS age–sex population structure, with an additional optional calibration step to better match the within-state spatial pattern implied by same-sex couples.

---

## What the estimates *mean* (interpretation)

At any target geography (e.g., a buffer around an address, a block group, a county), the model produces:

- **Expected LGB count**: \(\sum_{a,s} \text{Pop}_{a,s}\times \Pr(\text{LGB}\mid a,s,\text{state})\)
- **Expected LGB share**: Expected LGB count divided by total population in the same aggregation.

Key point: **these are not observed LGB counts**. They are **model-based expectations** under the assumptions listed below. They should be interpreted similarly to *small-area estimates* that combine survey rates with census population structure.

### Uncalibrated vs calibrated outputs

The current script produces both:

- **Uncalibrated expected LGB**: driven purely by age–sex composition and the state-specific LGB rates.
- **Calibrated expected LGB**: the uncalibrated county totals are *reallocated within state* to align with the county distribution of **same-sex couples**, using a multiplicative **calibration factor**.

Intuitively, calibration shifts expected LGB population **toward counties where same-sex couples are more prevalent** (relative to opposite-sex couples), without changing the overall within-state total implied by the demographic rates (unless you choose a calibration scheme that changes totals).

---

## Data sources

### 1) Household Pulse Survey (HPS) microdata
Used to estimate demographic LGB rates. The script expects multiple HPS CSV files and a person weight variable.

- Key fields used in the current implementation (names may differ depending on HPS vintage):  
  - Sexual orientation (`SEXUAL_ORIENTATION`)
  - Sex assigned at birth / reported sex (`EGENID_BIRTH`)
  - Birth year (`TBIRTH_YEAR`)
  - State FIPS (`EST_ST`)
  - Survey week (`WEEK`)
  - Person weight (`PWEIGHT`)

### 2) ACS 5-year estimates via `tidycensus`
Used for small-area population denominators and the couples-based calibration:

- **Age–sex counts (B01001)** at block-group level
- **Same-sex vs opposite-sex couple counts** (current script uses tables in the B09019 family/relationship series)

> Note: ACS couples tables reflect cohabiting couples/households as defined by ACS. They are **not** a count of LGB individuals and are imperfect proxies for spatial concentration.

---

## Methodology (high level)

### Step A — Estimate state × sex × age-bin LGB rates from HPS

1. Compute respondent age: `age = ref_year - TBIRTH_YEAR`.
2. Map individual ages into ACS-compatible bins (18–19, 20, 21, 22–24, ..., 85+).
3. Map sexual orientation responses into categories:
   - `LG` (lesbian/gay)  
   - `Bisexual`  
   - `Straight` (includes nonresponse/other in current script)

4. Within each `state × sex × age_bin × week`, compute weighted totals using `PWEIGHT`.
5. Average across survey weeks to reduce week-to-week sampling noise.
6. Convert to shares:
   - `lg_shr = LG / (LG + Straight + Bisexual)`
   - `lgb_shr = (LG + Bisexual) / (LG + Straight + Bisexual)`

**Output:** a table of \(\Pr(\text{LG})\) and \(\Pr(\text{LGB})\) by state, sex, age bin.

### Step B — Project expected LGB counts to block groups (or any aggregation)

For each block group:

1. Pull ACS B01001 age–sex counts.
2. Convert to a long format with `sex` and `acs_bin`.
3. Merge on the estimated state-specific LGB rates.
4. Compute expected LGB counts:
   - `expected_lgb = pop × lgb_shr` (and similarly for `LG`)

These expected counts can then be aggregated to:
- counties (for calibration),
- a buffer around an address (current script),
- a user-defined region/polygon.

### Step C — Spatial calibration using same-sex couples (within-state)

The script computes a county-level **couples intensity**:

- `ss = same-sex spouse + same-sex partner`
- `os = opposite-sex spouse + opposite-sex partner`
- `intensity = ss / max(os, 1)`

Then it:

1. Aggregates uncalibrated expected LGB counts to county and normalizes to a county share.
2. Creates a **target** county distribution by tilting that uncalibrated share using `intensity^gamma` (default `gamma = 1`).
3. Forms a multiplicative adjustment:
   - `calib_factor = target_share / uncal_share`
4. Applies `calib_factor` back to small areas within each county.

> The calibration changes the within-state spatial distribution of expected LGB counts. It **does not create new information** about individuals; it is a model choice intended to better reflect geographic clustering that is not captured by age–sex composition alone.

---

## Outputs and how to read them

The example function `estimate_lgb_for_address()` returns:

- **Population totals** in the chosen radius (e.g., `pop`, `pop_22_29`, `pop_m_22_29`)
- **Expected LG/LGB counts** overall and for ages 22–29 split by sex (uncalibrated)
- **A calibrated LGB total** (`est_lgb_cal`) and implied calibration factor (`local_cal_factor`)
- **Calibrated subgroup counts and shares** for ages 22–29 (e.g., `shr_m_22_29_lgb`)

Suggested interpretation:

- Use **`est_lgb_cal`** (and corresponding calibrated subgroup shares) if you believe the couples-based proxy helps approximate within-state spatial concentration.
- Use **`est_lgb`** (uncalibrated) if you want estimates driven strictly by demographic composition and state-level rates.

---

## Key assumptions (what must be true for the estimates to be valid)

1. **Transportability of rates:**  
   Within each state, LGB identification probabilities depend only on **sex and age bin** (as modeled), and those rates apply to every small area in that state.

2. **Alignment of definitions:**  
   The HPS sexual-orientation measure and the target concept of “LGB” are consistent over time and comparable to what users mean by LGB population.

3. **Survey weighting adequacy:**  
   The person weights (`PWEIGHT`) correctly reweight HPS respondents to represent the state population by sex and age.

4. **ACS measurement accuracy:**  
   ACS B01001 age–sex counts are accurate at block-group resolution (not always true for small populations).

5. **Calibration proxy validity (if used):**  
   Same-sex couples (ACS) are an informative proxy for where LGB individuals live, **after** controlling for age–sex composition.

6. **Temporal coherence:**  
   The “ref_year” used to compute age and the ACS vintage year refer to the same population period (or differences are negligible).

---

## Caveats and limitations (important)

- **Identity and underreporting:**  
  Sexual orientation is sensitive, and survey respondents may underreport or misreport. Rates may vary by context, cohort, and local norms.

- **Not the same as LGBT:**  
  The current workflow is **LGB** (sexual orientation) and does not estimate transgender/nonbinary populations.

- **Small-area uncertainty is not quantified (yet):**  
  The script currently reports point estimates, not confidence intervals. In reality, uncertainty comes from:
  - sampling error in HPS rates,
  - ACS sampling error at block-group level,
  - model specification uncertainty (especially calibration).

- **Same-sex couples ≠ LGB individuals:**  
  Couples data reflect household structure and reporting; they miss single adults, non-cohabiting relationships, and may be affected by differential household formation.

- **Calibration can amplify noise:**  
  Counties with small ACS couple denominators can generate extreme intensities. Consider:
  - smoothing intensity measures,
  - trimming/winsorizing calibration factors,
  - setting `gamma < 1`,
  - minimum-count thresholds.

- **Geocoding + privacy:**  
  The address-based function geocodes an address and returns estimates for a buffer. Do **not** publish point-level estimates that could be used to infer sensitive traits about individuals.

---

## Reproducibility notes

To run locally you will need:

- R packages: `tidycensus`, `tigris`, `sf`, `dplyr`, `tidyr`, `purrr`, `stringr`, `units`, `tidygeocoder`, plus `data.table` (for `fread`) and `readxl` if using the Excel address import.
- A Census API key (see `tidycensus::census_api_key()`).

**Recommended repository structure**
```
data/
  hps_raw/            # HPS CSVs (microdata)
  acs_cache/          # optional caching
outputs/
  rates/              # estimated state×sex×age rates
  bg_estimates/        # block-group expected counts
  maps/               # geo outputs for web app
R/
  lgb_estimation.R
  helpers.R
```

The current script contains **hard-coded local file paths**; replacing these with relative paths (or config/env vars) is recommended before publishing.

---

## Planned interactive mapping tool (GitHub Pages / Shiny / Leaflet)

An interactive tool can help users explore results safely at appropriate geographic resolution.

### Minimal viable features
- Choose **state** and **year** (ACS vintage).
- Display a **choropleth map** of:
  - expected LGB share,
  - expected LGB count,
  - optionally calibrated vs uncalibrated toggle.
- Hover tooltips with:
  - GEOID, total pop, expected LGB, expected share,
  - and warning banners for small denominators.
- Search bar to zoom to a place.
- Export/download:
  - aggregated results,
  - a static map image.

### Safe defaults
- Show results at **block-group** or **tract** only when population is sufficiently large (or aggregate to tract/county by default).
- Disable address lookup (or return only coarse results, e.g., within a ZIP or tract).
- Include a prominent **disclaimer** about sensitive inference.

### Technical options
- **R Shiny + leaflet**: quickest if you stay in R.
- **Quarto + Observable/leaflet**: good for GitHub Pages, static hosting.
- **Python (Panel/Folium) or JS (Mapbox/Deck.gl)**: good for a pure web front-end.

---

## Suggested citation / acknowledgement

If you use or build on this workflow, cite:
- U.S. Census Bureau Household Pulse Survey (HPS)
- American Community Survey (ACS) 5-year estimates
- `tidycensus` and `tigris` R packages

---

## License and ethics

This work produces modeled estimates for a sensitive characteristic. Users should:
- avoid releasing fine-grained maps that could increase risk for individuals or communities,
- document assumptions, and
- consult IRB / ethics guidance when necessary.
