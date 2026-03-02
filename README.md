# Gaydar: Estimating Local LGBTQ Populations

This repository contains an R workflow that:

1. **Estimates LGB identification rates** by **state × sex × age bin** using Household Pulse Survey (HPS) microdata (pooled across weekly files).
2. **Projects LGB population counts to small areas** using **ACS 5-year age–sex counts (table B01001)** at the **block-group** level.
3. **Crosswalks from sex to gender identity** using HPS estimates of gender identity within age–sexuality bins, enabling analysis of non-binary and transgender populations.
4. **Reweights those projections spatially** using the observed **concentration of same-sex couples** (ACS relationship tables) as a proxy for within-state geographic clustering.

The intended use is *descriptive*: producing **expected** LGBTQ population counts and shares rather than causal estimates.

---

## Data Sources

### 1. Household Pulse Survey (HPS)

Repeated cross-sectional microdata collected from Weeks 40–60, containing:

* Sexual orientation
* Gender identity
* Sex assigned at birth
* Year of birth
* State of residence
* Survey weights

The microdata allow the relationship between sexual orientation and gender identity to vary flexibly by age and geography.
For example, younger gay and lesbian individuals are more likely to identify as non-binary than older cohorts.

---

### 2. Gallup LGBTQ Identification

The [February 2026 Gallup report](https://news.gallup.com/poll/702206/lgbtq-identification-holds.aspx) provides the most recent nationally representative estimate of LGBTQ identification in the United States.

The report includes:

* Overall national LGBTQ identification
* LGBTQ identification by generation
* LGBTQ identification by generation × sex (for select cohorts)

Given uncertainty in survey measurement, these estimates are used as **soft calibration targets**, not hard constraints.
In practice, most state-level estimates fall within one to two percentage points of Gallup’s top-line estimate of 9.1%.

---

### 3. ACS Age Bins

HPS respondents are mapped into ACS-compatible age bins to ensure downstream compatibility with ACS population counts.

---

## Step 1: Measurement in HPS Microdata

From raw HPS responses:

* Sex is defined using sex assigned at birth.
* Age is inferred from birth year and mapped to generational categories.
* The LGBT indicator is defined as:

  * 1 if sexual orientation ≠ straight **or** gender identity = transgender.
* LGBT subcategories include:

  * Lesbian or Gay
  * Bisexual
  * Queer
  * Trans

Survey weights (`PWEIGHT`) are retained throughout.
At this stage, estimates are purely survey-based and reflect sampling variability.

---

## Step 2: Soft Calibration to Gallup Marginals

### Motivation

HPS and Gallup differ in:

* Survey mode
* Sampling variation
* Field timing

Rather than forcing exact agreement, this step nudges HPS-based probabilities toward Gallup benchmarks while preserving within-cell variation.

---

### Method

Each individual begins with a baseline probability:

[
p_{0i} = \text{clipped}(Y_i)
]

A logit adjustment model is estimated:

[
\text{logit}(p_i) = \text{logit}(p_{0i}) + X_i \theta
]

where:

* ( X_i ) includes generation effects and generation × sex interactions
* ( \theta ) minimizes:

  1. A deviation penalty (keeps estimates close to HPS baselines)
  2. A margin penalty (penalizes deviations from Gallup targets)

Optimization is performed using BFGS.

The tuning parameter ( \lambda ) controls calibration strength.
I set ( \lambda = 0.5 ), allowing moderate geographic concentration without assuming all LGBTQ individuals are clustered in a few destinations.

---

### Interpretation

* Targets are **soft**, not exact constraints
* Adjustments occur on the logit scale
* Output is an adjusted probability ( p_i^{\text{adj}} ) for each respondent

---

## Step 3: State–Cohort LGBT Rate Estimation

Using calibrated probabilities:

[
\widehat{P}(Y=1 \mid S, \text{State}, A)
========================================

\frac{\sum_i w_i p_i^{\text{adj}}}{\sum_i w_i}
]

Effective sample size is computed as:

[
n_{\text{eff}} = \frac{(\sum w_i)^2}{\sum w_i^2}
]

Weekly waves are averaged to reduce short-term noise.

---

## Step 4: Age-Local Dirichlet Smoothing of Gender Composition

This step estimates the distribution of gender identity
(cisgender, non-binary, transgender) within each:

* Sexual orientation category
* Sex assigned at birth
* Exact age (top-coded at 62)

Raw single-year age estimates can be noisy, so composition shares are regularized using **age-local Dirichlet shrinkage**.

---

### Notation

Fix a sexual orientation category ( \ell ) and sex ( s ). Let:

* ( a ) denote exact age
* ( g \in {\text{cis}, \text{non_binary}, \text{trans}} ) index gender categories
* ( n_{g,\ell,s,a} ) be the survey-weighted count
* ( n_{\ell,s,a} = \sum_g n_{g,\ell,s,a} )
* ( \hat{p}*{g,\ell,s,a} = n*{g,\ell,s,a} / n_{\ell,s,a} )

Effective sample size:

[
n_{\text{eff},\ell,s,a}
=======================

\frac{(\sum w_i)^2}{\sum w_i^2}
]

---

### Age-Local Prior Construction

For a focal age ( a ), all other ages with positive effective sample size are treated as neighbors.

Each neighboring age ( a' \neq a ) receives weight:

[
w_{a'}^{(a)}
============

n_{\text{eff},\ell,s,a'}
\cdot
\exp!\left(
-\frac{|a - a'|}{\tau}
\right)
]

where ( \tau > 0 ) controls the smoothness of age pooling.

The age-local prior mean is:

[
m_{g}^{(a)}
===========

\frac{
\sum_{a' \neq a}
w_{a'}^{(a)} \hat{p}*{g,\ell,s,a'}
}{
\sum*{a' \neq a}
w_{a'}^{(a)}
}
]

---

### Dirichlet Moment Matching

A Dirichlet prior is centered at
( \mathbf{m}^{(a)} = (m_g^{(a)}) ).

The concentration parameter is estimated via method of moments:

[
\alpha_0
========

\text{mean}_g
\left(
\frac{m_g (1 - m_g)}{v_g} - 1
\right)
]

where ( v_g ) is the weighted variance across neighbor ages.

The prior vector is:

[
\alpha_g = \alpha_0 , m_g
]

---

### Posterior Shrinkage

The smoothed estimate for focal age ( a ) is:

[
\tilde{p}_{g,\ell,s,a}
======================

\frac{
n_{\text{eff},\ell,s,a} \hat{p}*{g,\ell,s,a}
+
\alpha_g
}{
n*{\text{eff},\ell,s,a}
+
\alpha_0
}
]

Properties:

* Shrinkage adapts to effective sample size
* Borrowing is local in age, not across states
* Large cells remain close to raw estimates
* Output is a valid probability simplex

After smoothing, estimates are aggregated to ACS age bins using weighted averages.

---

## Step 6: Cross-State Shrinkage of LGBT Subcategory Shares

To stabilize state-level estimates of LGBT subcategory composition, a separate empirical Bayes procedure is applied across states.

For each sex ( s ), ACS age bin ( a ), and LGBT category ( \ell ):

[
\hat{\pi}_{\ell,s,a,j}
======================

P(\ell \mid Y=1, s, a, \text{state } j)
]

with effective sample size ( n_{\text{eff},j} ).

---

### Beta Prior Estimation

For each ( (s, a, \ell) ), a beta prior is fit across states using method of moments:

[
\pi_{\ell,s,a,j}
\sim
\text{Beta}(\alpha, \beta)
]

---

### Posterior Shrinkage

State-level estimates are shrunk toward the cross-state mean:

[
\tilde{\pi}_{\ell,s,a,j}
========================

\frac{
n_{\text{eff},j} \hat{\pi}*{\ell,s,a,j}
+
\alpha
}{
n*{\text{eff},j}
+
\alpha + \beta
}
]

Final shares are renormalized to sum to one within each state × sex × age bin.

---

If you want next:

* a **shortened README-safe version**
* or a **formal methods appendix** version
  just say the word.
