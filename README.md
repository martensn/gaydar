# Gaydar: Estimating Local LGBTQ Populations

Have you ever wondered how many LGBTQ people live within a ten mile radius of any address in the United States **or** what share of the people in your community are non-binary?
Gaydar answers this question (and any related inquiries about the spatial distribution of different members of the LGBTQ+ community). 
I focus on adult (18+) gay and lesbian, bisexual, queer, and transgender populations (sorry asexuals), split into three mutually exclusive and collectively exhaustive gender identities: men, women, and non-binary people. 

Note transgender people might also identify as gay, lesbian, bisexual, or queer, meaning the sub-population shares will not sum to the total because I don't double count LGBQ transgender people. 
Additionally, limitations within Household Pulse Survey microdata mean that I cannot measure non-binary transgender individuals. 

---

## Summary

[Methodology](https://martensn.shinyapps.io/gaydar/methods.html) offers a precise explanation of the process through which I derived the estimates, though I summarize the steps below:

1. Estimate LGBTQ identification rates by state × sex assigned at birth × age bin using pooled Household Pulse Survey (HPS) microdata, with survey-weighted estimation and soft calibration to recent Gallup national benchmarks.

2. Estimate gender identity composition conditional on sexual orientation, sex assigned at birth, and age, using HPS microdata to recover age-specific, distributions of different gender identities. I assume cisgender, non-binary, and transgender are mutually exclusive and collectively exhaustive.
Step 2 constructs a crosswalk from ACS age-sex counts to age-gender identity-sexuality estimates. 
   These composition estimates are stabilized using **age-local Dirichlet smoothing**, allowing gender identity shares to vary smoothly across cohorts while remaining responsive to well-identified differences by age. 
   
3. Estimate LGBTQ sub-population shares (i.e. lesbians and gays, bisexuals, queers, transgender people) by pooling nationally. State-level cells (think transgender people in any state or bisexuals in Wyoming) were too sparse to trust even with empirical Bayes shrinkage across states, so every state now applies the same national LGBTQ mix to its own state-level identification rate from step 1.

4. Project expected LGBTQ population counts to Census tract by combining estimated identification rates with ACS 5-year age–sex population counts (table B01001). Since age and sex vary locally, sexuality and gender identity will too. 

5. Reweight projected LGBTQ populations within states by blending two PUMA-level spatial concentration proxies:  same-sex couples and singles within the same age-sex cell.
Same-sex couples tend to be older, wealthier, better educated, suggesting their spatial distribution might systematically differ from younger populations. 
To remedy the bias among smae-sex couples relative to the broader LGBTQ population, I constructed a weighted average of the two PUMA-level concentration measures, using the national HPS marriage rate for a cell as the weight. Reweighting ensures the estimates reflect more than the mere age and sex composition of a Census tract
   If I naively rescaled tract-level LGBTQ populations based only on the PUMA-level signal, LGBTQ population shares would show implausibly large breaks around PUMA boundaries. Instead, I rescale the tract estimates using inverse distance weighting (IDW). IDW blends rescaling factors from every PUMA within a state, giving the most weight to nearby PUMAs. Thus, IDW generates distinct rescaling factors for tracts within the same PUMA based on their proximity to neighboring PUMAs.

6. In a small number of tracts (typically tracts near higher education institutions or other clusters of youth) the LGBTQ populations must be capped because the PUMA reweighing scales an already-high LGBTQ share above 100 percent. I rescale the population in these tracts, minutely reducing the state-level population estimates but improving the plausibility of local estimates.

---

## Data Sources

### Household Pulse Survey (HPS)

HPS was the first the Census data product to ask respondents about sexual orientation and gender identity, 
basing the design of questions on the results of a [2022 NASEM report](https://www.nationalacademies.org/publications/26424).
Designed to elicit responses to rapidly changing situations (i.e. food security early in the coronavirus pandemic),
the HPS surveyed different households each week. 
While topics varied between weeks, HPS consistently collected the following demographic data:

* Sexual orientation
* Gender identity
* Sex assigned at birth
* Year of birth
* Marital status
* State of residence
* Survey weights

The microdata allow the relationship between sexual orientation and gender identity to vary flexibly by age and geography. For example, younger gay and lesbian individuals are more likely to identify as non-binary than older cohorts.

---

### Gallup Polling

The [February 2026 Gallup report](https://news.gallup.com/poll/702206/lgbtq-identification-holds.aspx) provides the most recent nationally representative estimate of LGBTQ identification in the United States.

The report includes:

* Overall national LGBTQ identification
* LGBTQ identification by generation
* LGBTQ identification by generation × sex (for select cohorts)

Given uncertainty in survey measurement, these estimates are used as **soft calibration targets**, not hard constraints.
In practice, most state-level estimates fall within one to two percentage points of Gallup’s top-line estimate of 9 percent.

---

### American Community Survey (ACS)

Population projections are anchored to **ACS 5-year estimates**, which provide stable small-area population counts at fine geographic resolution.

Specifically, the pipeline uses:

* [Table B01001 (Sex by Age)](https://censusreporter.org/tables/B01001/) to obtain population counts by **sex assigned at birth** and granular age bins, measured at the **tract** level.
* [Table B09019 (Household Type by Relationship)](https://censusreporter.org/tables/B09019/) informs the same-sex-couples half of the spatial reweighting above, measured at the **PUMA** level.
* [Table B12002 (Sex by Marital Status by Age)](https://censusreporter.org/tables/B12002/) provides the singles-concentration signal blended into the same reweighting, also measured at the **PUMA** level.

Because the Household Pulse Survey (HPS) and ACS use different age groupings, HPS respondents are mapped into ACS-compatible age bins prior to estimation. This ensures that estimated LGBTQ identification rates and gender-composition shares can be directly applied to ACS population counts without interpolation or extrapolation.

---
