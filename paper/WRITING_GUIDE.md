# Writing guide — IJGGC CO₂ depressurisation manuscript

Notes distilled from (a) the IJGGC / Elsevier author guidelines and (b) three recent
Andreasen papers, to keep the manuscript consistent with the target journal and with the
author's established scientific voice. **Tone throughout is neutral, formal and precise —
no colourful language, no superlatives, no rhetoric.**

Style sources studied:
- Andreasen & Stegelmann (2025), *Open source pressure vessel blowdown modelling under
  partial phase equilibrium*, chemRxiv 2025-00xzc — **closest template** (same topic family,
  same lead author); its section structure maps directly onto this paper.
- Andreasen et al. (2024), *A framework for optimisation and techno-economic analysis of CO₂
  pressurisation*, chemRxiv 2024-hpnbt.
- Andreasen (2024/25), *Techno-economic analysis and optimisation of grey and green methanol
  synthesis*, chemRxiv 2024-50852.
- Andreasen (2026), *Open Source Tools for Pressure Vessel Thermo-Mechanical Response…*,
  J. Loss Prev. Process Ind. 103, 106088 (the published blowdown paper).

---

## 1. Journal requirements (IJGGC / Elsevier)

Verified against the IJGGC "Guide for Authors" (ISSN 1750-5836).

- **Reference style: Harvard (author–date).** In text (directly or parenthetically): single
  author = name + year; two authors = both names + year; three or more = first author + `et al.`
  + year. Examples: *"as demonstrated (Allan, 2020a, 2020b; Allan and Jones, 2019)"*,
  *"Kramer et al. (2023) have recently shown"*. Reference list **alphabetical then
  chronological**; same author+year distinguished by `a`, `b`, `c`. **Abbreviate journal names
  per the LTWA.** List format example:
  *"Van der Geer, J., Handgraaf, T., Lupton, R.A., 2020. The art of writing a scientific article.
  J. Sci. Commun. 163, 51–59. https://doi.org/10.1016/j.sc.2020.0037."* (article-number variant
  also given). All handled by `elsarticle-harv.bst` + `\citep{}` / `\citet{}` (already set in
  `main.tex`).
- **Structure:** IMRaD — Introduction, Methods (here "Model"), Results and discussion,
  Conclusions. Frontmatter: title, authors + affiliations, abstract, keywords, highlights.
- **Abstract:** single unstructured paragraph; concise; no citations/abbreviations.
- **Highlights: REQUIRED at submission.** A short set of bullet points capturing the novel
  results and any new methods; **≤ 85 characters each, including spaces** (Elsevier standard: 3–5
  bullets).
- **Graphical abstract:** encouraged at submission (the phase decision tree / model map is a
  strong candidate — one concise pictorial summary).
- **Keywords:** ~4–6, separated by `\sep`.
- **CRediT authorship contribution statement — required** (roles: Conceptualization, Data
  curation, Formal analysis, Funding acquisition, Investigation, Methodology, etc.).
- **Declarations — required:** Declaration of competing interests (even if none); Funding
  sources; Declaration of generative-AI use (if any).
- **Data statement / availability + data linking** — state where code/data live (HydDown
  repository; the open CARDICE and SINTEF datasets, with DOIs).
- **Units:** SI throughout; bar acceptable for pressure in this field. Symbols italic, units
  roman; number + space + unit (e.g. `5.18 bar`, `216.6 K`). Provide a nomenclature list.
- **Figures/tables:** numbered, referenced in order; **self-contained, descriptive captions**;
  vector where possible; legible in grayscale. Equations numbered.

---

## 2. Voice and tone

- **Neutral, formal, impersonal.** State facts and observations plainly; let the numbers carry
  the claim. No adjectives of enthusiasm.
- **First person plural "we" — only for the authors' own contributions, choices and paper
  structure.** e.g. *"In the present paper we propose a rigorous model…"*, *"we include…"*.
- **Passive / impersonal for everything else** — descriptions, method steps, and especially
  results: *"comparison … is made"*, *"It is also noted that …"*, *"the pressure … is predicted
  very well"*, *"Data have been sourced from ref. …"*.
- **British / `-ise` spelling** consistently: modelling, depressurisation, pressurisation,
  behaviour, optimisation, utilise, visualise, characterise, vapour, sulphur, bench-mark.
  This is the convention in the two topically-relevant papers (blowdown 2025-00xzc and CO₂
  pressurisation 2024-hpnbt). The methanol paper (2024-50852) uses American spelling
  ("optimization", "utilizing") because it is an ACS-style, co-authored manuscript — do **not**
  follow that here; IJGGC is an Elsevier journal and the CO₂/blowdown work is British.
- **CO₂** as `CO$_2$` (and `\texorpdfstring` in any heading/bookmark if hyperref complains).

## 3. Claims and hedging

Claims are **measured and honest**, always. Model–data agreement is described precisely, and
shortfalls are stated openly and quantified.

- **Prefer:** "compare well", "predicted very well", "good/adequate predictive capability",
  "within the experimental range", "slight over-/under-prediction", "approximately",
  "reasonable", "it is noted that", "it appears that … may be", "captures … with adequate
  accuracy".
- **Avoid:** "excellent", "perfect", "dramatic(ally)", "surprisingly", "clearly", "obviously",
  "significantly better" (unless quantified), "novel/first" without qualification, exclamation.
- Quantify discrepancies with magnitude and location: *"a slight over-prediction, most
  pronounced at intermediate times of 200–500 s"*; *"does not decline as rapidly as the
  experiments, but from around 60 s the temperature is within the observed range"*.
- Mechanistic explanations are offered tentatively: *"It appears that slight differences in the
  heat-transfer modelling may be one of the reasons …"*.

## 4. Structure and section patterns (from the blowdown paper)

**Introduction** — three moves:
1. Define the topic and its industrial/safety importance (here: CO₂ intermediate storage and
   transport for CCS; depressurisation; dry-ice formation and low-temperature exposure).
2. Literature review: summarise each key prior model/experiment in 2–4 sentences — its
   approach and findings — **and note its limitation/gap**. e.g. *"Witkowski and Majkut (Year)
   investigated 13 … concepts … Some of their findings include … However, the authors do not
   include …"*.
3. Purpose and contribution: *"The purpose of the present paper is to …"*; position as an
   open-source, rigorously validated alternative; state objectives (often an enumerated list).

**Model / Methods** — a general overview first (with a schematic figure of the heat/mass-transfer
processes), then one numbered subsection per sub-model: thermodynamics/phase equilibrium →
mass & energy balance → discharge → heat transfer → (dry-ice / below-triple) → code
implementation. State assumptions explicitly; number all equations; cite the open-source
packages used (CoolProp, fluids, thermopack for the offline solid table).

**Results and discussion** — 
- Open with the validation intent: *"In order to investigate the quality of the presented …
  model, comparison with relevant experiments … is made."*
- List the cases/experiments with **data provenance** (*"Data have been sourced from ref. …"*)
  and a summary **Table** with a detailed caption (initial conditions, geometry, Cd, assumptions).
- One subsection per case/regime (here: low/medium-pressure CARDICE; dense-phase SINTEF).
- Report agreement precisely, note discrepancies honestly, compare to other models/tools where
  relevant, and give hedged mechanistic reasoning.

**Conclusions** — concise restatement of what was done, the validated capability, the key
quantitative findings, and the honest limitations (e.g. 0-D stratification).

## 5. Citations, figures, tables, reproducibility

- Weave citations into the sentence: *"a much more recent published model by Park et al.
  (Year) share many similarities …"* / *"For a comprehensive review … the recent study of
  Shafiq et al. (Year) is referred to."*
- Reference figures/tables as *"As seen from Figure X …"* / *"summarised in Table Y …"*.
- Captions are self-contained: define every symbol and every assumption in the caption itself.
- Emphasise **openness and reproducibility** — the tool is open source and free; datasets are
  open; cite the exact packages and data DOIs. This is a recurring, deliberate theme.

## 6. Ready-made phrasings (author's voice)

- Opening: *"Rapid depressurisation of pressure vessels … also referred to as blowdown, is an
  essential part of the plant process safety measures."*
- Contribution: *"The purpose of the present paper is to present a newly developed … model
  which has been made available to the public as open source and free of charge."*
- Validation intent: *"In order to investigate the quality of the presented … model,
  comparison with relevant experiments from the literature is made."*
- Agreement: *"As seen from Figure X, the calculations compare well with the experimental
  results."* / *"the time-dependent pressure and the vessel wall temperature are predicted
  very well."*
- Honest shortfall: *"The calculated … does not decline as rapidly as the experiments, but
  from around … s the temperature is within the experimentally observed range."*
- Hedged mechanism: *"It appears that slight differences in the heat-transfer modelling may be
  one of the reasons for these subtle differences between models."*

## 7. This paper — working section map (mirrors the blowdown paper)

1. Introduction — CCS CO₂ transport/intermediate storage; dry-ice and cold-temperature hazards;
   prior CO₂-release models and experiments (CARDICE, SINTEF, Hammer, Munkejord) with gaps;
   purpose (open-source, thermopack-free, validated against two campaigns).
2. Experimental basis — CARDICE (Ineris 2 m³ sphere) and Høydalsvik/Munkejord (SINTEF dense-phase).
3. Model — thermodynamics & the solid-CO₂ lever; discharge (HEM + non-equilibrium); in-vessel
   dry ice (two-zone plateau/descent); heat & mass transfer.
4. Results and validation — low/medium-pressure (CARDICE); dense-phase (SINTEF).
5. Discussion — 0-D limitations (stratification), model choices, applicability.
6. Conclusions.

Much of the physics/experimental prose can be adapted from `docs/techref/` (already written in
a matching neutral register) into manuscript form — but re-cast it in the paper's voice
(measured claims, Harvard citations, British spelling) rather than copying verbatim.
