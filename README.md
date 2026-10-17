# Drug-Likeness Screener

[![Streamlit App](https://static.streamlit.io/badges/streamlit_badge_black_white.svg)](https://drug-likeness-screener-nd9ybhbd8rcnmg3nxdxxme.streamlit.app/)
![Python](https://img.shields.io/badge/python-3.x-blue)
![RDKit](https://img.shields.io/badge/RDKit-cheminformatics-green)

An interactive cheminformatics tool that screens chemical libraries against **Lipinski's Rule of Five**. It computes key molecular descriptors with RDKit, flags rule violations for every compound, and visualizes candidate properties so weak candidates can be filtered out early.

**Live demo:** https://drug-likeness-screener-nd9ybhbd8rcnmg3nxdxxme.streamlit.app/

<!-- TODO: add a screenshot. Save it as docs/screenshot.png, then uncomment the next line.
![App screenshot](docs/screenshot.png)
-->

---

## Why this exists

Before a compound goes into docking or wet-lab testing, it helps to know whether it has obvious oral drug-likeness problems. Lipinski's Rule of Five is the standard first-pass filter. This app automates that check for a whole library instead of one molecule at a time, which makes early lead-compound triage faster and reproducible.

## Features

- Evaluates a library of compounds against Lipinski's Rule of Five
- Computes key molecular descriptors using RDKit
- Identifies and counts rule violations per compound
- Visualizes candidate properties to compare compounds at a glance
- Runs in the browser through Streamlit, with nothing to install for the live version

<!-- TODO: confirm against app.py and add or remove bullets:
     - accepted input format (SMILES column? CSV upload? pasted text?)
     - any extra descriptors beyond MW / logP / HBD / HBA (TPSA, rotatable bonds, ...)
     - whether results can be downloaded as CSV
     - which chart types are used -->

## The Rule of Five

| Property | Lipinski threshold |
|---|---|
| Molecular weight | ≤ 500 Da |
| Calculated logP | ≤ 5 |
| Hydrogen-bond donors (HBD) | ≤ 5 |
| Hydrogen-bond acceptors (HBA) | ≤ 10 |

A compound is generally considered to have a poor oral drug-likeness profile when it breaks more than one of these criteria.

## How it works

1. Each molecule is read into RDKit (typically from a SMILES string).
2. Descriptors are calculated for every molecule.
3. Each descriptor is compared against its Lipinski threshold.
4. Violations are counted per compound and compounds are flagged accordingly.
5. Results are displayed in tables and charts in the Streamlit interface.

<!-- TODO: edit these steps so they match app.py exactly. -->

## Run it locally

```bash
git clone https://github.com/kannegantimeenakshi5-max/Drug-likeness-screener.git
cd Drug-likeness-screener
pip install -r requirements.txt
streamlit run app.py
```

The repo also includes a `.devcontainer` folder, so it can be opened directly in a GitHub Codespace or VS Code Dev Container.

## Limitations

- **Rule-based only.** This is a physicochemical filter. It does not predict absorption, metabolism, toxicity, or activity.
- **Built for oral small molecules.** Natural products, macrocycles, peptides, and many antibiotics legitimately break the Rule of Five and still work as drugs. A violation is a flag to interpret, not an automatic rejection.
- **Calculated values.** logP and the other descriptors are RDKit estimates, not measured values.

## Possible extensions

- Additional filters such as Veber rules (TPSA, rotatable bonds)
- PAINS and structural-alert screening using RDKit's filter catalogs
- CSV export of screening results
- Batch upload of larger libraries

## Tech stack

Python · RDKit · Streamlit

## Author

**Kanneganti Meenakshi Naga Prasanna**
B.Tech Bioinformatics, Vignan University
[GitHub](https://github.com/kannegantimeenakshi5-max)
