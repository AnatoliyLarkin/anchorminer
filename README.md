# AnchorMiner — MHC-I Anchor Position Predictor

Anchor residues are amino acids within a peptide that form direct contacts with the MHC-I binding groove. They serve a structural docking function rather than defining T cell recognition, making their identification useful for immunogenicity modeling, epitope masking, and feature engineering in machine learning pipelines.

AnchorMiner identifies canonical and non-canonical anchor positions for a given peptide and MHC-I allele using allele-specific Position Weight Matrices (PWMs) built from NetMHCpan-4.2 predictions on a 600,000-peptide background sampled from the human SwissProt proteome.

For a detailed description of the method and usage examples, see `notebooks/Demo.ipynb`.

---

## Requirements

See `requirements.txt`. Install with:

```bash
pip install -r requirements.txt
```

---

## Quick Start
```python
from src.predict_anchors import predict_anchors

# Without visualization
result = predict_anchors('KILDGVFAV', 'HLA-A02:01', viz=False)

# With visualization — generates KL plot and sequence logo
result = predict_anchors('KILDGVFAV', 'HLA-A02:01', viz=True)

# Access results
print(result['coords'])   # anchor position indices (0-indexed)
print(result['verbose'])  # per-position details
```

**Output:** dict with keys:
- `coords` — list of anchor position indices (0-indexed)
- `verbose` — list of dicts per anchor position with KL divergence, canonical/non-canonical classification, biochemical cluster, cumulative cluster frequency, and current amino acid
- `KL` — KL divergence plot (only if viz=True)
- `logo` — sequence logo with anchor positions highlighted in red (canonical) or blue (non-canonical) (only if viz=True)

---

## Repository Structure
```
AnchorMiner/
├── src/
│   ├── AnchorMiner.py          # Core class
│   ├── predict_anchors.py      # Wrapper function
│   └── standartize_hlas.py     # HLA normalization utility
├── datasets/
│   ├── PWM/                    # Precomputed PWM matrices (.npy)
│   ├── PPM/                    # Precomputed PPM matrices (.npy)
│   ├── KL/                     # Precomputed KL divergence vectors (.npy)
│   └── AnchorMiner_supported_alleles.csv
├── notebooks/
│   ├── Demo.ipynb              # Tutorial and usage examples
│   ├── Testing.ipynb           # Sanity checks and structural validation
│   └── AnchorMiner_DevelopmentNotebook.ipynb
├── requirements.txt
└── README.md
```


---

## Supported Alleles and Peptide Lengths

146 MHC-I alleles covering 92% of worldwide population. Peptide lengths 9–12. Visit `datasets/Anchor_Miner_supportedalleles.csv` to get a dataframe of allowed HLAs and lengths.

---

## Notes

- Homo sapiens HLA allele format: `HLA-X##:##` e.g. `HLA-A02:01`. Use `src/standartize_hlas.normalize_allele()` to normalize non-standard inputs. Only MHCI alleles are currently supported
- Mus Musculus MHC format: `H:2:Kd`.
