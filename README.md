# AnchorMiner — MHC-I Anchor Position Predictor

Anchor residues are amino acids within a peptide that form direct contacts with the MHC-I binding groove. They serve a structural docking function rather than defining T cell recognition, making their identification useful for immunogenicity modeling, epitope masking, and feature engineering in machine learning pipelines.

AnchorMiner identifies canonical and non-canonical anchor positions for a given peptide and MHC-I allele using allele-specific Position Weight Matrices (PWMs) built from NetMHCpan-4.2 predictions on a peptide background sampled from the human and mouse SwissProt proteomes.

For a detailed description of the method and usage examples, see `notebooks/Demo.ipynb`.

---

## Requirements

See `requirements.txt`. Install with:

```bash
pip install -r requirements.txt
pip install anchorminer
```

---

## Quick Start
```python
from anchorminer import predict_anchors

# Without visualization
result = predict_anchors('KILDGVFAV', 'HLA-A*02:01', viz=False)

# With visualization — generates KL plot and sequence logo
result = predict_anchors('KILDGVFAV', 'HLA-A*02:01', viz=True)

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

## Supported Alleles and Peptide Lengths

 - Homo sapiens: 146 MHC-I alleles covering 92% of worldwide population. Peptide lengths 9–12. 
 - Mus musculus: 10 alleles are supported.
 - Check `data/Anchor_Miner_supportedalleles.csv` to get a dataframe of allowed HLAs and lengths.

---

## Input formats

- Homo sapiens HLA allele format: `HLA-X*##:##` e.g. `HLA-A*02:01`.
- Mus Musculus MHC format: `H2-K*b`.
- mhcgnomes (https://github.com/pirl-unc/mhcgnomes) string formatting or manual data curation is strongly recommended


## Performance
 ~ 400 epitope-HLA pairs / sec