import re

def normalize_allele(raw):

    DROP_PATTERNS = ['h-2', 'hla-?', 'undetermined', 'hla-class', 'hla-bx', 'allele und', 'class i']
    DROP_SEROLOGY = {'B12', 'B15', 'B5', 'B40', 'A28', 'A10', 'A9', 'B70', 'B17', 'B14', 'B24', 'B62', 'B61', 'B63', 'B60'}
    
    MHC1_LOCI = {'A', 'B', 'C', 'E'}
    MHC2_LOCI = {'DRB1', 'DRB3', 'DRB4', 'DRB5', 'DQA1', 'DQB1', 'DPA1', 'DPB1'}
    VALID_LOCI = MHC1_LOCI | MHC2_LOCI

    """
    Normalize a single HLA allele string to HLA-X*FF:FF format.
    Returns None if invalid/ambiguous, or a list if input contains multiple alleles.
    """
    if not isinstance(raw, str):
        return None
 
    s = raw.strip()
    if not s or s.lower() == 'nan':
        return None
 
    # Drop known garbage
    for pat in DROP_PATTERNS:
        if pat in s.lower():
            return None
 
    # Split haplotypes (DQA1_05__DQB1_02 style)
    if '__' in s:
        parts = [normalize_allele(p) for p in s.split('__')]
        return [p for p in parts if p and not isinstance(p, list)]
 
    # Split comma-separated multi-allele entries
    if ',' in s:
        parts = [normalize_allele(p.strip()) for p in s.split(',')]
        flat = []
        for p in parts:
            if isinstance(p, list):
                flat.extend(p)
            elif p:
                flat.append(p)
        return flat if flat else None
 
    # Strip trailing mutation notes e.g. "HLA-B*27:05 L81A mutant"
    s = s.split()[0]
 
    # Strip trailing ?
    s = s.rstrip('?')
 
    # Remove parenthetical e.g. (HLA-B*4427)
    s = re.sub(r'\(.*?\)', '', s).strip()
 
    # Normalize prefix: remove HLA- or HLA_ 
    s = re.sub(r'^HLA[-_]?', '', s, flags=re.IGNORECASE)
 
    # Cw -> C (old serological notation)
    s = re.sub(r'^Cw', 'C', s, flags=re.IGNORECASE)
 
    # Handle mouse (shouldn't reach here but safety net)
    if s.startswith('H-2') or s.startswith('h-2'):
        return None
 
    # Extract locus and number string
    # Handles: A*02:01, A0201, A02, A2, DRB1*01:01, DRB1_01_01
    match = re.match(r'^([A-Za-z]+\d?)\*?([\d:_\.]+)?$', s)
    if not match:
        return None
 
    locus = match.group(1).upper()
    numstr = match.group(2) or ''
 
    if locus not in VALID_LOCI:
        return None
 
    # Broad serology drop (ambiguous splits like B15, B40)
    if not numstr and locus in DROP_SEROLOGY:
        return None
 
    # Strip non-digit separators from number string
    numstr = re.sub(r'[:_\.]', '', numstr)
 
    # Truncate to 4 digits max (handles typos like 02304 -> 0230)
    numstr = numstr[:4]
 
    if len(numstr) == 0:
        # Pure locus only e.g. "HLA-A" — too ambiguous
        return None
    elif len(numstr) == 1:
        # e.g. A2 -> 02:01
        field1 = f"0{numstr}"
        field2 = "01"
    elif len(numstr) == 2:
        # e.g. A02 -> 02:01
        field1 = numstr
        field2 = "01"
    elif len(numstr) == 3:
        # e.g. A021 -> ambiguous, treat as 02:10 (rare, flag it)
        field1 = numstr[:2]
        field2 = numstr[2:] + "0"
    else:
        # 4 digits: 0201 -> 02:01
        field1 = numstr[:2]
        field2 = numstr[2:]
 
    return f"HLA-{locus}{field1}:{field2}"
