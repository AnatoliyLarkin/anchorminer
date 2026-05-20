import re

def normalize_allele(raw):
    DROP_PATTERNS = ['undetermined', 'hla-class', 'hla-bx', 'allele und', 'class i']
    DROP_SEROLOGY = {'B12', 'B15', 'B5', 'B40', 'A28', 'A10', 'A9', 'B70', 'B17', 'B14', 'B24', 'B62', 'B61', 'B63', 'B60'}

    MHC1_LOCI = {'A', 'B', 'C', 'E'}
    MHC2_LOCI = {'DRB1', 'DRB3', 'DRB4', 'DRB5', 'DQA1', 'DQB1', 'DPA1', 'DPB1'}
    VALID_LOCI = MHC1_LOCI | MHC2_LOCI

    MOUSE_ALLELES = {
        'D': {'b', 'd', 'q'},
        'K': {'b', 'd', 'k', 'q'},
        'L': {'d', 'q'},
        'Qa1': {''},
        'Qa2': {''},
    }

    if not isinstance(raw, str):
        return None

    s = raw.strip()
    if not s or s.lower() == 'nan':
        return None

    for pat in DROP_PATTERNS:
        if pat in s.lower():
            return None

    if '__' in s:
        parts = [normalize_allele(p) for p in s.split('__')]
        return [p for p in parts if p and not isinstance(p, list)]

    if ',' in s:
        parts = [normalize_allele(p.strip()) for p in s.split(',')]
        flat = []
        for p in parts:
            if isinstance(p, list):
                flat.extend(p)
            elif p:
                flat.append(p)
        return flat if flat else None

    # Mouse H-2 alleles
    mouse_match = re.match(
        r'^H[-_]2[-_]([A-Za-z]+\d?)([a-z]?)$',
        s,
        flags=re.IGNORECASE
    )
    if mouse_match:
        locus_raw = mouse_match.group(1)
        suffix    = mouse_match.group(2).lower()

        locus_norm = locus_raw[0].upper() + locus_raw[1:].lower()
        locus_norm = re.sub(r'Qa(\d)', lambda m: f'Qa{m.group(1)}', locus_norm)

        allowed = MOUSE_ALLELES.get(locus_norm)
        if allowed is None:
            return None

        if suffix not in allowed:
            return None

        if suffix:
            return f'H-2-{locus_norm}{suffix}'
        else:
            return f'H-2-{locus_norm}'

    # HLA alleles
    s = s.split()[0]
    s = s.rstrip('?')
    s = re.sub(r'\(.*?\)', '', s).strip()
    s = re.sub(r'^HLA[-_]?', '', s, flags=re.IGNORECASE)
    s = re.sub(r'^Cw', 'C', s, flags=re.IGNORECASE)

    match = re.match(r'^([A-Za-z]+\d?)\*?([\d:_\.]+)?$', s)
    if not match:
        return None

    locus = match.group(1).upper()
    numstr = match.group(2) or ''

    if locus not in VALID_LOCI:
        return None

    if not numstr and locus in DROP_SEROLOGY:
        return None

    numstr = re.sub(r'[:_\.]', '', numstr)
    numstr = numstr[:4]

    if len(numstr) == 0:
        return None
    elif len(numstr) == 1:
        field1 = f"0{numstr}"
        field2 = "01"
    elif len(numstr) == 2:
        field1 = numstr
        field2 = "01"
    elif len(numstr) == 3:
        field1 = numstr[:2]
        field2 = numstr[2:] + "0"
    else:
        field1 = numstr[:2]
        field2 = numstr[2:]

    return f"HLA-{locus}{field1}:{field2}"