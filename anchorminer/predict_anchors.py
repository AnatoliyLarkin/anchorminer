from .AnchorMiner import AnchorMiner


def predict_anchors(pept, hla, viz=False):
    try:
        AM = AnchorMiner(pept, hla, viz)
    except Exception as e:
        raise ValueError(f'Could not initialize anchor miner: {e}')
    
    return AM.run_anchor_miner()

