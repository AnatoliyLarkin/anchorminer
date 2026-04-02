import numpy as np
import sys
import subprocess
from src.AnchorMiner import AnchorMiner


def predict_anchors(pept, hla, threshold, viz):



    try:
        AM = AnchorMiner(pept,hla,threshold,viz)
    except Exception as e:
        raise ValueError(f'Could not initialize anchor miner: {e}')
    
    return(AM.run_anchor_miner())
    



if __name__ == '__main__':
    predict_anchors()