import numpy as np
from pathlib import Path
import logomaker
import matplotlib.pyplot as plt
import pandas as pd

class AnchorMiner:
    """Predicts MHC-I anchor positions for a given peptide and HLA allele.
 
    Uses precomputed PPM matrices and KL divergence vectors to identify
    canonical and non-canonical anchor positions. Canonical anchors are
    defined as P2 (index 1) and PΩ (last position).
 
    Attributes:
        alphabet: List of 20 standard amino acids in fixed order.
        peptide: Input peptide sequence.
        HLA: HLA allele string in format HLA-X##:##.
        threshold: KL divergence threshold for anchor calling.
        viz: True to include vizualizations (KL divergence per residue plot; motif logo); False otherwise
 
    Examples:
        ::
            am = AnchorMiner('KLYDWWWWKKK', 'HLA-A03:01', 0.3, 'True')
            result = am.run_anchor_miner()
    """
 
    def __init__(self, peptide, HLA, threshold,viz):
        """Initializes AnchorMiner with peptide, HLA, threshold and mode.
 
        Args:
            peptide: Peptide sequence string using single-letter amino acid
                     notation. Must consist of standard amino acids only.
            HLA: HLA allele string. Must contain 'HLA' and ':'.
                 Example: 'HLA-A02:01'.
            threshold: KL divergence threshold for anchor detection.
                       Float in range [0, 1].
            viz: Boolean. If True, generates KL divergence plot and sequence
                logo with anchor positions highlighted. Red = canonical,
                blue = non-canonical.
        """
        self.alphabet = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I',
                         'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']
        
        self.cluster_significance_threshold = 0.3

        self.clusters = {
        'aliphatic' : ['I','L', 'V','A'],
        'sulfur' : ['M','C'],
        'phenylalanine' : ['F'],
        'tryptophane' : ['W'],
        'prolyne' : ['P'],
        'glycine' : ['G'],
        'tyrosine' : ['Y'],
        'uncharged' : ['T', 'S'],
        'basic' : ['H','K','R'],
        'acidic' : ['E','D'],
        'amide' : ['Q','N'],

        }
        self.peptide = peptide
        self.HLA = HLA
        self.threshold = threshold
        self.viz = viz

 
    def run_anchor_miner(self):
        """Runs the full anchor prediction pipeline.
 
        Validates input, loads precomputed data, predicts anchors,
        and returns output in the format specified by self.mode.
 
        Returns:
            If mode='v': list of dicts, one per anchor position.
            If mode='p': list of int anchor position indices (0-indexed).
            If mode=None: tuple (coords, verbose).
 
        Raises:
            ValueError: If input validation fails or data cannot be loaded.
 
        Examples:
            ::
                am = AnchorMiner('KILDGVFAV', 'HLA-A02:01', 0.3, 'v')
                result = am.run_anchor_miner()
                print(result)
        """
        self.validate_input()
        if not self.validated_flag:
            raise ValueError('Invalid data')
        self.fetch_data()
        self.predict_anchors()

        if self.viz == True:
            self.build_KL_plot()
            self.build_logo()

        return self.make_output()
 
    def validate_input(self):
        """Validates peptide, HLA, threshold, and checks PPM file availability.
 
        Sets self.validated_flag to True if all checks pass.
 
        Raises:
            ValueError: If peptide is not a string.
            ValueError: If HLA format is invalid.
            ValueError: If threshold is out of [0, 1] range.
            ValueError: If peptide contains non-standard amino acids.
            ValueError: If HLA allele is not found in PPM dataset.
            ValueError: If peptide length has no precomputed PPM for this HLA.
            ValueError: If PPM dataset directory cannot be accessed.
 
        Examples:
            ::
                am = AnchorMiner('KILDGVFAV', 'HLA-A02:01', 0.3, 'v')
                am.validate_input()
                print(am.validated_flag)  # True
        """
        self.validated_flag = False
 
        if not isinstance(self.peptide, str) or len(self.peptide) == 0:
            raise ValueError('Peptide must be a non-empty string.')
 
        if (not isinstance(self.HLA, str)
                or 'HLA' not in self.HLA
                or ':' not in self.HLA):
            raise ValueError(
                'HLA must be a string in format HLA-X##:##. '
                'Example: HLA-A02:01'
            )
 
        if (not isinstance(self.threshold, float)
                or self.threshold < 0
                or self.threshold > 1):
            raise ValueError(
                'Threshold must be a number'
            )
 
        for aa in self.peptide:
            if aa not in self.alphabet:
                raise ValueError(
                    f'Unknown amino acid: {aa}. '  # fixed f-string bug
                    'Use single-letter notation of 20 standard amino acids.'
                )
 
        try:
            pwm_dir = Path('../datasets/PPM')
            files = [f.name for f in pwm_dir.glob('*.npy')]
        except Exception as e:
            raise ValueError(
                f'Could not access PPM dataset directory: {e}. '
                'Ensure repository structure is intact.'
            )
 
        allowed_lengths = []
        for fname in files:
            parts = fname.replace('.npy', '').rsplit('_', 1)
            if len(parts) == 2 and parts[0].replace('PPM-', '', 1) == self.HLA:
                try:
                    allowed_lengths.append(int(parts[1]))
                except ValueError:
                    continue
 
        if len(allowed_lengths) == 0:
            raise ValueError(
                f'HLA allele {self.HLA} not found in PPM dataset. '
                'Check normalization or run NetMHCpan for this allele.'
            )
 
        if len(self.peptide) not in allowed_lengths:
            raise ValueError(
                f'No precomputed PWM for peptide length {len(self.peptide)} '
                f'with HLA {self.HLA}. Available lengths: {sorted(allowed_lengths)}'
            )
 
        self.validated_flag = True
 
    def fetch_data(self):
        """Loads precomputed PPM and KL divergence arrays for HLA + peptide length.
 
        Sets self.PPM (shape 20 x L) and self.KL (shape L,).
 
        Raises:
            ValueError: If .npy files cannot be loaded.
 
        Examples:
            ::
                am = AnchorMiner('KILDGVFAV', 'HLA-A02:01', 0.3, 'v')
                am.validate_input()
                am.fetch_data()
                print(am.PWM.shape)  # (20, 9)
        """
        combination = f'{self.HLA}_{len(self.peptide)}'
        try:
            self.PWM = np.load(f'../datasets/PWM/PWM-{combination}.npy')
            self.KL = np.load(f'../datasets/KL/KL-{combination}.npy')
            self.PPM = np.load(f'../datasets/PPM/PPM-{combination}.npy')
        except Exception as e:
            raise ValueError(
                f'Could not load data for {combination}: {e}. '
                'Ensure PWM and KL files exist and are not corrupted.'
            )
 
    def predict_anchors(self):
        """Identifies anchor positions using KL threshold and cumulative frequency for aminoacid clusters.
 
        For each position with KL > threshold, fetches the cluster of aminoacids, sharing physical-chemical properties with peptide[position]; 
        If cumulative frequency for cluster > cumulative cluster threshold --> label this residue as anchor.
        Classifies anchors as canonical (P2, PΩ) or non-canonical.
 
        Sets self.res_verbose and self.coords.
 
        Examples:
            ::
                am = AnchorMiner('KILDGVFAV', 'HLA-A02:01', 0.3, 'v')
                am.validate_input()
                am.fetch_data()
                am.predict_anchors()
                print(am.coords)  # e.g. [1, 8]
        """
        canonical_positions = {1, len(self.peptide) - 1}
        anchor_positions = np.where(self.KL > self.threshold)[0]
 
        results = []
        anchor_coords = []



        def find_cluster_for_aminoacid(aa):
            ''' 
            Auxillary function, returns cluster for desired aminoacid.

            Args:
                aa (str) - aminoacid in 1 letter notation
            
            Returns:
                cluster (str) - cluster, associated with queried aminoacid

            Example:
                find_cluster_for_aminoacid('L')
                >aliphatic
        
            '''

            for i in self.clusters:
                if aa in self.clusters[i]:
                    return(i)
            
            raise ValueError('Internal error: could not assign aminoacid to cluster')
        

        for pos in anchor_positions:

            cluster_id = find_cluster_for_aminoacid(self.peptide[pos])


            cluster_freq = sum([self.PPM[self.alphabet.index(i), pos] for i in self.clusters[cluster_id]])

            if cluster_freq >= self.cluster_significance_threshold:
                is_anchor = True

                if pos in canonical_positions:

                    is_canonical = True

                else:

                    is_canonical = False


    
                results.append({
                    'position': int(pos),
                    'kl': float(round(self.KL[pos], 4)),
                    'type': 'canonical' if is_canonical else 'non-canonical',
                    'cluster': f'{cluster_id}({self.clusters[cluster_id]})',
                    'cumulative cluster frequency': cluster_freq,
                    'current_aa': self.peptide[pos],
                    'result': 'anchor' if is_anchor else 'not anchor'
                })

            else:
                is_anchor = False
 
            if is_anchor:
                anchor_coords.append(int(pos))
 
        self.res_verbose = results
        self.coords = anchor_coords


    def build_KL_plot(self):
        """Builds KL divergence per-position plot with threshold line.

        Generates a line plot of KL divergence across all peptide positions
        with a horizontal dashed line at the anchor detection threshold.
        Anchor candidate positions are those above the threshold line.
        Stores result in self.KL_plot.

        Examples:
            ::
                am = AnchorMiner('KILDGVFAV', 'HLA-A02:01', 0.3, True)
                am.validate_input()
                am.fetch_data()
                am.build_KL_plot()
                am.KL_plot.show()
        """


        fig, ax = plt.subplots(figsize=(6, 3)) 

        ax.plot(range(len(self.peptide)), self.KL, 'b', label='KL divergence')
        ax.axhline(y=self.threshold, color='r', linestyle='--', label='threshold')
        ax.set_ylabel('KL divergence')
        ax.set_xlabel('Aminoacid residue index (0-based)')
        ax.set_title(f'KL divergence vs threshold for {self.HLA} and peptide of len {len(self.peptide)}', fontsize=10)
        ax.legend()
        ax.grid()

        self.KL_plot = fig   


    def build_logo(self):

        """Builds sequence logo from PPM with anchor positions highlighted.

        Constructs an information-content logo using logomaker. Canonical
        anchor positions (P2, PΩ) are highlighted in red, non-canonical
        anchors in blue. Requires predict_anchors() to have been called first
        so that self.coords is populated.
        Stores result in self.logoplot.

        Examples:
            ::
                am = AnchorMiner('KILDGVFAV', 'HLA-A02:01', 0.3, True)
                am.validate_input()
                am.fetch_data()
                am.predict_anchors()
                am.build_logo()
                am.logoplot.fig.show()
        """

        canonical_positions = {1, len(self.peptide) - 1}  #P2, Pomega (0 based notation)
        

        ppm_matrix = pd.DataFrame({self.alphabet[i]: self.PPM[i] for i in range(len(self.alphabet))})


        information_mat = logomaker.transform_matrix(ppm_matrix, 
                                    from_type='probability', 
                                    to_type='information')
        
        logo = logomaker.Logo(information_mat, figsize=(10,5))

        for i in self.coords:
            if i in canonical_positions:
                logo.highlight_position(i, color = 'red', alpha = 0.2)
            else:
                logo.highlight_position(i, color = 'blue', alpha = 0.2)

        logo.ax.set_xlabel('Aminoacid residue index (0-based)',fontsize=14)
        logo.ax.set_ylabel("Information content (bits)", labelpad=-1,fontsize=14)
        logo.ax.set_title(f'Motif logo for {self.HLA}-binding epitopes of len {len(self.peptide)} with anchors highligted')

        self.logoplot = logo


    def make_output(self):
        """Returns prediction results in the format specified by self.mode.
 
        Returns:
            output (dict):
                ['coords'] - anchor residues coordinates
                ['verbose] - detailed information per position
            if viz == True:
                ['KL'] - KL divergence per residue plot
                ['logo'] - motif logo with anchor positions highlighted
  
        Examples:
            ::
                am = AnchorMiner('KILDGVFAV', 'HLA-A02:01', 0.3, 'p')
                am.run_anchor_miner()
                # Returns: [1, 8]
        """

        output = dict()

        output['coords'] = self.coords
        output['verbose'] = self.res_verbose
        if self.viz == True:
            output['KL'] = self.KL_plot
            output['logo'] = self.logoplot

        return(output)

