"""
Core algorithms for primer design and analysis
"""

import math
from typing import Dict, List, Tuple, Optional
from Bio.Seq import Seq
from Bio.SeqUtils import molecular_weight
import re

class PrimerAnalyzer:
    """Core algorithms for primer design and quality assessment"""
    
    # Nearest neighbor thermodynamic parameters (ΔH and ΔS)
    # Values in kcal/mol for ΔH and cal/(mol·K) for ΔS
    NN_PARAMS = {
        'AA': {'dH': -7.6, 'dS': -21.3}, 'AT': {'dH': -7.2, 'dS': -20.4},
        'AC': {'dH': -8.5, 'dS': -22.7}, 'AG': {'dH': -7.8, 'dS': -21.0},
        'TA': {'dH': -7.2, 'dS': -21.3}, 'TT': {'dH': -7.6, 'dS': -21.3},
        'TC': {'dH': -8.2, 'dS': -22.2}, 'TG': {'dH': -8.5, 'dS': -22.7},
        'CA': {'dH': -8.5, 'dS': -22.7}, 'CT': {'dH': -7.8, 'dS': -21.0},
        'CC': {'dH': -8.0, 'dS': -19.9}, 'CG': {'dH': -10.6, 'dS': -27.2},
        'GA': {'dH': -8.2, 'dS': -22.2}, 'GT': {'dH': -8.4, 'dS': -22.4},
        'GC': {'dH': -9.8, 'dS': -24.4}, 'GG': {'dH': -8.0, 'dS': -19.9}
    }
    
    # Terminal penalties
    TERMINAL_AT = {'dH': 2.3, 'dS': 4.1}
    TERMINAL_GC = {'dH': 0.1, 'dS': -2.8}
    
    def __init__(self):
        self.salt_correction = 0.368  # Default salt correction factor
        
    def calculate_melting_temperature(self, sequence: str, 
                                    primer_conc: float = 0.25, 
                                    salt_conc: float = 0.05) -> float:
        """
        Calculate melting temperature using nearest neighbor method
        
        Args:
            sequence: DNA sequence string
            primer_conc: Primer concentration in μM
            salt_conc: Salt concentration in M
            
        Returns:
            Melting temperature in Celsius
        """
        sequence = sequence.upper().replace('U', 'T')
        
        if len(sequence) < 2:
            return 0.0
            
        # Calculate thermodynamic parameters
        delta_h = 0.0  # kcal/mol
        delta_s = 0.0  # cal/(mol·K)
        
        # Add nearest neighbor contributions
        for i in range(len(sequence) - 1):
            dinucleotide = sequence[i:i+2]
            if dinucleotide in self.NN_PARAMS:
                delta_h += self.NN_PARAMS[dinucleotide]['dH']
                delta_s += self.NN_PARAMS[dinucleotide]['dS']
        
        # Add terminal penalties
        if sequence[0] in 'AT':
            delta_h += self.TERMINAL_AT['dH']
            delta_s += self.TERMINAL_AT['dS']
        else:
            delta_h += self.TERMINAL_GC['dH']
            delta_s += self.TERMINAL_GC['dS']
            
        if sequence[-1] in 'AT':
            delta_h += self.TERMINAL_AT['dH']
            delta_s += self.TERMINAL_AT['dS']
        else:
            delta_h += self.TERMINAL_GC['dH']
            delta_s += self.TERMINAL_GC['dS']
        
        # Convert to appropriate units
        delta_h *= 1000  # Convert to cal/mol
        
        # Calculate Tm with salt correction
        R = 1.987  # Gas constant in cal/(mol·K)
        
        # Basic Tm calculation
        tm_kelvin = delta_h / (delta_s + R * math.log(primer_conc * 1e-6 / 4))
        
        # Salt correction (von Ahsen et al., 2001)
        tm_celsius = tm_kelvin - 273.15
        tm_corrected = tm_celsius + 16.6 * math.log10(salt_conc) + 0.62
        
        return round(tm_corrected, 1)
    
    def calculate_gc_content(self, sequence: str) -> float:
        """Calculate GC content percentage"""
        if not sequence:
            return 0.0
        sequence = sequence.upper()
        gc_count = sequence.count('G') + sequence.count('C')
        total_bases = len(sequence)
        return round((gc_count / total_bases) * 100, 1)
    
    def calculate_primer_quality_score(self, sequence: str) -> Dict[str, float]:
        """
        Calculate comprehensive primer quality score
        
        Returns:
            Dictionary with individual scores and overall quality
        """
        sequence = sequence.upper().replace('U', 'T')
        length = len(sequence)
        
        scores = {
            'length_score': 0.0,
            'gc_score': 0.0,
            'tm_score': 0.0,
            'complexity_score': 0.0,
            'terminal_score': 0.0,
            'overall_score': 0.0
        }
        
        # Length score (optimal 18-25 bp)
        if 18 <= length <= 25:
            scores['length_score'] = 100.0
        elif 16 <= length < 18 or 25 < length <= 30:
            scores['length_score'] = 80.0
        elif 14 <= length < 16 or 30 < length <= 35:
            scores['length_score'] = 60.0
        else:
            scores['length_score'] = 20.0
        
        # GC content score (optimal 40-60%)
        gc_content = self.calculate_gc_content(sequence)
        if 40 <= gc_content <= 60:
            scores['gc_score'] = 100.0
        elif 30 <= gc_content < 40 or 60 < gc_content <= 70:
            scores['gc_score'] = 80.0
        elif 20 <= gc_content < 30 or 70 < gc_content <= 80:
            scores['gc_score'] = 60.0
        else:
            scores['gc_score'] = 20.0
        
        # Tm score (optimal 55-65°C)
        tm = self.calculate_melting_temperature(sequence)
        if 55 <= tm <= 65:
            scores['tm_score'] = 100.0
        elif 50 <= tm < 55 or 65 < tm <= 70:
            scores['tm_score'] = 80.0
        elif 45 <= tm < 50 or 70 < tm <= 75:
            scores['tm_score'] = 60.0
        else:
            scores['tm_score'] = 20.0
        
        # Complexity score (avoid runs and repeats)
        complexity_penalty = 0
        
        # Check for runs of same nucleotide (>3 consecutive)
        for base in 'ATCG':
            if base * 4 in sequence:
                complexity_penalty += 20
        
        # Check for dinucleotide repeats
        for i in range(len(sequence) - 5):
            if sequence[i:i+2] == sequence[i+2:i+4] == sequence[i+4:i+6]:
                complexity_penalty += 15
        
        scores['complexity_score'] = max(0, 100 - complexity_penalty)
        
        # Terminal score (prefer G/C at 3' end, avoid runs at 3' end)
        terminal_penalty = 0
        
        # Check 3' end stability
        if sequence[-1] in 'AT':
            terminal_penalty += 10
        
        # Check for 3' end runs
        if sequence[-3:] in ['AAA', 'TTT', 'CCC', 'GGG']:
            terminal_penalty += 20
        
        scores['terminal_score'] = max(0, 100 - terminal_penalty)
        
        # Calculate overall score (weighted average)
        weights = {
            'length_score': 0.2,
            'gc_score': 0.25,
            'tm_score': 0.25,
            'complexity_score': 0.2,
            'terminal_score': 0.1
        }
        
        scores['overall_score'] = sum(
            scores[key] * weights[key] for key in weights
        )
        
        return {k: round(v, 1) for k, v in scores.items()}
    
    def detect_primer_dimers(self, primer1: str, primer2: str, 
                           min_overlap: int = 4) -> List[Dict]:
        """
        Detect potential primer dimers between two primers
        
        Args:
            primer1: First primer sequence
            primer2: Second primer sequence
            min_overlap: Minimum overlap length to consider
            
        Returns:
            List of potential dimer structures
        """
        primer1 = primer1.upper().replace('U', 'T')
        primer2 = primer2.upper().replace('U', 'T')
        
        dimers = []
        
        # Check for 3' complementarity (most critical)
        for i in range(min_overlap, min(len(primer1), len(primer2)) + 1):
            # Check primer1 3' end vs primer2 3' end
            seq1_3prime = primer1[-i:]
            seq2_3prime = primer2[-i:]
            
            # Reverse complement check
            complement = str(Seq(seq2_3prime).reverse_complement())
            
            matches = sum(1 for a, b in zip(seq1_3prime, complement) if a == b)
            
            if matches >= min_overlap:
                stability_score = self._calculate_dimer_stability(seq1_3prime, complement)
                
                dimers.append({
                    'type': '3prime_dimer',
                    'primer1_region': seq1_3prime,
                    'primer2_region': seq2_3prime,
                    'matches': matches,
                    'total_length': i,
                    'match_percentage': (matches / i) * 100,
                    'stability_score': stability_score,
                    'risk_level': self._assess_dimer_risk(matches, i, stability_score)
                })
        
        # Check for hairpin structures in individual primers
        for primer_name, primer_seq in [('primer1', primer1), ('primer2', primer2)]:
            hairpins = self._detect_hairpins(primer_seq)
            for hairpin in hairpins:
                hairpin['primer'] = primer_name
                dimers.append(hairpin)
        
        return sorted(dimers, key=lambda x: x.get('stability_score', 0), reverse=True)
    
    def _calculate_dimer_stability(self, seq1: str, seq2: str) -> float:
        """Calculate thermodynamic stability of dimer formation"""
        if len(seq1) != len(seq2):
            return 0.0
        
        # Simple stability calculation based on GC content and length
        gc_content = sum(1 for i, (a, b) in enumerate(zip(seq1, seq2)) 
                        if a == b and a in 'GC') / len(seq1)
        
        match_ratio = sum(1 for a, b in zip(seq1, seq2) if a == b) / len(seq1)
        
        return round(gc_content * match_ratio * len(seq1) * 2, 1)
    
    def _assess_dimer_risk(self, matches: int, total_length: int, 
                          stability_score: float) -> str:
        """Assess the risk level of dimer formation"""
        match_percentage = (matches / total_length) * 100
        
        if stability_score > 15 or (matches >= 4 and match_percentage > 75):
            return 'HIGH'
        elif stability_score > 8 or (matches >= 3 and match_percentage > 60):
            return 'MEDIUM'
        else:
            return 'LOW'
    
    def _detect_hairpins(self, sequence: str, min_stem_length: int = 4) -> List[Dict]:
        """Detect hairpin structures in a single primer"""
        hairpins = []
        seq_len = len(sequence)
        
        for i in range(seq_len - min_stem_length * 2):
            for j in range(i + min_stem_length * 2, seq_len):
                # Check for complementarity
                region1 = sequence[i:i + min_stem_length]
                region2 = sequence[j:j + min_stem_length]
                
                complement = str(Seq(region2).reverse_complement())
                matches = sum(1 for a, b in zip(region1, complement) if a == b)
                
                if matches >= min_stem_length:
                    loop_size = j - i - min_stem_length
                    
                    if 3 <= loop_size <= 10:  # Reasonable loop size
                        stability = self._calculate_dimer_stability(region1, complement)
                        
                        hairpins.append({
                            'type': 'hairpin',
                            'stem_sequence': region1,
                            'loop_size': loop_size,
                            'matches': matches,
                            'stability_score': stability,
                            'risk_level': self._assess_dimer_risk(matches, min_stem_length, stability)
                        })
        
        return hairpins
    
    def analyze_primer_sequence(self, sequence: str) -> Dict:
        """
        Comprehensive analysis of a single primer sequence
        
        Returns:
            Complete analysis including all metrics
        """
        analysis = {
            'sequence': sequence.upper(),
            'length': len(sequence),
            'gc_content': self.calculate_gc_content(sequence),
            'melting_temperature': self.calculate_melting_temperature(sequence),
            'molecular_weight': round(molecular_weight(Seq(sequence), 'DNA'), 1),
            'quality_scores': self.calculate_primer_quality_score(sequence),
            'hairpins': self._detect_hairpins(sequence)
        }
        
        return analysis