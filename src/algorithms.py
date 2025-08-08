"""
Core algorithms for primer design and analysis
"""

import math
from typing import Dict, List, Tuple, Optional
from Bio.Seq import Seq
from Bio.SeqUtils import molecular_weight
import re
import itertools

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
        return round((gc_count / len(sequence)) * 100, 1)
    
    def calculate_primer_quality_score(self, sequence: str) -> Dict[str, float]:
        """
        Calculate comprehensive primer quality scores
        
        Args:
            sequence: DNA sequence string
            
        Returns:
            Dictionary with quality scores for different aspects
        """
        if not sequence:
            return {
                'length_score': 0.0,
                'gc_score': 0.0,
                'tm_score': 0.0,
                'complexity_score': 0.0,
                'terminal_score': 0.0,
                'overall_score': 0.0
            }
        
        sequence = sequence.upper()
        length = len(sequence)
        
        # Length score (optimal: 18-25 bp)
        if 18 <= length <= 25:
            length_score = 1.0
        elif 15 <= length <= 30:
            length_score = 0.8
        elif 12 <= length <= 35:
            length_score = 0.6
        else:
            length_score = 0.3
        
        # GC content score (optimal: 40-60%)
        gc_content = self.calculate_gc_content(sequence)
        if 40 <= gc_content <= 60:
            gc_score = 1.0
        elif 30 <= gc_content <= 70:
            gc_score = 0.8
        elif 20 <= gc_content <= 80:
            gc_score = 0.6
        else:
            gc_score = 0.3
        
        # Melting temperature score (optimal: 55-65°C)
        tm = self.calculate_melting_temperature(sequence)
        if 55 <= tm <= 65:
            tm_score = 1.0
        elif 50 <= tm <= 70:
            tm_score = 0.8
        elif 45 <= tm <= 75:
            tm_score = 0.6
        else:
            tm_score = 0.3
        
        # Sequence complexity score
        complexity_score = self._calculate_complexity_score(sequence)
        
        # Terminal stability score
        terminal_score = self._calculate_terminal_score(sequence)
        
        # Overall score (weighted average)
        overall_score = (
            length_score * 0.2 +
            gc_score * 0.2 +
            tm_score * 0.25 +
            complexity_score * 0.2 +
            terminal_score * 0.15
        )
        
        return {
            'length_score': round(length_score, 3),
            'gc_score': round(gc_score, 3),
            'tm_score': round(tm_score, 3),
            'complexity_score': round(complexity_score, 3),
            'terminal_score': round(terminal_score, 3),
            'overall_score': round(overall_score, 3)
        }
    
    def _calculate_complexity_score(self, sequence: str) -> float:
        """Calculate sequence complexity score"""
        # Check for repetitive patterns
        repetitive_score = 1.0
        
        # Check for homopolymers (runs of same nucleotide)
        for base in 'ATCG':
            max_run = max(len(list(g)) for k, g in itertools.groupby(sequence) if k == base)
            if max_run > 4:
                repetitive_score -= 0.2 * (max_run - 4)
        
        # Check for dinucleotide repeats
        dinucleotide_repeats = re.findall(r'(..)\1{2,}', sequence)
        if dinucleotide_repeats:
            repetitive_score -= 0.1 * len(dinucleotide_repeats)
        
        return max(0.0, repetitive_score)
    
    def _calculate_terminal_score(self, sequence: str) -> float:
        """Calculate terminal stability score"""
        if len(sequence) < 2:
            return 0.0
        
        # Prefer G/C at 3' end for better extension
        terminal_score = 1.0
        
        if sequence[-1] in 'GC':
            terminal_score += 0.2
        elif sequence[-1] in 'AT':
            terminal_score -= 0.1
        
        # Avoid G at 5' end (can cause secondary structure)
        if sequence[0] == 'G':
            terminal_score -= 0.1
        
        return max(0.0, min(1.0, terminal_score))
    
    def detect_primer_dimers(self, primer1: str, primer2: str, 
                           min_overlap: int = 4) -> List[Dict]:
        """
        Detect potential primer-dimer formations
        
        Args:
            primer1: First primer sequence
            primer2: Second primer sequence
            min_overlap: Minimum overlap length to consider
            
        Returns:
            List of dimer structures found
        """
        dimers = []
        
        # Check all possible overlaps
        for overlap_len in range(min_overlap, min(len(primer1), len(primer2)) + 1):
            # Check 3' end of primer1 with 3' end of primer2
            if primer1[-overlap_len:] == self._reverse_complement(primer2[-overlap_len:]):
                dimers.append({
                    'type': '3-3',
                    'overlap_length': overlap_len,
                    'overlap_sequence': primer1[-overlap_len:],
                    'stability_score': self._calculate_dimer_stability(
                        primer1[-overlap_len:], 
                        primer2[-overlap_len:]
                    ),
                    'risk_level': self._assess_dimer_risk(overlap_len, overlap_len, 0.8)
                })
            
            # Check 3' end of primer1 with 5' end of primer2
            if primer1[-overlap_len:] == self._reverse_complement(primer2[:overlap_len]):
                dimers.append({
                    'type': '3-5',
                    'overlap_length': overlap_len,
                    'overlap_sequence': primer1[-overlap_len:],
                    'stability_score': self._calculate_dimer_stability(
                        primer1[-overlap_len:], 
                        primer2[:overlap_len]
                    ),
                    'risk_level': self._assess_dimer_risk(overlap_len, overlap_len, 0.6)
                })
        
        return dimers
    
    def _calculate_dimer_stability(self, seq1: str, seq2: str) -> float:
        """Calculate stability score for dimer formation"""
        # Simple scoring based on GC content and length
        gc_content = self.calculate_gc_content(seq1)
        length = len(seq1)
        
        # Higher GC content and longer length = more stable
        stability = (gc_content / 100) * (length / 20)
        return min(1.0, stability)
    
    def _assess_dimer_risk(self, matches: int, total_length: int, 
                          stability_score: float) -> str:
        """Assess risk level of dimer formation"""
        risk_score = (matches / total_length) * stability_score
        
        if risk_score > 0.7:
            return 'high'
        elif risk_score > 0.4:
            return 'medium'
        else:
            return 'low'
    
    def _detect_hairpins(self, sequence: str, min_stem_length: int = 4) -> List[Dict]:
        """Detect potential hairpin structures"""
        hairpins = []
        
        # Check for palindromic sequences that could form hairpins
        for i in range(len(sequence) - min_stem_length):
            for j in range(i + min_stem_length, len(sequence)):
                stem_length = j - i
                if stem_length > 8:  # Limit stem length
                    continue
                
                stem_seq = sequence[i:j]
                if stem_seq == self._reverse_complement(stem_seq):
                    hairpins.append({
                        'start': i,
                        'end': j,
                        'stem_length': stem_length,
                        'stem_sequence': stem_seq
                    })
        
        return hairpins
    
    def _reverse_complement(self, sequence: str) -> str:
        """Get reverse complement of sequence"""
        complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
        return ''.join(complement.get(base, base) for base in reversed(sequence))
    
    def analyze_primer_sequence(self, sequence: str) -> Dict:
        """
        Comprehensive primer sequence analysis
        
        Args:
            sequence: DNA sequence string
            
        Returns:
            Dictionary with analysis results
        """
        if not sequence:
            return {}
        
        sequence = sequence.upper()
        
        return {
            'sequence': sequence,
            'length': len(sequence),
            'gc_content': self.calculate_gc_content(sequence),
            'melting_temperature': self.calculate_melting_temperature(sequence),
            'quality_scores': self.calculate_primer_quality_score(sequence),
            'hairpins': self._detect_hairpins(sequence),
            'molecular_weight': molecular_weight(Seq(sequence), 'DNA')
        }