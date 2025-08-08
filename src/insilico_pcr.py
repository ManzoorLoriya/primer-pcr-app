import re
from Bio.Seq import Seq
from Bio.SeqUtils import MeltingTemp as mt
# from Bio.SeqUtils.GC import gc_fraction
import logging
from collections import namedtuple
import numpy as np

def calculate_gc_content(sequence):
    """Calculate GC content percentage of a sequence."""
    gc_count = sequence.upper().count('G') + sequence.upper().count('C')
    return (gc_count / len(sequence)) * 100 if sequence else 0

# Define data structures
AmpliconResult = namedtuple('AmpliconResult', [
    'forward_start', 'forward_end', 'reverse_start', 'reverse_end',
    'amplicon_size', 'amplicon_sequence', 'template_id', 'strand'
])

class InSilicoPCR:
    """
    Performs in-silico PCR simulation to predict amplification products
    """
    
    def __init__(self, max_product_size=5000, min_product_size=50):
        self.max_product_size = max_product_size
        self.min_product_size = min_product_size
        self.logger = logging.getLogger(__name__)
    
    def simulate_pcr(self, templates, forward_primer, reverse_primer, 
                    max_mismatches=2, min_tm_diff=5):
        """
        Simulate PCR amplification across multiple template sequences
        
        Args:
            templates (dict): Template sequences {id: sequence}
            forward_primer (str): Forward primer sequence
            reverse_primer (str): Reverse primer sequence
            max_mismatches (int): Maximum allowed mismatches
            min_tm_diff (float): Minimum Tm difference for binding
            
        Returns:
            list: Predicted amplification products
        """
        self.logger.info(f"Simulating PCR with {len(templates)} templates")
        
        all_amplicons = []
        
        for template_id, template_seq in templates.items():
            # Find primer binding sites
            forward_sites = self._find_binding_sites(
                template_seq, forward_primer, max_mismatches, 'forward'
            )
            reverse_sites = self._find_binding_sites(
                template_seq, reverse_primer, max_mismatches, 'reverse'
            )
            
            # Generate amplicons for all valid primer pairs
            amplicons = self._generate_amplicons(
                template_seq, forward_sites, reverse_sites, template_id
            )
            
            all_amplicons.extend(amplicons)
        
        # Filter by size constraints
        valid_amplicons = [
            amp for amp in all_amplicons 
            if self.min_product_size <= amp.amplicon_size <= self.max_product_size
        ]
        
        self.logger.info(f"Found {len(valid_amplicons)} valid amplicons")
        return valid_amplicons
    
    def _find_binding_sites(self, template, primer, max_mismatches, direction):
        """
        Find all potential primer binding sites in template
        
        Args:
            template (str): Template sequence
            primer (str): Primer sequence
            max_mismatches (int): Maximum allowed mismatches
            direction (str): 'forward' or 'reverse'
            
        Returns:
            list: Binding sites with positions and mismatch counts
        """
        binding_sites = []
        primer_len = len(primer)
        
        # For reverse primer, work with reverse complement
        if direction == 'reverse':
            search_primer = str(Seq(primer).reverse_complement())
        else:
            search_primer = primer
        
        # Sliding window search
        for i in range(len(template) - primer_len + 1):
            template_segment = template[i:i + primer_len]
            mismatches = self._count_mismatches(search_primer, template_segment)
            
            if mismatches <= max_mismatches:
                binding_sites.append({
                    'position': i,
                    'mismatches': mismatches,
                    'sequence': template_segment,
                    'direction': direction
                })
        
        return binding_sites
    
    def _count_mismatches(self, seq1, seq2):
        """Count mismatches between two sequences"""
        if len(seq1) != len(seq2):
            return len(seq1)  # Treat as complete mismatch
        
        return sum(1 for a, b in zip(seq1, seq2) if a != b)
    
    def _generate_amplicons(self, template, forward_sites, reverse_sites, template_id):
        """
        Generate amplicons from primer binding sites
        
        Args:
            template (str): Template sequence
            forward_sites (list): Forward primer binding sites
            reverse_sites (list): Reverse primer binding sites
            template_id (str): Template identifier
            
        Returns:
            list: Amplicon results
        """
        amplicons = []
        
        for fwd_site in forward_sites:
            for rev_site in reverse_sites:
                # Ensure forward site comes before reverse site
                if fwd_site['position'] >= rev_site['position']:
                    continue
                
                # Calculate amplicon size
                amplicon_size = rev_site['position'] - fwd_site['position']
                
                # Check size constraints
                if not (self.min_product_size <= amplicon_size <= self.max_product_size):
                    continue
                
                # Extract amplicon sequence
                amplicon_seq = template[fwd_site['position']:rev_site['position']]
                
                # Create amplicon result
                amplicon = AmpliconResult(
                    forward_start=fwd_site['position'],
                    forward_end=fwd_site['position'] + len(fwd_site['sequence']),
                    reverse_start=rev_site['position'],
                    reverse_end=rev_site['position'] + len(rev_site['sequence']),
                    amplicon_size=amplicon_size,
                    amplicon_sequence=amplicon_seq,
                    template_id=template_id,
                    strand='+'
                )
                
                amplicons.append(amplicon)
        
        return amplicons
    
    def analyze_amplicon_specificity(self, amplicons, target_template_ids=None):
        """
        Analyze specificity of predicted amplicons
        
        Args:
            amplicons (list): List of amplicon results
            target_template_ids (list): List of target template IDs
            
        Returns:
            dict: Specificity analysis results
        """
        if not amplicons:
            return {
                'total_amplicons': 0,
                'target_amplicons': 0,
                'off_target_amplicons': 0,
                'specificity_score': 0.0
            }
        
        # Separate target and off-target amplicons
        target_amplicons = []
        off_target_amplicons = []
        
        for amplicon in amplicons:
            if target_template_ids and amplicon.template_id in target_template_ids:
                target_amplicons.append(amplicon)
            else:
                off_target_amplicons.append(amplicon)
        
        # Calculate specificity metrics
        total_amplicons = len(amplicons)
        target_count = len(target_amplicons)
        off_target_count = len(off_target_amplicons)
        
        # Specificity score (higher is better)
        if total_amplicons == 0:
            specificity_score = 0.0
        else:
            specificity_score = target_count / total_amplicons
        
        return {
            'total_amplicons': total_amplicons,
            'target_amplicons': target_count,
            'off_target_amplicons': off_target_count,
            'specificity_score': specificity_score,
            'size_distribution': self._calculate_size_distribution(amplicons),
            'problematic_products': self._identify_problematic_products(
                off_target_amplicons, target_amplicons
            )
        }
    
    def _calculate_size_distribution(self, amplicons):
        """Calculate size distribution of amplicons"""
        if not amplicons:
            return {}
        
        sizes = [amp.amplicon_size for amp in amplicons]
        
        return {
            'min_size': min(sizes),
            'max_size': max(sizes),
            'mean_size': sum(sizes) / len(sizes),
            'size_ranges': {
                'small (<200bp)': len([s for s in sizes if s < 200]),
                'medium (200-500bp)': len([s for s in sizes if 200 <= s <= 500]),
                'large (>500bp)': len([s for s in sizes if s > 500])
            }
        }
    
    def _identify_problematic_products(self, non_target_amplicons, target_amplicons):
        """Identify problematic amplification products"""
        problematic = []
        
        # Check for off-target products similar in size to targets
        if target_amplicons:
            target_sizes = [amp.amplicon_size for amp in target_amplicons]
            target_size_range = (min(target_sizes) * 0.8, max(target_sizes) * 1.2)
            
            for amp in non_target_amplicons:
                if target_size_range[0] <= amp.amplicon_size <= target_size_range[1]:
                    problematic.append({
                        'type': 'size_conflict',
                        'amplicon': amp,
                        'issue': f'Off-target product ({amp.amplicon_size}bp) similar to target size'
                    })
        
        # Check for very small or very large products
        for amp in non_target_amplicons:
            if amp.amplicon_size < 100:
                problematic.append({
                    'type': 'small_product',
                    'amplicon': amp,
                    'issue': f'Very small off-target product ({amp.amplicon_size}bp)'
                })
            elif amp.amplicon_size > 2000:
                problematic.append({
                    'type': 'large_product',
                    'amplicon': amp,
                    'issue': f'Very large off-target product ({amp.amplicon_size}bp)'
                })
        
        return problematic
    
    def predict_pcr_efficiency(self, amplicons, primer_tm_forward, primer_tm_reverse):
        """
        Predict PCR efficiency based on primer properties and amplicon characteristics
        
        Args:
            amplicons (list): List of amplicon results
            primer_tm_forward (float): Forward primer melting temperature
            primer_tm_reverse (float): Reverse primer melting temperature
            
        Returns:
            dict: Efficiency predictions
        """
        if not amplicons:
            return {
                'overall_efficiency': 'unknown',
                'efficiency_factors': [],
                'recommendations': []
            }
        
        efficiency_factors = []
        recommendations = []
        
        # Analyze primer Tm difference
        tm_diff = abs(primer_tm_forward - primer_tm_reverse)
        if tm_diff <= 2:
            efficiency_factors.append({
                'factor': 'primer_tm_difference',
                'status': 'good',
                'value': f'{tm_diff:.1f}°C',
                'description': 'Primer Tm difference is optimal'
            })
        elif tm_diff <= 5:
            efficiency_factors.append({
                'factor': 'primer_tm_difference',
                'status': 'acceptable',
                'value': f'{tm_diff:.1f}°C',
                'description': 'Primer Tm difference is acceptable'
            })
        else:
            efficiency_factors.append({
                'factor': 'primer_tm_difference',
                'status': 'poor',
                'value': f'{tm_diff:.1f}°C',
                'description': 'Large Tm difference may reduce efficiency'
            })
            recommendations.append('Consider redesigning primers to have similar Tm')
        
        # Analyze amplicon sizes
        sizes = [amp.amplicon_size for amp in amplicons]
        avg_size = sum(sizes) / len(sizes)
        
        if 100 <= avg_size <= 500:
            efficiency_factors.append({
                'factor': 'amplicon_size',
                'status': 'good',
                'value': f'{avg_size:.0f}bp',
                'description': 'Optimal amplicon size range'
            })
        elif avg_size < 100:
            efficiency_factors.append({
                'factor': 'amplicon_size',
                'status': 'poor',
                'value': f'{avg_size:.0f}bp',
                'description': 'Very small amplicons may be difficult to detect'
            })
            recommendations.append('Consider increasing amplicon size')
        else:
            efficiency_factors.append({
                'factor': 'amplicon_size',
                'status': 'acceptable',
                'value': f'{avg_size:.0f}bp',
                'description': 'Large amplicons may have reduced efficiency'
            })
        
        # Determine overall efficiency
        good_factors = len([f for f in efficiency_factors if f['status'] == 'good'])
        poor_factors = len([f for f in efficiency_factors if f['status'] == 'poor'])
        
        if poor_factors == 0 and good_factors >= len(efficiency_factors) * 0.7:
            overall_efficiency = 'high'
        elif poor_factors <= 1:
            overall_efficiency = 'medium'
        else:
            overall_efficiency = 'low'
        
        return {
            'overall_efficiency': overall_efficiency,
            'efficiency_factors': efficiency_factors,
            'recommendations': recommendations
        }

class PCRSimulationReport:
    """
    Generates comprehensive reports from in-silico PCR simulations
    """
    
    def __init__(self):
        self.logger = logging.getLogger(__name__)
    
    def generate_report(self, amplicons, specificity_analysis, efficiency_predictions):
        """
        Generate a comprehensive PCR simulation report
        
        Args:
            amplicons (list): List of amplicon results
            specificity_analysis (dict): Specificity analysis results
            efficiency_predictions (dict): Efficiency prediction results
            
        Returns:
            dict: Comprehensive report
        """
        return {
            'summary': {
                'total_amplicons': len(amplicons),
                'specificity_score': specificity_analysis.get('specificity_score', 0),
                'efficiency_rating': efficiency_predictions.get('overall_efficiency', 'unknown')
            },
            'amplicon_details': self._format_amplicon_details(amplicons),
            'specificity_analysis': specificity_analysis,
            'efficiency_analysis': efficiency_predictions,
            'recommendations': self._generate_recommendations(
                specificity_analysis, efficiency_predictions
            )
        }
    
    def _format_amplicon_details(self, amplicons):
        """Format amplicon details for reporting"""
        if not amplicons:
            return []
        
        return [
            {
                'template_id': amp.template_id,
                'size': amp.amplicon_size,
                'position': f"{amp.forward_start}-{amp.reverse_end}",
                'sequence_preview': amp.amplicon_sequence[:50] + "..." if len(amp.amplicon_sequence) > 50 else amp.amplicon_sequence
            }
            for amp in amplicons
        ]
    
    def _summarize_efficiency(self, efficiency_predictions):
        """Summarize efficiency predictions"""
        factors = efficiency_predictions.get('efficiency_factors', [])
        
        return {
            'overall_rating': efficiency_predictions.get('overall_efficiency', 'unknown'),
            'factor_count': len(factors),
            'good_factors': len([f for f in factors if f['status'] == 'good']),
            'poor_factors': len([f for f in factors if f['status'] == 'poor'])
        }
    
    def _generate_recommendations(self, specificity_analysis, efficiency_predictions):
        """Generate recommendations based on analysis results"""
        recommendations = []
        
        # Specificity-based recommendations
        specificity_score = specificity_analysis.get('specificity_score', 0)
        if specificity_score < 0.5:
            recommendations.append("Low specificity detected. Consider redesigning primers.")
        elif specificity_score < 0.8:
            recommendations.append("Moderate specificity. Validate with experimental testing.")
        
        # Efficiency-based recommendations
        efficiency_recs = efficiency_predictions.get('recommendations', [])
        recommendations.extend(efficiency_recs)
        
        # General recommendations
        if not recommendations:
            recommendations.append("Primer pair appears suitable for PCR amplification.")
        
        return recommendations