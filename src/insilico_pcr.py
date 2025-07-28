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
                binding_site = {
                    'position': i,
                    'end_position': i + primer_len - 1,
                    'mismatches': mismatches,
                    'sequence': template_segment,
                    'tm': mt.Tm_NN(template_segment),
                    'gc_content': calculate_gc_content(template_segment),
                    'direction': direction
                }
                binding_sites.append(binding_site)
        
        return binding_sites
    
    def _count_mismatches(self, seq1, seq2):
        """
        Count mismatches between two sequences
        
        Args:
            seq1 (str): First sequence
            seq2 (str): Second sequence
            
        Returns:
            int: Number of mismatches
        """
        if len(seq1) != len(seq2):
            return float('inf')
        
        return sum(1 for a, b in zip(seq1.upper(), seq2.upper()) if a != b)
    
    def _generate_amplicons(self, template, forward_sites, reverse_sites, template_id):
        """
        Generate amplicons from primer binding sites
        
        Args:
            template (str): Template sequence
            forward_sites (list): Forward primer binding sites
            reverse_sites (list): Reverse primer binding sites
            template_id (str): Template identifier
            
        Returns:
            list: Generated amplicons
        """
        amplicons = []
        
        for f_site in forward_sites:
            for r_site in reverse_sites:
                # Check if primers are oriented correctly
                if f_site['position'] < r_site['position']:
                    # Forward primer upstream of reverse primer
                    amplicon_start = f_site['position']
                    amplicon_end = r_site['end_position']
                    amplicon_size = amplicon_end - amplicon_start + 1
                    
                    # Extract amplicon sequence
                    amplicon_seq = template[amplicon_start:amplicon_end + 1]
                    
                    amplicon = AmpliconResult(
                        forward_start=f_site['position'],
                        forward_end=f_site['end_position'],
                        reverse_start=r_site['position'],
                        reverse_end=r_site['end_position'],
                        amplicon_size=amplicon_size,
                        amplicon_sequence=amplicon_seq,
                        template_id=template_id,
                        strand='+'
                    )
                    amplicons.append(amplicon)
        
        return amplicons
    
    def analyze_amplicon_specificity(self, amplicons, target_template_ids=None):
        """
        Analyze amplicon specificity and identify non-specific products
        
        Args:
            amplicons (list): List of amplicons
            target_template_ids (list): Expected target template IDs
            
        Returns:
            dict: Specificity analysis results
        """
        # Group amplicons by size
        size_groups = {}
        for amplicon in amplicons:
            size = amplicon.amplicon_size
            if size not in size_groups:
                size_groups[size] = []
            size_groups[size].append(amplicon)
        
        # Identify target vs non-target amplicons
        target_amplicons = []
        non_target_amplicons = []
        
        for amplicon in amplicons:
            if target_template_ids and amplicon.template_id in target_template_ids:
                target_amplicons.append(amplicon)
            else:
                non_target_amplicons.append(amplicon)
        
        # Calculate specificity metrics
        total_amplicons = len(amplicons)
        target_count = len(target_amplicons)
        non_target_count = len(non_target_amplicons)
        
        specificity_ratio = (target_count / total_amplicons * 100) if total_amplicons > 0 else 0
        
        analysis = {
            'total_amplicons': total_amplicons,
            'target_amplicons': target_count,
            'non_target_amplicons': non_target_count,
            'specificity_percentage': specificity_ratio,
            'size_distribution': self._calculate_size_distribution(amplicons),
            'non_specific_products': self._identify_problematic_products(
                non_target_amplicons, target_amplicons
            )
        }
        
        return analysis
    
    def _calculate_size_distribution(self, amplicons):
        """
        Calculate size distribution of amplicons
        
        Args:
            amplicons (list): List of amplicons
            
        Returns:
            dict: Size distribution statistics
        """
        if not amplicons:
            return {}
        
        sizes = [amp.amplicon_size for amp in amplicons]
        
        return {
            'mean_size': np.mean(sizes),
            'median_size': np.median(sizes),
            'std_size': np.std(sizes),
            'min_size': min(sizes),
            'max_size': max(sizes),
            'size_range': max(sizes) - min(sizes)
        }
    
    def _identify_problematic_products(self, non_target_amplicons, target_amplicons):
        """
        Identify potentially problematic non-specific products
        
        Args:
            non_target_amplicons (list): Non-target amplicons
            target_amplicons (list): Target amplicons
            
        Returns:
            list: Problematic products with risk assessment
        """
        if not target_amplicons:
            return []
        
        target_sizes = [amp.amplicon_size for amp in target_amplicons]
        target_size_range = (min(target_sizes), max(target_sizes))
        
        problematic_products = []
        
        for amplicon in non_target_amplicons:
            risk_level = "Low"
            risk_factors = []
            
            # Check size similarity to target products
            if target_size_range[0] <= amplicon.amplicon_size <= target_size_range[1]:
                risk_level = "High"
                risk_factors.append("Size overlaps with target products")
            elif abs(amplicon.amplicon_size - np.mean(target_sizes)) < 100:
                risk_level = "Medium"
                risk_factors.append("Size close to target products")
            
            # Check for abundant templates
            if "chromosome" in amplicon.template_id.lower():
                risk_factors.append("Chromosomal template")
                if risk_level == "Low":
                    risk_level = "Medium"
            
            product_info = {
                'template_id': amplicon.template_id,
                'size': amplicon.amplicon_size,
                'risk_level': risk_level,
                'risk_factors': risk_factors
            }
            
            problematic_products.append(product_info)
        
        # Sort by risk level
        risk_order = {"High": 3, "Medium": 2, "Low": 1}
        problematic_products.sort(
            key=lambda x: risk_order.get(x['risk_level'], 0), 
            reverse=True
        )
        
        return problematic_products
    
    def predict_pcr_efficiency(self, amplicons, primer_tm_forward, primer_tm_reverse):
        """
        Predict PCR efficiency based on amplicon characteristics
        
        Args:
            amplicons (list): List of amplicons
            primer_tm_forward (float): Forward primer melting temperature
            primer_tm_reverse (float): Reverse primer melting temperature
            
        Returns:
            dict: Efficiency predictions for each amplicon
        """
        efficiency_predictions = []
        
        for amplicon in amplicons:
            # Calculate GC content of amplicon
            gc_content = calculate_gc_content(amplicon.amplicon_sequence)
            
            # Estimate amplicon Tm
            amplicon_tm = mt.Tm_NN(amplicon.amplicon_sequence)
            
            # Predict efficiency based on multiple factors
            efficiency_score = 100  # Start with perfect efficiency
            
            # Size penalty
            if amplicon.amplicon_size > 1000:
                efficiency_score -= (amplicon.amplicon_size - 1000) * 0.01
            
            # GC content penalty
            optimal_gc = 50
            gc_penalty = abs(gc_content - optimal_gc) * 0.5
            efficiency_score -= gc_penalty
            
            # Tm difference penalty
            avg_primer_tm = (primer_tm_forward + primer_tm_reverse) / 2
            tm_diff_penalty = abs(amplicon_tm - avg_primer_tm) * 0.1
            efficiency_score -= tm_diff_penalty
            
            # Ensure efficiency is between 0 and 100
            efficiency_score = max(0, min(100, efficiency_score))
            
            prediction = {
                'amplicon': amplicon,
                'efficiency_score': efficiency_score,
                'gc_content': gc_content,
                'amplicon_tm': amplicon_tm,
                'size_factor': min(1.0, 1000 / amplicon.amplicon_size),
                'predicted_yield': 'High' if efficiency_score > 80 else 
                                 'Medium' if efficiency_score > 60 else 'Low'
            }
            
            efficiency_predictions.append(prediction)
        
        return efficiency_predictions


class PCRSimulationReport:
    """
    Generates detailed reports for PCR simulation results
    """
    
    def __init__(self):
        self.logger = logging.getLogger(__name__)
    
    def generate_report(self, amplicons, specificity_analysis, efficiency_predictions):
        """
        Generate comprehensive PCR simulation report
        
        Args:
            amplicons (list): List of amplicons
            specificity_analysis (dict): Specificity analysis results
            efficiency_predictions (list): Efficiency predictions
            
        Returns:
            dict: Comprehensive report
        """
        report = {
            'summary': {
                'total_products': len(amplicons),
                'expected_products': specificity_analysis.get('target_amplicons', 0),
                'non_specific_products': specificity_analysis.get('non_target_amplicons', 0),
                'specificity_score': specificity_analysis.get('specificity_percentage', 0)
            },
            'product_details': self._format_amplicon_details(amplicons),
            'size_analysis': specificity_analysis.get('size_distribution', {}),
            'efficiency_analysis': self._summarize_efficiency(efficiency_predictions),
            'recommendations': self._generate_recommendations(
                specificity_analysis, efficiency_predictions
            )
        }
        
        return report
    
    def _format_amplicon_details(self, amplicons):
        """
        Format amplicon details for report
        
        Args:
            amplicons (list): List of amplicons
            
        Returns:
            list: Formatted amplicon details
        """
        details = []
        
        for i, amplicon in enumerate(amplicons, 1):
            detail = {
                'product_id': f"Product_{i}",
                'template': amplicon.template_id,
                'size': amplicon.amplicon_size,
                'forward_position': f"{amplicon.forward_start}-{amplicon.forward_end}",
                'reverse_position': f"{amplicon.reverse_start}-{amplicon.reverse_end}",
                'gc_content': round(calculate_gc_content(amplicon.amplicon_sequence), 2)
            }
            details.append(detail)
        
        return details
    
    def _summarize_efficiency(self, efficiency_predictions):
        """
        Summarize efficiency predictions
        
        Args:
            efficiency_predictions (list): Efficiency predictions
            
        Returns:
            dict: Efficiency summary
        """
        if not efficiency_predictions:
            return {}
        
        scores = [pred['efficiency_score'] for pred in efficiency_predictions]
        
        return {
            'average_efficiency': np.mean(scores),
            'efficiency_range': f"{min(scores):.1f} - {max(scores):.1f}",
            'high_efficiency_products': len([s for s in scores if s > 80]),
            'medium_efficiency_products': len([s for s in scores if 60 <= s <= 80]),
            'low_efficiency_products': len([s for s in scores if s < 60])
        }
    
    def _generate_recommendations(self, specificity_analysis, efficiency_predictions):
        """
        Generate recommendations based on analysis
        
        Args:
            specificity_analysis (dict): Specificity analysis
            efficiency_predictions (list): Efficiency predictions
            
        Returns:
            list: List of recommendations
        """
        recommendations = []
        
        # Specificity recommendations
        specificity_score = specificity_analysis.get('specificity_percentage', 0)
        if specificity_score < 70:
            recommendations.append(
                "Low specificity detected. Consider primer redesign or "
                "optimization of PCR conditions."
            )
        elif specificity_score < 90:
            recommendations.append(
                "Moderate specificity. Validate results with gel electrophoresis."
            )
        
        # Efficiency recommendations
        if efficiency_predictions:
            avg_efficiency = np.mean([p['efficiency_score'] for p in efficiency_predictions])
            if avg_efficiency < 60:
                recommendations.append(
                    "Low predicted efficiency. Consider optimizing annealing temperature "
                    "or primer concentrations."
                )
        
        # Size recommendations
        size_dist = specificity_analysis.get('size_distribution', {})
        if size_dist.get('size_range', 0) > 1000:
            recommendations.append(
                "Wide size range of products. Consider more specific primers."
            )
        
        if not recommendations:
            recommendations.append("PCR design looks good. Proceed with experimental validation.")
        
        return recommendations


# Example usage
if __name__ == "__main__":
    # Configure logging
    logging.basicConfig(level=logging.INFO)
    
    # Example template sequences
    templates = {
        "target_gene": "ATGGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGC",
        "off_target_1": "ATGGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGT"
    }
    
    # Example primers
    forward_primer = "ATGGCTAGCTAGCTAGCTAG"
    reverse_primer = "CTAGCTAGCTAGCTAGCCAT"
    
    # Run simulation
    pcr_sim = InSilicoPCR()
    amplicons = pcr_sim.simulate_pcr(templates, forward_primer, reverse_primer)
    
    # Analyze results
    specificity = pcr_sim.analyze_amplicon_specificity(amplicons, ["target_gene"])
    efficiency = pcr_sim.predict_pcr_efficiency(amplicons, 60.0, 58.0)
    
    # Generate report
    reporter = PCRSimulationReport()
    report = reporter.generate_report(amplicons, specificity, efficiency)
    
    print("PCR Simulation Report:")
    print(f"Total products: {report['summary']['total_products']}")
    print(f"Specificity: {report['summary']['specificity_score']:.1f}%")
    print(f"Recommendations: {report['recommendations']}")