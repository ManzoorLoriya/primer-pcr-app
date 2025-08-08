"""
Primer3 integration and primer design functionality
"""

import primer3
from typing import Dict, List, Optional, Tuple
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from .algorithms import PrimerAnalyzer

class PrimerDesigner:
    """Main class for primer design using Primer3 with custom analysis"""
    
    def __init__(self):
        self.analyzer = PrimerAnalyzer()
        self.default_params = {
            'PRIMER_TASK': 'generic',
            'PRIMER_PICK_LEFT_PRIMER': 1,
            'PRIMER_PICK_INTERNAL_OLIGO': 0,
            'PRIMER_PICK_RIGHT_PRIMER': 1,
            'PRIMER_OPT_SIZE': 20,
            'PRIMER_MIN_SIZE': 18,
            'PRIMER_MAX_SIZE': 25,
            'PRIMER_OPT_TM': 60.0,
            'PRIMER_MIN_TM': 57.0,
            'PRIMER_MAX_TM': 63.0,
            'PRIMER_MIN_GC': 20.0,
            'PRIMER_MAX_GC': 80.0,
            'PRIMER_OPT_GC_PERCENT': 50.0,
            'PRIMER_MAX_POLY_X': 4,
            'PRIMER_SALT_MONOVALENT': 50.0,
            'PRIMER_DNA_CONC': 50.0,
            'PRIMER_MAX_NS_ACCEPTED': 0,
            'PRIMER_MAX_SELF_ANY': 12.0,
            'PRIMER_MAX_SELF_END': 8.0,
            'PRIMER_PAIR_MAX_COMPL_ANY': 12.0,
            'PRIMER_PAIR_MAX_COMPL_END': 8.0,
            'PRIMER_PRODUCT_SIZE_RANGE': [[100, 1000]]
        }
    
    def design_primers(self, sequence: str, target_region: Optional[Tuple[int, int]] = None,
                      custom_params: Optional[Dict] = None, num_return: int = 5) -> Dict:
        """
        Design primers for a given sequence using Primer3
        
        Args:
            sequence: Target DNA sequence
            target_region: Optional tuple (start, length) for specific targeting
            custom_params: Custom primer3 parameters
            num_return: Number of primer pairs to return
            
        Returns:
            Dictionary containing primer design results
        """
        # Prepare parameters
        params = self.default_params.copy()
        if custom_params:
            params.update(custom_params)
        
        params['PRIMER_NUM_RETURN'] = num_return
        
        # Set up sequence and target
        seq_args = {
            'SEQUENCE_ID': 'target_sequence',
            'SEQUENCE_TEMPLATE': sequence.upper()
        }
        
        if target_region:
            seq_args['SEQUENCE_TARGET'] = [target_region[0], target_region[1]]
        
        try:
            # Run Primer3
            primer_results = primer3.bindings.design_primers(seq_args, params)
            
            # Process and enhance results
            processed_results = self._process_primer3_results(
                primer_results, sequence, params
            )
            
            return processed_results
            
        except Exception as e:
            return {
                'success': False,
                'error': f"Primer design failed: {str(e)}",
                'primer_pairs': []
            }
    
    def _process_primer3_results(self, results: Dict, sequence: str, 
                                params: Dict) -> Dict:
        """Process and enhance Primer3 results with custom analysis"""
        
        processed = {
            'success': True,
            'num_pairs_found': results.get('PRIMER_PAIR_NUM_RETURNED', 0),
            'primer_pairs': [],
            'sequence_info': {
                'length': len(sequence),
                'gc_content': self.analyzer.calculate_gc_content(sequence)
            }
        }
        
        # Extract primer pairs
        for i in range(processed['num_pairs_found']):
            pair_data = self._extract_primer_pair_data(results, i, sequence)
            if pair_data:
                processed['primer_pairs'].append(pair_data)
        
        return processed
    
    def _extract_primer_pair_data(self, results: Dict, pair_index: int, 
                                 sequence: str) -> Optional[Dict]:
        """Extract and analyze a single primer pair from Primer3 results"""
        
        try:
            # Extract primer sequences
            left_seq = results.get(f'PRIMER_LEFT_{pair_index}_SEQUENCE', '')
            right_seq = results.get(f'PRIMER_RIGHT_{pair_index}_SEQUENCE', '')
            
            if not left_seq or not right_seq:
                return None
            
            # Get positions
            left_pos = results.get(f'PRIMER_LEFT_{pair_index}', [0, 0])[0]
            right_pos = results.get(f'PRIMER_RIGHT_{pair_index}', [0, 0])[0]
            
            # Analyze individual primers
            left_analysis = self.analyzer.analyze_primer(left_seq)
            right_analysis = self.analyzer.analyze_primer(right_seq)
            
            # Calculate Tm difference
            tm_diff = abs(left_analysis['melting_temperature'] - right_analysis['melting_temperature'])
            
            # Analyze dimer formation
            dimer_analysis = self.analyzer.analyze_dimer_formation(left_seq, right_seq)
            
            # Calculate overall pair quality
            pair_quality = self._calculate_pair_quality(
                left_analysis, right_analysis, tm_diff, dimer_analysis
            )
            
            return {
                'pair_index': pair_index,
                'left_primer': {
                    'sequence': left_seq,
                    'position': left_pos,
                    'analysis': left_analysis
                },
                'right_primer': {
                    'sequence': right_seq,
                    'position': right_pos,
                    'analysis': right_analysis
                },
                'product_size': right_pos - left_pos + len(right_seq),
                'tm_difference': tm_diff,
                'dimer_analysis': dimer_analysis,
                'pair_quality_score': pair_quality
            }
            
        except Exception as e:
            return None
    
    def _calculate_pair_quality(self, left_analysis: Dict, right_analysis: Dict,
                              tm_difference: float, dimer_analysis: List) -> float:
        """Calculate overall quality score for a primer pair"""
        
        # Base quality from individual primers
        left_quality = left_analysis['quality_scores']['overall_score']
        right_quality = right_analysis['quality_scores']['overall_score']
        
        # Tm difference penalty
        tm_penalty = max(0, (tm_difference - 3) * 5)
        
        # Dimer penalty
        dimer_penalty = sum(d['risk_score'] for d in dimer_analysis if d['risk_level'] == 'HIGH')
        
        # Calculate weighted average
        base_quality = (left_quality + right_quality) / 2
        final_quality = max(0, base_quality - tm_penalty - dimer_penalty)
        
        return round(final_quality, 1)
    
    def batch_design_primers(self, sequences: List[Dict], 
                          common_params: Optional[Dict] = None) -> List[Dict]:
        """Design primers for multiple sequences"""
        
        results = []
        for seq_data in sequences:
            result = self.design_primers(
                seq_data['sequence'],
                target_region=seq_data.get('target_region'),
                custom_params=common_params,
                num_return=seq_data.get('num_return', 5)
            )
            result['sequence_id'] = seq_data.get('id', 'unknown')
            results.append(result)
        
        return results
    
    def optimize_primers(self, sequence: str, failed_attempts: List[Dict]) -> Dict:
        """Optimize primer design based on previous failures"""
        
        # Analyze failure patterns
        common_issues = self._analyze_failures(failed_attempts)
        
        # Adjust parameters based on issues
        optimized_params = self.default_params.copy()
        
        if 'high_tm' in common_issues:
            optimized_params['PRIMER_MAX_TM'] = 58.0
            optimized_params['PRIMER_OPT_TM'] = 55.0
        
        if 'low_gc' in common_issues:
            optimized_params['PRIMER_MIN_GC'] = 30.0
            optimized_params['PRIMER_OPT_GC_PERCENT'] = 45.0
        
        if 'dimer_issues' in common_issues:
            optimized_params['PRIMER_MAX_SELF_ANY'] = 8.0
            optimized_params['PRIMER_MAX_SELF_END'] = 4.0
        
        # Try design with optimized parameters
        return self.design_primers(sequence, custom_params=optimized_params)
    
    def _analyze_failures(self, failed_attempts: List[Dict]) -> List[str]:
        """Analyze common issues in failed primer designs"""
        
        issues = []
        for attempt in failed_attempts:
            if attempt.get('error', '').lower().find('tm') != -1:
                issues.append('high_tm')
            if attempt.get('error', '').lower().find('gc') != -1:
                issues.append('low_gc')
            if attempt.get('error', '').lower().find('dimer') != -1:
                issues.append('dimer_issues')
        
        return list(set(issues))
    
    def get_primer_recommendations(self, analysis_result: Dict) -> List[str]:
        """Generate recommendations based on primer analysis"""
        
        recommendations = []
        
        if not analysis_result.get('success', False):
            recommendations.append("Primer design failed - check sequence quality and parameters")
            return recommendations
        
        if not analysis_result['primer_pairs']:
            recommendations.append("No suitable primer pairs found - try relaxing constraints")
            return recommendations
        
        # Analyze top pair
        top_pair = analysis_result['primer_pairs'][0]
        
        # Check Tm difference
        tm_diff = top_pair['tm_difference']
        if tm_diff > 5:
            recommendations.append(
                f"Large Tm difference ({tm_diff:.1f}°C) - consider adjusting primer lengths"
            )
        
        # Check dimer formation
        high_risk_dimers = [d for d in top_pair['dimer_analysis'] if d['risk_level'] == 'HIGH']
        if high_risk_dimers:
            recommendations.append(
                f"Found {len(high_risk_dimers)} high-risk primer dimers - consider redesigning"
            )
        
        # Check individual primer quality
        left_quality = top_pair['left_primer']['analysis']['quality_scores']['overall_score']
        right_quality = top_pair['right_primer']['analysis']['quality_scores']['overall_score']
        
        if left_quality < 70:
            recommendations.append("Left primer quality is suboptimal - check GC content and Tm")
        if right_quality < 70:
            recommendations.append("Right primer quality is suboptimal - check GC content and Tm")
        
        if not recommendations:
            recommendations.append("Primer design looks good - ready for validation!")
        
        return recommendations
    
    def multiplex_design(self, sequences, params=None):
        """Design primers for multiple targets simultaneously"""
        multiplex_params = params or self.default_params.copy()
        multiplex_params.update({
            'PRIMER_PICK_ANYWAY': 1,
            'PRIMER_MAX_TM_DIFF': 3.0,
            'PRIMER_PAIR_WT_TM_DIFF': 0.5
        })
        
        results = {}
        for seq_id, sequence in sequences.items():
            results[seq_id] = self.design_primers(sequence, custom_params=multiplex_params)
        
        # Cross-check specificity
        all_primers = [p for res in results.values() for pair in res['primer_pairs'] 
                      for p in [pair['left_primer']['sequence'], pair['right_primer']['sequence']]]
        
        for seq_id in results:
            for pair in results[seq_id]['primer_pairs']:
                cross_reactivity = self.analyzer.check_cross_reactivity(
                    pair['left_primer']['sequence'],
                    pair['right_primer']['sequence'],
                    all_primers
                )
                pair['cross_reactivity'] = cross_reactivity
        
        return results

    def optimize_primers(self, sequence, initial_results):
        """Provide optimization suggestions"""
        recommendations = []
        top_pair = initial_results['primer_pairs'][0]
        
        # Tm optimization
        tm_diff = abs(top_pair['left_primer']['analysis']['melting_temperature'] - 
                     top_pair['right_primer']['analysis']['melting_temperature'])
        if tm_diff > 3:
            recommendations.append(
                f"Adjust primer lengths to reduce Tm difference (current ΔTm={tm_diff:.1f}°C)"
            )
        
        # GC optimization
        for side in ['left_primer', 'right_primer']:
            gc = top_pair[side]['analysis']['gc_content']
            if gc < 40 or gc > 60:
                recommendations.append(
                    f"Optimize GC content of {side.replace('_', ' ')} (current: {gc:.1f}%)"
                )
        
        # Dimer resolution
        if any(d['risk_level'] == 'HIGH' for d in top_pair['dimer_analysis']):
            recommendations.append("Redesign primers to eliminate high-risk dimers")
        
        return recommendations

    def calculate_synthesis_cost(self, primer_sequence, scale='25nm', purification='standard'):
        """Calculate synthesis cost based on commercial pricing"""
        base_costs = {
            '25nm': 0.25,
            '100nm': 0.35,
            '250nm': 0.50
        }
        purification_costs = {
            'standard': 0,
            'hplc': 15,
            'page': 25
        }
        
        length = len(primer_sequence)
        return round(
            base_costs.get(scale, 0.35) * length + purification_costs.get(purification, 0),
            2
        )

    def search_primer_databases(self, sequence):
        """Check commercial databases for existing primers"""
        # Placeholder for actual API integrations
        return {
            'IDT': f"https://www.idtdna.com/pages/tools/analyze?sequence={sequence}",
            'ThermoFisher': f"https://www.thermofisher.com/order/oligo-builder/input?sequence={sequence}"
        }