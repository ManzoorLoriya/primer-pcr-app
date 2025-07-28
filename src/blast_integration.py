import requests
import time
import xml.etree.ElementTree as ET
from Bio.Blast import NCBIWWW, NCBIXML
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import logging

class BlastIntegration:
    """
    Handles BLAST searches for primer specificity checking
    """
    
    def __init__(self):
        self.base_url = "https://blast.ncbi.nlm.nih.gov/Blast.cgi"
        self.databases = {
            'nt': 'Nucleotide collection (nr/nt)',
            'refseq_rna': 'RefSeq RNA sequences',
            'refseq_genomic': 'RefSeq Genome Database',
            'human_genomic': 'Human Genome Database',
            'mouse_genomic': 'Mouse Genome Database'
        }
        self.logger = logging.getLogger(__name__)
    
    def blast_primer(self, primer_sequence, database='nt', organism=None, 
                    max_hits=100, expect_threshold=10):
        """
        Perform BLAST search for a primer sequence
        
        Args:
            primer_sequence (str): Primer sequence to search
            database (str): BLAST database to search against
            organism (str): Organism to limit search (optional)
            max_hits (int): Maximum number of hits to return
            expect_threshold (float): E-value threshold
            
        Returns:
            list: BLAST results with alignment information
        """
        try:
            # Submit BLAST search
            self.logger.info(f"Submitting BLAST search for primer: {primer_sequence[:20]}...")
            result_handle = NCBIWWW.qblast(
                program="blastn",
                database=database,
                sequence=primer_sequence,
                expect=expect_threshold,
                hitlist_size=max_hits,
                entrez_query=f'{organism}[ORGN]' if organism else None,
                format_type="XML"
            )
            
            # Parse results
            blast_records = NCBIXML.parse(result_handle)
            blast_record = next(blast_records)
            
            # Process alignments
            alignments = []
            for alignment in blast_record.alignments:
                for hsp in alignment.hsps:
                    alignment_data = {
                        'title': alignment.title,
                        'accession': alignment.accession,
                        'length': alignment.length,
                        'e_value': hsp.expect,
                        'bit_score': hsp.bits,
                        'score': hsp.score,
                        'identities': hsp.identities,
                        'gaps': hsp.gaps,
                        'align_length': hsp.align_length,
                        'query_start': hsp.query_start,
                        'query_end': hsp.query_end,
                        'subject_start': hsp.sbjct_start,
                        'subject_end': hsp.sbjct_end,
                        'query_seq': hsp.query,
                        'subject_seq': hsp.sbjct,
                        'match_seq': hsp.match,
                        'percent_identity': (hsp.identities / hsp.align_length) * 100
                    }
                    alignments.append(alignment_data)
            
            self.logger.info(f"Found {len(alignments)} alignments")
            return alignments
            
        except Exception as e:
            self.logger.error(f"BLAST search failed: {str(e)}")
            return []
    
    def batch_blast_primers(self, primer_list, database='nt', organism=None):
        """
        Perform BLAST searches for multiple primers
        
        Args:
            primer_list (list): List of primer sequences
            database (str): BLAST database
            organism (str): Organism filter
            
        Returns:
            dict: Results for each primer
        """
        results = {}
        
        for i, primer in enumerate(primer_list):
            self.logger.info(f"Processing primer {i+1}/{len(primer_list)}")
            
            # Add delay to respect NCBI rate limits
            if i > 0:
                time.sleep(3)
            
            results[primer] = self.blast_primer(
                primer, database, organism
            )
            
        return results
    
    def filter_off_targets(self, blast_results, identity_threshold=80, 
                          coverage_threshold=80):
        """
        Filter BLAST results to identify potential off-targets
        
        Args:
            blast_results (list): BLAST alignment results
            identity_threshold (float): Minimum identity % for off-target
            coverage_threshold (float): Minimum coverage % for off-target
            
        Returns:
            list: Filtered off-target results
        """
        off_targets = []
        
        for result in blast_results:
            # Calculate coverage
            query_coverage = ((result['query_end'] - result['query_start'] + 1) / 
                            len(result['query_seq'])) * 100
            
            # Filter based on thresholds
            if (result['percent_identity'] >= identity_threshold and 
                query_coverage >= coverage_threshold):
                
                result['query_coverage'] = query_coverage
                result['risk_level'] = self._assess_risk_level(
                    result['percent_identity'], query_coverage
                )
                off_targets.append(result)
        
        # Sort by risk level (highest first)
        off_targets.sort(key=lambda x: x['percent_identity'], reverse=True)
        
        return off_targets
    
    def _assess_risk_level(self, identity, coverage):
        """
        Assess the risk level of off-target binding
        
        Args:
            identity (float): Percent identity
            coverage (float): Query coverage
            
        Returns:
            str: Risk level (High, Medium, Low)
        """
        if identity >= 95 and coverage >= 90:
            return "High"
        elif identity >= 85 and coverage >= 80:
            return "Medium"
        else:
            return "Low"
    
    def get_organism_databases(self):
        """
        Return available organism-specific databases
        
        Returns:
            dict: Database options
        """
        return self.databases


class SpecificityAnalyzer:
    """
    Analyzes primer specificity and provides recommendations
    """
    
    def __init__(self):
        self.blast_client = BlastIntegration()
        self.logger = logging.getLogger(__name__)
    
    def analyze_primer_pair(self, forward_primer, reverse_primer, 
                           target_organism=None, database='nt'):
        """
        Analyze specificity for a primer pair
        
        Args:
            forward_primer (str): Forward primer sequence
            reverse_primer (str): Reverse primer sequence
            target_organism (str): Target organism
            database (str): BLAST database
            
        Returns:
            dict: Comprehensive specificity analysis
        """
        self.logger.info("Starting primer pair specificity analysis")
        
        # Analyze forward primer
        forward_results = self.blast_client.blast_primer(
            forward_primer, database, target_organism
        )
        forward_off_targets = self.blast_client.filter_off_targets(forward_results)
        
        # Analyze reverse primer
        reverse_results = self.blast_client.blast_primer(
            reverse_primer, database, target_organism
        )
        reverse_off_targets = self.blast_client.filter_off_targets(reverse_results)
        
        # Combine analysis
        analysis = {
            'forward_primer': {
                'sequence': forward_primer,
                'total_hits': len(forward_results),
                'off_targets': forward_off_targets,
                'specificity_score': self._calculate_specificity_score(forward_off_targets)
            },
            'reverse_primer': {
                'sequence': reverse_primer,
                'total_hits': len(reverse_results),
                'off_targets': reverse_off_targets,
                'specificity_score': self._calculate_specificity_score(reverse_off_targets)
            }
        }
        
        # Overall assessment
        analysis['overall_specificity'] = self._assess_overall_specificity(
            analysis['forward_primer']['specificity_score'],
            analysis['reverse_primer']['specificity_score']
        )
        
        return analysis
    
    def _calculate_specificity_score(self, off_targets):
        """
        Calculate specificity score based on off-targets
        
        Args:
            off_targets (list): List of off-target hits
            
        Returns:
            float: Specificity score (0-100)
        """
        if not off_targets:
            return 100.0
        
        # Weight off-targets by risk level
        penalty = 0
        for target in off_targets:
            if target['risk_level'] == 'High':
                penalty += 20
            elif target['risk_level'] == 'Medium':
                penalty += 10
            else:
                penalty += 5
        
        # Cap at 100
        score = max(0, 100 - penalty)
        return score
    
    def _assess_overall_specificity(self, forward_score, reverse_score):
        """
        Assess overall primer pair specificity
        
        Args:
            forward_score (float): Forward primer specificity score
            reverse_score (float): Reverse primer specificity score
            
        Returns:
            dict: Overall assessment
        """
        avg_score = (forward_score + reverse_score) / 2
        
        if avg_score >= 90:
            assessment = "Excellent"
        elif avg_score >= 80:
            assessment = "Good"
        elif avg_score >= 70:
            assessment = "Fair"
        else:
            assessment = "Poor"
        
        return {
            'score': avg_score,
            'assessment': assessment,
            'recommendation': self._get_recommendation(assessment)
        }
    
    def _get_recommendation(self, assessment):
        """
        Provide recommendations based on specificity assessment
        
        Args:
            assessment (str): Specificity assessment
            
        Returns:
            str: Recommendation text
        """
        recommendations = {
            "Excellent": "Primer pair shows excellent specificity. Proceed with confidence.",
            "Good": "Primer pair shows good specificity. Consider validation with gel electrophoresis.",
            "Fair": "Primer pair shows fair specificity. Optimize PCR conditions or consider redesign.",
            "Poor": "Primer pair shows poor specificity. Redesign recommended."
        }
        
        return recommendations.get(assessment, "Unable to assess specificity.")


# Example usage
if __name__ == "__main__":
    # Configure logging
    logging.basicConfig(level=logging.INFO)
    
    # Example primers
    forward_primer = "ATGGCTAGCTAGCTAGCTAG"
    reverse_primer = "CTAGCTAGCTAGCTAGCCAT"
    
    # Analyze specificity
    analyzer = SpecificityAnalyzer()
    results = analyzer.analyze_primer_pair(
        forward_primer, 
        reverse_primer, 
        target_organism="Homo sapiens"
    )
    
    print("Specificity Analysis Results:")
    print(f"Forward Primer Specificity: {results['forward_primer']['specificity_score']}")
    print(f"Reverse Primer Specificity: {results['reverse_primer']['specificity_score']}")
    print(f"Overall Assessment: {results['overall_specificity']['assessment']}")
    print(f"Recommendation: {results['overall_specificity']['recommendation']}")