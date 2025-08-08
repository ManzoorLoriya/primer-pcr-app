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
        
        for primer in primer_list:
            try:
                results[primer] = self.blast_primer(primer, database, organism)
                time.sleep(1)  # Rate limiting
            except Exception as e:
                self.logger.error(f"BLAST failed for primer {primer}: {str(e)}")
                results[primer] = []
        
        return results
    
    def filter_off_targets(self, blast_results, identity_threshold=80, 
                          coverage_threshold=80):
        """
        Filter BLAST results to identify potential off-targets
        
        Args:
            blast_results (list): BLAST alignment results
            identity_threshold (float): Minimum identity percentage
            coverage_threshold (float): Minimum coverage percentage
            
        Returns:
            list: Filtered off-target results
        """
        off_targets = []
        
        for result in blast_results:
            identity = result['percent_identity']
            coverage = (result['align_length'] / len(result['query_seq'])) * 100
            
            if identity >= identity_threshold and coverage >= coverage_threshold:
                risk_level = self._assess_risk_level(identity, coverage)
                result['risk_level'] = risk_level
                off_targets.append(result)
        
        return sorted(off_targets, key=lambda x: x['e_value'])
    
    def _assess_risk_level(self, identity, coverage):
        """Assess risk level based on identity and coverage"""
        if identity >= 95 and coverage >= 90:
            return 'high'
        elif identity >= 85 and coverage >= 80:
            return 'medium'
        else:
            return 'low'
    
    def get_organism_databases(self):
        """Get available organism-specific databases"""
        return self.databases

class SpecificityAnalyzer:
    """
    Analyzes primer specificity using BLAST results
    """
    
    def __init__(self):
        self.blast_integration = BlastIntegration()
    
    def analyze_primer_pair(self, forward_primer, reverse_primer, 
                           target_organism=None, database='nt'):
        """
        Analyze specificity of a primer pair
        
        Args:
            forward_primer (str): Forward primer sequence
            reverse_primer (str): Reverse primer sequence
            target_organism (str): Target organism
            database (str): BLAST database
            
        Returns:
            dict: Specificity analysis results
        """
        # BLAST both primers
        forward_results = self.blast_integration.blast_primer(
            forward_primer, database, target_organism
        )
        reverse_results = self.blast_integration.blast_primer(
            reverse_primer, database, target_organism
        )
        
        # Filter off-targets
        forward_off_targets = self.blast_integration.filter_off_targets(forward_results)
        reverse_off_targets = self.blast_integration.filter_off_targets(reverse_results)
        
        # Calculate specificity scores
        forward_score = self._calculate_specificity_score(forward_off_targets)
        reverse_score = self._calculate_specificity_score(reverse_off_targets)
        
        # Overall assessment
        overall_assessment = self._assess_overall_specificity(forward_score, reverse_score)
        
        return {
            'forward_primer': {
                'sequence': forward_primer,
                'off_targets': forward_off_targets,
                'specificity_score': forward_score
            },
            'reverse_primer': {
                'sequence': reverse_primer,
                'off_targets': reverse_off_targets,
                'specificity_score': reverse_score
            },
            'overall_assessment': overall_assessment,
            'recommendation': self._get_recommendation(overall_assessment)
        }
    
    def _calculate_specificity_score(self, off_targets):
        """Calculate specificity score based on off-targets"""
        if not off_targets:
            return 1.0  # Perfect specificity
        
        # Penalize based on number and quality of off-targets
        total_penalty = 0
        
        for off_target in off_targets:
            # Higher penalty for high-risk off-targets
            if off_target['risk_level'] == 'high':
                total_penalty += 0.3
            elif off_target['risk_level'] == 'medium':
                total_penalty += 0.15
            else:
                total_penalty += 0.05
            
            # Additional penalty for very high identity
            if off_target['percent_identity'] >= 98:
                total_penalty += 0.1
        
        return max(0.0, 1.0 - total_penalty)
    
    def _assess_overall_specificity(self, forward_score, reverse_score):
        """Assess overall specificity of primer pair"""
        avg_score = (forward_score + reverse_score) / 2
        
        if avg_score >= 0.9:
            return 'excellent'
        elif avg_score >= 0.7:
            return 'good'
        elif avg_score >= 0.5:
            return 'fair'
        else:
            return 'poor'
    
    def _get_recommendation(self, assessment):
        """Get recommendation based on specificity assessment"""
        recommendations = {
            'excellent': 'Primer pair has excellent specificity. Ready for use.',
            'good': 'Primer pair has good specificity. Minor optimization may be beneficial.',
            'fair': 'Primer pair has fair specificity. Consider redesigning primers.',
            'poor': 'Primer pair has poor specificity. Redesign recommended.'
        }
        return recommendations.get(assessment, 'Assessment unclear.')