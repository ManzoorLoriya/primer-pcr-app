import logging
from concurrent.futures import ThreadPoolExecutor
import time
from dataclasses import dataclass
from typing import List, Dict, Optional
from src.blast_integration import BlastIntegration, SpecificityAnalyzer
from src.insilico_pcr import InSilicoPCR, PCRSimulationReport
from Bio import SeqIO
from Bio.Seq import Seq
import json

@dataclass
class PrimerPair:
    """Data class for primer pair information"""
    forward: str
    reverse: str
    name: str = ""
    target_region: str = ""
    expected_size: int = 0

@dataclass
class SpecificityResult:
    """Data class for comprehensive specificity results"""
    primer_pair: PrimerPair
    blast_results: Dict
    pcr_simulation: Dict
    overall_score: float
    recommendations: List[str]
    processing_time: float

class ComprehensiveSpecificityChecker:
    """
    Integrated specificity checker combining BLAST and in-silico PCR
    """
    
    def __init__(self, max_threads=2):
        self.blast_analyzer = SpecificityAnalyzer()
        self.pcr_simulator = InSilicoPCR()
        self.report_generator = PCRSimulationReport()
        self.max_threads = max_threads
        self.logger = logging.getLogger(__name__)
        
        # Default organism databases
        self.organism_databases = {
            'human': {
                'blast_db': 'nt',
                'organism_filter': 'Homo sapiens',
                'common_templates': self._get_human_templates()
            },
            'mouse': {
                'blast_db': 'nt',
                'organism_filter': 'Mus musculus',
                'common_templates': self._get_mouse_templates()
            },
            'ecoli': {
                'blast_db': 'nt',
                'organism_filter': 'Escherichia coli',
                'common_templates': self._get_ecoli_templates()
            }
        }
    
    def check_specificity(self, primer_pair: PrimerPair, 
                         organism: str = 'human',
                         custom_templates: Optional[Dict[str, str]] = None,
                         blast_enabled: bool = True,
                         pcr_sim_enabled: bool = True) -> SpecificityResult:
        """
        Perform comprehensive specificity checking
        
        Args:
            primer_pair: PrimerPair object with forward/reverse sequences
            organism: Target organism ('human', 'mouse', 'ecoli', or 'custom')
            custom_templates: Custom template sequences for PCR simulation
            blast_enabled: Whether to perform BLAST analysis
            pcr_sim_enabled: Whether to perform PCR simulation
            
        Returns:
            SpecificityResult: Comprehensive analysis results
        """
        start_time = time.time()
        self.logger.info(f"Starting specificity check for {primer_pair.name}")
        
        # Initialize results
        blast_results = {}
        pcr_simulation = {}
        
        # Get organism configuration
        org_config = self.organism_databases.get(organism, {})
        templates = custom_templates or org_config.get('common_templates', {})
        
        # Run analyses
        if blast_enabled:
            self.logger.info("Running BLAST analysis...")
            blast_results = self._run_blast_analysis(primer_pair, org_config)
        
        if pcr_sim_enabled and templates:
            self.logger.info("Running PCR simulation...")
            pcr_simulation = self._run_pcr_simulation(primer_pair, templates)
        
        # Calculate overall score
        overall_score = self._calculate_overall_score(blast_results, pcr_simulation)
        
        # Generate recommendations
        recommendations = self._generate_comprehensive_recommendations(
            blast_results, pcr_simulation, overall_score
        )
        
        processing_time = time.time() - start_time
        
        result = SpecificityResult(
            primer_pair=primer_pair,
            blast_results=blast_results,
            pcr_simulation=pcr_simulation,
            overall_score=overall_score,
            recommendations=recommendations,
            processing_time=processing_time
        )
        
        self.logger.info(f"Specificity check completed in {processing_time:.2f} seconds")
        return result
    
    def batch_check_specificity(self, primer_pairs: List[PrimerPair],
                               organism: str = 'human',
                               custom_templates: Optional[Dict[str, str]] = None) -> List[SpecificityResult]:
        """
        Perform batch specificity checking for multiple primer pairs
        
        Args:
            primer_pairs: List of PrimerPair objects
            organism: Target organism
            custom_templates: Custom template sequences
            
        Returns:
            List[SpecificityResult]: Results for all primer pairs
        """
        self.logger.info(f"Starting batch specificity check for {len(primer_pairs)} primer pairs")
        
        results = []
        
        # Use ThreadPoolExecutor for parallel processing
        with ThreadPoolExecutor(max_workers=self.max_threads) as executor:
            # Submit all tasks
            futures = []
            for primer_pair in primer_pairs:
                future = executor.submit(
                    self.check_specificity,
                    primer_pair,
                    organism,
                    custom_templates
                )
                futures.append((future, primer_pair))
            
            # Collect results
            for future, primer_pair in futures:
                try:
                    result = future.result(timeout=300)  # 5 minute timeout
                    results.append(result)
                except Exception as e:
                    self.logger.error(f"Error processing {primer_pair.name}: {str(e)}")
                    # Create error result
                    error_result = SpecificityResult(
                        primer_pair=primer_pair,
                        blast_results={'error': str(e)},
                        pcr_simulation={'error': str(e)},
                        overall_score=0.0,
                        recommendations=[f"Analysis failed: {str(e)}"],
                        processing_time=0.0
                    )
                    results.append(error_result)
        
        return results
    
    def _run_blast_analysis(self, primer_pair: PrimerPair, org_config: Dict) -> Dict:
        """Run BLAST analysis for primer pair"""
        try:
            blast_db = org_config.get('blast_db', 'nt')
            organism_filter = org_config.get('organism_filter')
            
            analysis = self.blast_analyzer.analyze_primer_pair(
                primer_pair.forward,
                primer_pair.reverse,
                organism_filter,
                blast_db
            )
            
            return analysis
            
        except Exception as e:
            self.logger.error(f"BLAST analysis failed: {str(e)}")
            return {'error': str(e)}
    
    def _run_pcr_simulation(self, primer_pair: PrimerPair, templates: Dict[str, str]) -> Dict:
        """Run PCR simulation for primer pair"""
        try:
            # Simulate PCR
            amplicons = self.pcr_simulator.simulate_pcr(
                templates,
                primer_pair.forward,
                primer_pair.reverse
            )
            
            # Analyze specificity
            target_templates = [primer_pair.target_region] if primer_pair.target_region else None
            specificity_analysis = self.pcr_simulator.analyze_amplicon_specificity(
                amplicons, target_templates
            )
            
            # Predict efficiency
            efficiency_predictions = self.pcr_simulator.predict_pcr_efficiency(
                amplicons, 60.0, 60.0  # Default Tm values
            )
            
            # Generate report
            report = self.report_generator.generate_report(
                amplicons, specificity_analysis, efficiency_predictions
            )
            
            return {
                'amplicons': amplicons,
                'specificity_analysis': specificity_analysis,
                'efficiency_predictions': efficiency_predictions,
                'report': report
            }
            
        except Exception as e:
            self.logger.error(f"PCR simulation failed: {str(e)}")
            return {'error': str(e)}
    
    def _calculate_overall_score(self, blast_results: Dict, pcr_simulation: Dict) -> float:
        """Calculate overall specificity score"""
        scores = []
        
        # BLAST specificity score
        if 'overall_specificity' in blast_results:
            blast_score = blast_results['overall_specificity']['score']
            scores.append(blast_score * 0.6)  # 60% weight
        
        # PCR simulation specificity score
        if 'specificity_analysis' in pcr_simulation:
            pcr_score = pcr_simulation['specificity_analysis'].get('specificity_percentage', 0)
            scores.append(pcr_score * 0.4)  # 40% weight
        
        return sum(scores) if scores else 0.0
    
    def _generate_comprehensive_recommendations(self, blast_results: Dict, 
                                             pcr_simulation: Dict, 
                                             overall_score: float) -> List[str]:
        """Generate comprehensive recommendations"""
        recommendations = []
        
        # Overall score recommendations
        if overall_score >= 90:
            recommendations.append("Excellent primer specificity. Proceed with confidence.")
        elif overall_score >= 80:
            recommendations.append("Good primer specificity. Consider experimental validation.")
        elif overall_score >= 70:
            recommendations.append("Moderate specificity. Optimize PCR conditions or consider redesign.")
        else:
            recommendations.append("Poor specificity. Primer redesign strongly recommended.")
        
        # BLAST-specific recommendations
        if 'overall_specificity' in blast_results:
            blast_rec = blast_results['overall_specificity'].get('recommendation', '')
            if blast_rec and blast_rec not in recommendations:
                recommendations.append(f"BLAST analysis: {blast_rec}")
        
        # PCR simulation recommendations
        if 'report' in pcr_simulation and 'recommendations' in pcr_simulation['report']:
            pcr_recs = pcr_simulation['report']['recommendations']
            for rec in pcr_recs:
                if rec not in recommendations:
                    recommendations.append(f"PCR simulation: {rec}")
        
        # Error handling
        if 'error' in blast_results:
            recommendations.append(f"BLAST analysis failed: {blast_results['error']}")
        
        if 'error' in pcr_simulation:
            recommendations.append(f"PCR simulation failed: {pcr_simulation['error']}")
        
        return recommendations
    
    def export_results(self, results: List[SpecificityResult], 
                      output_format: str = 'json',
                      filename: str = 'specificity_results') -> str:
        """
        Export results to file
        
        Args:
            results: List of SpecificityResult objects
            output_format: 'json', 'csv', or 'html'
            filename: Output filename (without extension)
            
        Returns:
            str: Path to exported file
        """
        if output_format == 'json':
            return self._export_json(results, filename)
        elif output_format == 'csv':
            return self._export_csv(results, filename)
        elif output_format == 'html':
            return self._export_html(results, filename)
        else:
            raise ValueError(f"Unsupported output format: {output_format}")
    
    def _export_json(self, results: List[SpecificityResult], filename: str) -> str:
        """Export results to JSON"""
        output_file = f"{filename}.json"
        
        json_data = []
        for result in results:
            # Convert to JSON-serializable format
            data = {
                'primer_name': result.primer_pair.name,
                'forward_primer': result.primer_pair.forward,
                'reverse_primer': result.primer_pair.reverse,
                'overall_score': result.overall_score,
                'recommendations': result.recommendations,
                'processing_time': result.processing_time,
                'blast_summary': self._extract_blast_summary(result.blast_results),
                'pcr_summary': self._extract_pcr_summary(result.pcr_simulation)
            }
            json_data.append(data)
        
        with open(output_file, 'w') as f:
            json.dump(json_data, f, indent=2)
        
        return output_file
    
    def _export_csv(self, results: List[SpecificityResult], filename: str) -> str:
        """Export results to CSV"""
        import csv
        
        output_file = f"{filename}.csv"
        
        with open(output_file, 'w', newline='') as f:
            writer = csv.writer(f)
            
            # Header
            writer.writerow([
                'Primer Name', 'Forward Primer', 'Reverse Primer',
                'Overall Score', 'BLAST Specificity', 'PCR Specificity',
                'Processing Time', 'Main Recommendation'
            ])
            
            # Data rows
            for result in results:
                blast_score = self._extract_blast_score(result.blast_results)
                pcr_score = self._extract_pcr_score(result.pcr_simulation)
                main_rec = result.recommendations[0] if result.recommendations else "No recommendations"
                
                writer.writerow([
                    result.primer_pair.name,
                    result.primer_pair.forward,
                    result.primer_pair.reverse,
                    f"{result.overall_score:.1f}",
                    f"{blast_score:.1f}",
                    f"{pcr_score:.1f}",
                    f"{result.processing_time:.2f}s",
                    main_rec
                ])
        
        return output_file
    
    def _export_html(self, results: List[SpecificityResult], filename: str) -> str:
        """Export results to HTML report"""
        output_file = f"{filename}.html"
        
        html_content = """
        <!DOCTYPE html>
        <html>
        <head>
            <title>Primer Specificity Analysis Report</title>
            <style>
                body { font-family: Arial, sans-serif; margin: 20px; }
                table { border-collapse: collapse; width: 100%; }
                th, td { border: 1px solid #ddd; padding: 8px; text-align: left; }
                th { background-color: #f2f2f2; }
                .high-score { background-color: #d4edda; }
                .medium-score { background-color: #fff3cd; }
                .low-score { background-color: #f8d7da; }
            </style>
        </head>
        <body>
            <h1>Primer Specificity Analysis Report</h1>
            <table>
                <tr>
                    <th>Primer Name</th>
                    <th>Forward Primer</th>
                    <th>Reverse Primer</th>
                    <th>Overall Score</th>
                    <th>Assessment</th>
                    <th>Key Recommendations</th>
                </tr>
        """
        
        for result in results:
            score_class = ("high-score" if result.overall_score >= 80 else
                          "medium-score" if result.overall_score >= 60 else
                          "low-score")
            
            assessment = ("Excellent" if result.overall_score >= 90 else
                         "Good" if result.overall_score >= 80 else
                         "Fair" if result.overall_score >= 60 else
                         "Poor")
            
            recommendations_text = "; ".join(result.recommendations[:2])  # Top 2 recommendations
            
            html_content += f"""
                <tr class="{score_class}">
                    <td>{result.primer_pair.name}</td>
                    <td>{result.primer_pair.forward}</td>
                    <td>{result.primer_pair.reverse}</td>
                    <td>{result.overall_score:.1f}</td>
                    <td>{assessment}</td>
                    <td>{recommendations_text}</td>
                </tr>
            """
        
        html_content += """
            </table>
        </body>
        </html>
        """
        
        with open(output_file, 'w') as f:
            f.write(html_content)
        
        return output_file
    
    def _extract_blast_summary(self, blast_results: Dict) -> Dict:
        """Extract BLAST summary for export"""
        if 'overall_specificity' in blast_results:
            return {
                'score': blast_results['overall_specificity']['score'],
                'assessment': blast_results['overall_specificity']['assessment']
            }
        return {'error': blast_results.get('error', 'Unknown error')}
    
    def _extract_pcr_summary(self, pcr_simulation: Dict) -> Dict:
        """Extract PCR simulation summary for export"""
        if 'report' in pcr_simulation:
            return pcr_simulation['report']['summary']
        return {'error': pcr_simulation.get('error', 'Unknown error')}
    
    def _extract_blast_score(self, blast_results: Dict) -> float:
        """Extract BLAST specificity score"""
        return blast_results.get('overall_specificity', {}).get('score', 0.0)
    
    def _extract_pcr_score(self, pcr_simulation: Dict) -> float:
        """Extract PCR specificity score"""
        return pcr_simulation.get('specificity_analysis', {}).get('specificity_percentage', 0.0)
    
    def _get_human_templates(self) -> Dict[str, str]:
        """Get common human template sequences for simulation"""
        # In a real implementation, these would be loaded from a database
        return {
            "GAPDH": "ATGGGGAAGGTGAAGGTCGGAGTCAACGGATTTGGTCGTATTGGGCGCCTGGTCACCAGGGCTGCTTTTAACTCTGGTAAAGTGGATATTGTTGCCATCAATGACCCCTT",
            "ACTB": "ATGGATGATGATATCGCCGCGCTCGTCGTCGACAACGGCTCCGGCATGTGCAAGGCCGGCTTCGCGGGCGACGATGCCCCCCGGGCCGTCTTCCCCTCCATCGTGGGGCGCC"
        }
    
    def _get_mouse_templates(self) -> Dict[str, str]:
        """Get common mouse template sequences for simulation"""
        return {
            "Gapdh": "ATGGGGAAGGTGAAGGTCGGAGTCAACGGATTTGGTCGTATTGGGCGCCTGGTCACCAGGGCTGCTTTTAACTCTGGTAAAGTGGATATTGTTGCCATCAATGACCCCTT",
            "Actb": "ATGGATGATGATATCGCCGCGCTCGTCGTCGACAACGGCTCCGGCATGTGCAAGGCCGGCTTCGCGGGCGACGATGCCCCCCGGGCCGTCTTCCCCTCCATCGTGGGGCGCC"
        }
    
    def _get_ecoli_templates(self) -> Dict[str, str]:
        """Get common E. coli template sequences for simulation"""
        return {
            "16S_rRNA": "AGAGTTTGATCCTGGCTCAGATTGAACGCTGGCGGCAGGCCTAACACATGCAAGTCGAACGGTAACAGGAAGAAGCTTGCTCTTTGCTGACGAGTGGCGGACGGGTGAGT",
            "recA": "ATGAGTGTAGAGAAAGCGGTTGTAAATGAAGAGTTAGGTGTAGGTGGCGACGAAGGTGGCTATCGTGAAGCGATGGTTGAAGAAGATGCAACCGAAGAAGGCGAAATTT"
        }


# Example usage and testing
if __name__ == "__main__":
    # Configure logging
    logging.basicConfig(level=logging.INFO)
    
    # Create test primer pairs
    primer_pairs = [
        PrimerPair(
            forward="ATGGCTAGCTAGCTAGCTAG",
            reverse="CTAGCTAGCTAGCTAGCCAT",
            name="Test_Primer_1",
            target_region="GAPDH"
        ),
        PrimerPair(
            forward="GCTAGCTAGCTAGCTAGCTA",
            reverse="TAGCTAGCTAGCTAGCTAG",
            name="Test_Primer_2",
            target_region="ACTB"
        )
    ]
    
    # Initialize checker
    checker = ComprehensiveSpecificityChecker()
    
    # Run single specificity check
    print("Running single primer specificity check...")
    result = checker.check_specificity(primer_pairs[0], organism='human')
    print(f"Overall score: {result.overall_score:.1f}")
    print(f"Recommendations: {result.recommendations}")
    
    # Run batch specificity check
    print("\nRunning batch specificity check...")
    batch_results = checker.batch_check_specificity(primer_pairs, organism='human')
    
    # Export results
    print("\nExporting results...")
    json_file = checker.export_results(batch_results, 'json', 'test_results')
    csv_file = checker.export_results(batch_results, 'csv', 'test_results')
    html_file = checker.export_results(batch_results, 'html', 'test_results')
    
    print(f"Results exported to: {json_file}, {csv_file}, {html_file}")