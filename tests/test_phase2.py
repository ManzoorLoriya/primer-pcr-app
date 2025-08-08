<<<<<<< HEAD
#!/usr/bin/env python3
"""
Quick test for Phase 2: BLAST specificity and in-silico PCR.
"""

import logging
from src.blast_integration import SpecificityAnalyzer
from src.insilico_pcr import InSilicoPCR, PCRSimulationReport

if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)

    # Example primers - redesigned for better specificity
    forward = "GCTAGCTAGCTAGCTAGCTA"  # More unique sequence
    reverse = "TAGCTAGCTAGCTAGCTAGC"  # More unique sequence
    organism = "Homo sapiens"  # or None

    # 1) Specificity analysis
    spec_analyzer = SpecificityAnalyzer()
    spec_results = spec_analyzer.analyze_primer_pair(
        forward_primer=forward,
        reverse_primer=reverse,
        target_organism=organism,
        database='nt'
    )
    print("=== SPECIFICITY ANALYSIS ===")
    print(f"Forward specificity score: {spec_results['forward_primer']['specificity_score']:.1f}")
    print(f"Reverse specificity score: {spec_results['reverse_primer']['specificity_score']:.1f}")
    print(f"Overall assessment: {spec_results['overall_specificity']['assessment']}")
    print(f"Recommendation: {spec_results['overall_specificity']['recommendation']}\n")

    # 2) In-silico PCR simulation
    # Define template sequences with more realistic content
    templates = {
        "target_gene": "A" * 1000 + "GCTAGCTAGCTAGCTAGCTA" + "G" * 500 + "TAGCTAGCTAGCTAGCTAGC" + "T" * 1000,
        "off_target_chromosome": "C" * 500 + "GCTAGCTAGCTAGCTAGCTA" + "C" * 600 + "TAGCTAGCTAGCTAGCTAGC" + "C" * 500,
        "unrelated_sequence": "TTTTGGGGCCCCAAAA"  # no binding sites
    }

    pcr_sim = InSilicoPCR(max_product_size=2000, min_product_size=50)
    amplicons = pcr_sim.simulate_pcr(
        templates=templates,
        forward_primer=forward,
        reverse_primer=reverse,
        max_mismatches=1
    )

    print("=== IN-SILICO PCR AMPLICONS ===")
    print(f"Total amplicons found: {len(amplicons)}")
    for amp in amplicons:
        print(f"  {amp.template_id}: size={amp.amplicon_size} bp, "
              f"positions Fwd({amp.forward_start}-{amp.forward_end}) "
              f"Rev({amp.reverse_start}-{amp.reverse_end})")

    # 3) Specificity analysis on amplicons
    spec_analysis = pcr_sim.analyze_amplicon_specificity(
        amplicons=amplicons,
        target_template_ids=["target_gene"]
    )

    # 4) Efficiency predictions
    # For simplicity, use primer Tm = 60 °C for both
    eff_preds = pcr_sim.predict_pcr_efficiency(
        amplicons=amplicons,
        primer_tm_forward=60.0,
        primer_tm_reverse=60.0
    )

    # 5) Generate a full report
    reporter = PCRSimulationReport()
    report = reporter.generate_report(
        amplicons=amplicons,
        specificity_analysis=spec_analysis,
        efficiency_predictions=eff_preds
    )

    print("\n=== PCR SIMULATION REPORT SUMMARY ===")
    print(f"Specificity: {report['summary']['specificity_score']:.1f}%")
    print(f"Non-specific products: {report['summary']['non_specific_products']}")
    print("Recommendations:")
    for rec in report['recommendations']:
        print(f"  - {rec}")
=======
#!/usr/bin/env python3
"""
Quick test for Phase 2: BLAST specificity and in-silico PCR.
"""

import logging
from src.blast_integration import SpecificityAnalyzer
from src.insilico_pcr import InSilicoPCR, PCRSimulationReport

if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)

    # Example primers - redesigned for better specificity
    forward = "GCTAGCTAGCTAGCTAGCTA"  # More unique sequence
    reverse = "TAGCTAGCTAGCTAGCTAGC"  # More unique sequence
    organism = "Homo sapiens"  # or None

    # 1) Specificity analysis
    spec_analyzer = SpecificityAnalyzer()
    spec_results = spec_analyzer.analyze_primer_pair(
        forward_primer=forward,
        reverse_primer=reverse,
        target_organism=organism,
        database='nt'
    )
    print("=== SPECIFICITY ANALYSIS ===")
    print(f"Forward specificity score: {spec_results['forward_primer']['specificity_score']:.1f}")
    print(f"Reverse specificity score: {spec_results['reverse_primer']['specificity_score']:.1f}")
    print(f"Overall assessment: {spec_results['overall_specificity']['assessment']}")
    print(f"Recommendation: {spec_results['overall_specificity']['recommendation']}\n")

    # 2) In-silico PCR simulation
    # Define template sequences with more realistic content
    templates = {
        "target_gene": "A" * 1000 + "GCTAGCTAGCTAGCTAGCTA" + "G" * 500 + "TAGCTAGCTAGCTAGCTAGC" + "T" * 1000,
        "off_target_chromosome": "C" * 500 + "GCTAGCTAGCTAGCTAGCTA" + "C" * 600 + "TAGCTAGCTAGCTAGCTAGC" + "C" * 500,
        "unrelated_sequence": "TTTTGGGGCCCCAAAA"  # no binding sites
    }

    pcr_sim = InSilicoPCR(max_product_size=2000, min_product_size=50)
    amplicons = pcr_sim.simulate_pcr(
        templates=templates,
        forward_primer=forward,
        reverse_primer=reverse,
        max_mismatches=1
    )

    print("=== IN-SILICO PCR AMPLICONS ===")
    print(f"Total amplicons found: {len(amplicons)}")
    for amp in amplicons:
        print(f"  {amp.template_id}: size={amp.amplicon_size} bp, "
              f"positions Fwd({amp.forward_start}-{amp.forward_end}) "
              f"Rev({amp.reverse_start}-{amp.reverse_end})")

    # 3) Specificity analysis on amplicons
    spec_analysis = pcr_sim.analyze_amplicon_specificity(
        amplicons=amplicons,
        target_template_ids=["target_gene"]
    )

    # 4) Efficiency predictions
    # For simplicity, use primer Tm = 60 °C for both
    eff_preds = pcr_sim.predict_pcr_efficiency(
        amplicons=amplicons,
        primer_tm_forward=60.0,
        primer_tm_reverse=60.0
    )

    # 5) Generate a full report
    reporter = PCRSimulationReport()
    report = reporter.generate_report(
        amplicons=amplicons,
        specificity_analysis=spec_analysis,
        efficiency_predictions=eff_preds
    )

    print("\n=== PCR SIMULATION REPORT SUMMARY ===")
    print(f"Specificity: {report['summary']['specificity_score']:.1f}%")
    print(f"Non-specific products: {report['summary']['non_specific_products']}")
    print("Recommendations:")
    for rec in report['recommendations']:
        print(f"  - {rec}")
>>>>>>> 75f1c6a80d6b193bf8a1fa04c4fb03060d16c47f
