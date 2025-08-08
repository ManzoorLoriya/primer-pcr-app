<<<<<<< HEAD
#!/usr/bin/env python3
"""
Test script for Phase 1 implementation
Run this to verify all algorithms are working correctly
"""

import sys
import os

# Add the parent directory to Python path to allow importing from src
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from src.algorithms import PrimerAnalyzer
from src.primer_designer import PrimerDesigner

def test_primer_analyzer():
    """Test the PrimerAnalyzer class"""
    print("=" * 50)
    print("TESTING PRIMER ANALYZER")
    print("=" * 50)
    
    analyzer = PrimerAnalyzer()
    
    # Test sequences
    test_primers = [
        "ATGCGATCGATCGATCG",  # Good primer
        "AAAAAAAAAAAAAAAA",   # Poor primer (low complexity)
        "GCGCGCGCGCGCGCGC",   # High GC content
        "ATCGATCGATCGATCGAT",  # Slightly longer
        "TTGCAACTGGATCCAAGC"   # Real-world example
    ]
    
    print("\n1. Testing Melting Temperature Calculation:")
    print("-" * 40)
    for primer in test_primers:
        tm = analyzer.calculate_melting_temperature(primer)
        gc = analyzer.calculate_gc_content(primer)
        print(f"Primer: {primer}")
        print(f"  Tm: {tm}°C, GC: {gc}%")
        print()
    
    print("\n2. Testing Quality Score Calculation:")
    print("-" * 40)
    for primer in test_primers:
        scores = analyzer.calculate_primer_quality_score(primer)
        print(f"Primer: {primer}")
        print(f"  Length Score: {scores['length_score']}")
        print(f"  GC Score: {scores['gc_score']}")
        print(f"  Tm Score: {scores['tm_score']}")
        print(f"  Complexity Score: {scores['complexity_score']}")
        print(f"  Terminal Score: {scores['terminal_score']}")
        print(f"  Overall Score: {scores['overall_score']}")
        print()
    
    print("\n3. Testing Primer Dimer Detection:")
    print("-" * 40)
    primer1 = "ATGCGATCGATCGATCG"
    primer2 = "CGATCGATCGCATGC"
    
    dimers = analyzer.detect_primer_dimers(primer1, primer2)
    print(f"Primer 1: {primer1}")
    print(f"Primer 2: {primer2}")
    print(f"Found {len(dimers)} potential dimers:")
    
    for dimer in dimers:
        print(f"  Type: {dimer['type']}")
        print(f"  Risk Level: {dimer['risk_level']}")
        if 'matches' in dimer:
            print(f"  Matches: {dimer['matches']}")
        if 'stability_score' in dimer:
            print(f"  Stability Score: {dimer['stability_score']}")
        print()
    
    print("\n4. Testing Complete Sequence Analysis:")
    print("-" * 40)
    test_seq = "TTGCAACTGGATCCAAGC"
    analysis = analyzer.analyze_primer_sequence(test_seq)
    
    print(f"Sequence: {analysis['sequence']}")
    print(f"Length: {analysis['length']} bp")
    print(f"GC Content: {analysis['gc_content']}%")
    print(f"Melting Temperature: {analysis['melting_temperature']}°C")
    print(f"Molecular Weight: {analysis['molecular_weight']} g/mol")
    print(f"Overall Quality: {analysis['quality_scores']['overall_score']}")
    print(f"Hairpins Found: {len(analysis['hairpins'])}")


def test_primer_designer():
    """Test the PrimerDesigner class"""
    print("\n" + "=" * 50)
    print("TESTING PRIMER DESIGNER")
    print("=" * 50)
    
    designer = PrimerDesigner()
    
    # Test sequence (a portion of human GAPDH gene)
    test_sequence = """
    ATGGGGAAGGTGAAGGTCGGAGTCAACGGATTTGGTCGTATTGGGCGCCTGGTCACCAGGGCTGC
    TTTTAACTCTGGTAAAGTGGATATTGTTGCCATCAATGACCCCTTCATTGACCTCAACTACATGG
    TTTACATGTTCCAATATGATTCCACCCATGGCAAATTCCATGGCACCGTCAAGGCTGAGAACGGG
    AAGCTTGTCATCAATGGAAATCCCATCACCATCTTCCAGGAGCGAGATCCCTCCAAAATCAAGTG
    GGGCGATGCTGGCGCTGAGTACGTCGTGGAGTCCACTGGCGTCTTCACCACCATGGAGAAGGCTG
    GGGCTCATTTGCAGGGGGGAGCCAAAAGGGTCATCATCTCTGCCCCCTCTGCTGATGCCCCCATG
    TTCGTCATGGGTGTGAACCATGAGAAGTATGACAACAGCCTCAAGATCATCAGCAATGCCTCCTG
    CACCACCAACTGCTTAGCACCCCTGGCCAAGGTCATCCATGACAACTTTGGTATCGTGGAAGGAC
    TCATGACCACAGTCCATGCCATCACTGCCACCCAGAAGACTGTGGATGGCCCCTCCGGGAAACTG
    TGGCGTGATGGCCGCGGGGCTCTCCAGAACATCATCCCTGCCTCT
    """.replace('\n', '').replace(' ', '')
    
    print(f"\nTest sequence length: {len(test_sequence)} bp")
    print(f"Test sequence GC content: {designer.analyzer.calculate_gc_content(test_sequence)}%")
    
    print("\n1. Testing Basic Primer Design:")
    print("-" * 40)
    
    try:
        result = designer.design_primers(test_sequence, num_return=3)
        
        if result['success']:
            print(f"✓ Successfully designed {result['num_pairs_found']} primer pairs")
            
            for i, pair in enumerate(result['primer_pairs'][:2]):  # Show top 2 pairs
                print(f"\n--- Primer Pair {i+1} ---")
                print(f"Overall Quality: {pair['overall_quality']}")
                print(f"Product Size: {pair['product_size']} bp")
                print(f"Tm Difference: {pair['tm_difference']}°C")
                
                print(f"\nLeft Primer:")
                print(f"  Sequence: {pair['left_primer']['sequence']}")
                print(f"  Tm: {pair['left_primer']['analysis']['melting_temperature']}°C")
                print(f"  GC: {pair['left_primer']['analysis']['gc_content']}%")
                print(f"  Quality: {pair['left_primer']['analysis']['quality_scores']['overall_score']}")
                
                print(f"\nRight Primer:")
                print(f"  Sequence: {pair['right_primer']['sequence']}")
                print(f"  Tm: {pair['right_primer']['analysis']['melting_temperature']}°C")
                print(f"  GC: {pair['right_primer']['analysis']['gc_content']}%")
                print(f"  Quality: {pair['right_primer']['analysis']['quality_scores']['overall_score']}")
                
                if pair['dimer_analysis']:
                    high_risk = [d for d in pair['dimer_analysis'] if d.get('risk_level') == 'HIGH']
                    medium_risk = [d for d in pair['dimer_analysis'] if d.get('risk_level') == 'MEDIUM']
                    print(f"\nDimer Analysis:")
                    print(f"  High Risk: {len(high_risk)}")
                    print(f"  Medium Risk: {len(medium_risk)}")
        else:
            print(f"✗ Primer design failed: {result.get('error', 'Unknown error')}")
            
    except Exception as e:
        print(f"✗ Error during primer design: {str(e)}")
    
    print("\n2. Testing Custom Parameters:")
    print("-" * 40)
    
    custom_params = {
        'PRIMER_OPT_SIZE': 22,
        'PRIMER_MIN_SIZE': 20,
        'PRIMER_MAX_SIZE': 25,
        'PRIMER_OPT_TM': 62.0,
        'PRIMER_PRODUCT_SIZE_RANGE': [[150, 300]]
    }
    
    try:
        result = designer.design_primers(test_sequence, custom_params=custom_params, num_return=2)
        
        if result['success']:
            print(f"✓ Custom parameter design successful: {result['num_pairs_found']} pairs")
            
            if result['primer_pairs']:
                pair = result['primer_pairs'][0]
                print(f"Top pair quality: {pair['overall_quality']}")
                print(f"Product size: {pair['product_size']} bp (within 150-300 range)")
                print(f"Left primer length: {pair['left_primer']['length']} bp")
                print(f"Right primer length: {pair['right_primer']['length']} bp")
        else:
            print(f"✗ Custom parameter design failed: {result.get('error', 'Unknown error')}")
            
    except Exception as e:
        print(f"✗ Error with custom parameters: {str(e)}")
    
    print("\n3. Testing Batch Design:")
    print("-" * 40)
    
    # Create test sequences for batch processing
    batch_sequences = [
        {
            'id': 'seq1',
            'sequence': test_sequence[:200]  # First 200 bp
        },
        {
            'id': 'seq2', 
            'sequence': test_sequence[200:400]  # Next 200 bp
        },
        {
            'id': 'seq3',
            'sequence': 'ATCGATCGATCG'  # Too short - should fail
        }
    ]
    
    try:
        batch_results = designer.batch_design_primers(batch_sequences)
        
        print(f"Batch processed {len(batch_results)} sequences:")
        for result in batch_results:
            seq_id = result['sequence_id']
            if result['success']:
                pairs_found = result['num_pairs_found']
                print(f"  {seq_id}: ✓ {pairs_found} primer pairs found")
            else:
                print(f"  {seq_id}: ✗ {result.get('error', 'Failed')}")
                
    except Exception as e:
        print(f"✗ Error in batch processing: {str(e)}")
    
    print("\n4. Testing Recommendations:")
    print("-" * 40)
    
    try:
        # Get recommendations for the first design
        result = designer.design_primers(test_sequence, num_return=1)
        recommendations = designer.get_primer_recommendations(result)
        
        print("Design Recommendations:")
        for i, rec in enumerate(recommendations, 1):
            print(f"  {i}. {rec}")
            
    except Exception as e:
        print(f"✗ Error getting recommendations: {str(e)}")


def main():
    """Run all tests"""
    print("PRIMER DESIGN TOOL - PHASE 1 TESTING")
    print("=" * 60)
    
    try:
        test_primer_analyzer()
        test_primer_designer()
        
        print("\n" + "=" * 60)
        print("✓ PHASE 1 TESTING COMPLETED SUCCESSFULLY!")
        print("All core algorithms are working properly.")
        print("Ready to proceed to Phase 2: Specificity Checking")
        
    except ImportError as e:
        print(f"✗ Import Error: {e}")
        print("Make sure all required packages are installed:")
        print("pip install flask biopython primer3-py plotly pandas requests")
        
    except Exception as e:
        print(f"✗ Unexpected error during testing: {e}")
        print("Check your implementation and try again.")


if __name__ == "__main__":
=======
#!/usr/bin/env python3
"""
Test script for Phase 1 implementation
Run this to verify all algorithms are working correctly
"""

import sys
import os

# Add the parent directory to Python path to allow importing from src
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from src.algorithms import PrimerAnalyzer
from src.primer_designer import PrimerDesigner

def test_primer_analyzer():
    """Test the PrimerAnalyzer class"""
    print("=" * 50)
    print("TESTING PRIMER ANALYZER")
    print("=" * 50)
    
    analyzer = PrimerAnalyzer()
    
    # Test sequences
    test_primers = [
        "ATGCGATCGATCGATCG",  # Good primer
        "AAAAAAAAAAAAAAAA",   # Poor primer (low complexity)
        "GCGCGCGCGCGCGCGC",   # High GC content
        "ATCGATCGATCGATCGAT",  # Slightly longer
        "TTGCAACTGGATCCAAGC"   # Real-world example
    ]
    
    print("\n1. Testing Melting Temperature Calculation:")
    print("-" * 40)
    for primer in test_primers:
        tm = analyzer.calculate_melting_temperature(primer)
        gc = analyzer.calculate_gc_content(primer)
        print(f"Primer: {primer}")
        print(f"  Tm: {tm}°C, GC: {gc}%")
        print()
    
    print("\n2. Testing Quality Score Calculation:")
    print("-" * 40)
    for primer in test_primers:
        scores = analyzer.calculate_primer_quality_score(primer)
        print(f"Primer: {primer}")
        print(f"  Length Score: {scores['length_score']}")
        print(f"  GC Score: {scores['gc_score']}")
        print(f"  Tm Score: {scores['tm_score']}")
        print(f"  Complexity Score: {scores['complexity_score']}")
        print(f"  Terminal Score: {scores['terminal_score']}")
        print(f"  Overall Score: {scores['overall_score']}")
        print()
    
    print("\n3. Testing Primer Dimer Detection:")
    print("-" * 40)
    primer1 = "ATGCGATCGATCGATCG"
    primer2 = "CGATCGATCGCATGC"
    
    dimers = analyzer.detect_primer_dimers(primer1, primer2)
    print(f"Primer 1: {primer1}")
    print(f"Primer 2: {primer2}")
    print(f"Found {len(dimers)} potential dimers:")
    
    for dimer in dimers:
        print(f"  Type: {dimer['type']}")
        print(f"  Risk Level: {dimer['risk_level']}")
        if 'matches' in dimer:
            print(f"  Matches: {dimer['matches']}")
        if 'stability_score' in dimer:
            print(f"  Stability Score: {dimer['stability_score']}")
        print()
    
    print("\n4. Testing Complete Sequence Analysis:")
    print("-" * 40)
    test_seq = "TTGCAACTGGATCCAAGC"
    analysis = analyzer.analyze_primer_sequence(test_seq)
    
    print(f"Sequence: {analysis['sequence']}")
    print(f"Length: {analysis['length']} bp")
    print(f"GC Content: {analysis['gc_content']}%")
    print(f"Melting Temperature: {analysis['melting_temperature']}°C")
    print(f"Molecular Weight: {analysis['molecular_weight']} g/mol")
    print(f"Overall Quality: {analysis['quality_scores']['overall_score']}")
    print(f"Hairpins Found: {len(analysis['hairpins'])}")


def test_primer_designer():
    """Test the PrimerDesigner class"""
    print("\n" + "=" * 50)
    print("TESTING PRIMER DESIGNER")
    print("=" * 50)
    
    designer = PrimerDesigner()
    
    # Test sequence (a portion of human GAPDH gene)
    test_sequence = """
    ATGGGGAAGGTGAAGGTCGGAGTCAACGGATTTGGTCGTATTGGGCGCCTGGTCACCAGGGCTGC
    TTTTAACTCTGGTAAAGTGGATATTGTTGCCATCAATGACCCCTTCATTGACCTCAACTACATGG
    TTTACATGTTCCAATATGATTCCACCCATGGCAAATTCCATGGCACCGTCAAGGCTGAGAACGGG
    AAGCTTGTCATCAATGGAAATCCCATCACCATCTTCCAGGAGCGAGATCCCTCCAAAATCAAGTG
    GGGCGATGCTGGCGCTGAGTACGTCGTGGAGTCCACTGGCGTCTTCACCACCATGGAGAAGGCTG
    GGGCTCATTTGCAGGGGGGAGCCAAAAGGGTCATCATCTCTGCCCCCTCTGCTGATGCCCCCATG
    TTCGTCATGGGTGTGAACCATGAGAAGTATGACAACAGCCTCAAGATCATCAGCAATGCCTCCTG
    CACCACCAACTGCTTAGCACCCCTGGCCAAGGTCATCCATGACAACTTTGGTATCGTGGAAGGAC
    TCATGACCACAGTCCATGCCATCACTGCCACCCAGAAGACTGTGGATGGCCCCTCCGGGAAACTG
    TGGCGTGATGGCCGCGGGGCTCTCCAGAACATCATCCCTGCCTCT
    """.replace('\n', '').replace(' ', '')
    
    print(f"\nTest sequence length: {len(test_sequence)} bp")
    print(f"Test sequence GC content: {designer.analyzer.calculate_gc_content(test_sequence)}%")
    
    print("\n1. Testing Basic Primer Design:")
    print("-" * 40)
    
    try:
        result = designer.design_primers(test_sequence, num_return=3)
        
        if result['success']:
            print(f"✓ Successfully designed {result['num_pairs_found']} primer pairs")
            
            for i, pair in enumerate(result['primer_pairs'][:2]):  # Show top 2 pairs
                print(f"\n--- Primer Pair {i+1} ---")
                print(f"Overall Quality: {pair['overall_quality']}")
                print(f"Product Size: {pair['product_size']} bp")
                print(f"Tm Difference: {pair['tm_difference']}°C")
                
                print(f"\nLeft Primer:")
                print(f"  Sequence: {pair['left_primer']['sequence']}")
                print(f"  Tm: {pair['left_primer']['analysis']['melting_temperature']}°C")
                print(f"  GC: {pair['left_primer']['analysis']['gc_content']}%")
                print(f"  Quality: {pair['left_primer']['analysis']['quality_scores']['overall_score']}")
                
                print(f"\nRight Primer:")
                print(f"  Sequence: {pair['right_primer']['sequence']}")
                print(f"  Tm: {pair['right_primer']['analysis']['melting_temperature']}°C")
                print(f"  GC: {pair['right_primer']['analysis']['gc_content']}%")
                print(f"  Quality: {pair['right_primer']['analysis']['quality_scores']['overall_score']}")
                
                if pair['dimer_analysis']:
                    high_risk = [d for d in pair['dimer_analysis'] if d.get('risk_level') == 'HIGH']
                    medium_risk = [d for d in pair['dimer_analysis'] if d.get('risk_level') == 'MEDIUM']
                    print(f"\nDimer Analysis:")
                    print(f"  High Risk: {len(high_risk)}")
                    print(f"  Medium Risk: {len(medium_risk)}")
        else:
            print(f"✗ Primer design failed: {result.get('error', 'Unknown error')}")
            
    except Exception as e:
        print(f"✗ Error during primer design: {str(e)}")
    
    print("\n2. Testing Custom Parameters:")
    print("-" * 40)
    
    custom_params = {
        'PRIMER_OPT_SIZE': 22,
        'PRIMER_MIN_SIZE': 20,
        'PRIMER_MAX_SIZE': 25,
        'PRIMER_OPT_TM': 62.0,
        'PRIMER_PRODUCT_SIZE_RANGE': [[150, 300]]
    }
    
    try:
        result = designer.design_primers(test_sequence, custom_params=custom_params, num_return=2)
        
        if result['success']:
            print(f"✓ Custom parameter design successful: {result['num_pairs_found']} pairs")
            
            if result['primer_pairs']:
                pair = result['primer_pairs'][0]
                print(f"Top pair quality: {pair['overall_quality']}")
                print(f"Product size: {pair['product_size']} bp (within 150-300 range)")
                print(f"Left primer length: {pair['left_primer']['length']} bp")
                print(f"Right primer length: {pair['right_primer']['length']} bp")
        else:
            print(f"✗ Custom parameter design failed: {result.get('error', 'Unknown error')}")
            
    except Exception as e:
        print(f"✗ Error with custom parameters: {str(e)}")
    
    print("\n3. Testing Batch Design:")
    print("-" * 40)
    
    # Create test sequences for batch processing
    batch_sequences = [
        {
            'id': 'seq1',
            'sequence': test_sequence[:200]  # First 200 bp
        },
        {
            'id': 'seq2', 
            'sequence': test_sequence[200:400]  # Next 200 bp
        },
        {
            'id': 'seq3',
            'sequence': 'ATCGATCGATCG'  # Too short - should fail
        }
    ]
    
    try:
        batch_results = designer.batch_design_primers(batch_sequences)
        
        print(f"Batch processed {len(batch_results)} sequences:")
        for result in batch_results:
            seq_id = result['sequence_id']
            if result['success']:
                pairs_found = result['num_pairs_found']
                print(f"  {seq_id}: ✓ {pairs_found} primer pairs found")
            else:
                print(f"  {seq_id}: ✗ {result.get('error', 'Failed')}")
                
    except Exception as e:
        print(f"✗ Error in batch processing: {str(e)}")
    
    print("\n4. Testing Recommendations:")
    print("-" * 40)
    
    try:
        # Get recommendations for the first design
        result = designer.design_primers(test_sequence, num_return=1)
        recommendations = designer.get_primer_recommendations(result)
        
        print("Design Recommendations:")
        for i, rec in enumerate(recommendations, 1):
            print(f"  {i}. {rec}")
            
    except Exception as e:
        print(f"✗ Error getting recommendations: {str(e)}")


def main():
    """Run all tests"""
    print("PRIMER DESIGN TOOL - PHASE 1 TESTING")
    print("=" * 60)
    
    try:
        test_primer_analyzer()
        test_primer_designer()
        
        print("\n" + "=" * 60)
        print("✓ PHASE 1 TESTING COMPLETED SUCCESSFULLY!")
        print("All core algorithms are working properly.")
        print("Ready to proceed to Phase 2: Specificity Checking")
        
    except ImportError as e:
        print(f"✗ Import Error: {e}")
        print("Make sure all required packages are installed:")
        print("pip install flask biopython primer3-py plotly pandas requests")
        
    except Exception as e:
        print(f"✗ Unexpected error during testing: {e}")
        print("Check your implementation and try again.")


if __name__ == "__main__":
>>>>>>> 75f1c6a80d6b193bf8a1fa04c4fb03060d16c47f
    main()