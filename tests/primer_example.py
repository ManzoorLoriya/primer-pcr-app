from src.primer_designer import PrimerDesigner
from src.specificity_checker import ComprehensiveSpecificityChecker, PrimerPair
from src.insilico_pcr import InSilicoPCR

def main():
    # Initialize the primer designer
    designer = PrimerDesigner()
    
    # Example sequence (replace this with your target sequence)
    target_sequence = "CTGCAGACCCAGCACCAGGGAAATGATCCAGAAATTGCAACCTCAGCCCCCTGGCCATCTGCTGATGCCACCACCCCCAGGTCCCTAATGGGCCTGGTGGCAGAGTTTGGGAAGATGGGCTCAGGGCTATATAAAGTCCACAAGGACCTAAGAGCCCCCAGTGCTGCTGGGCCAGCTGTATTCTGAGGTGGTCAGCACACAGGTCTGTGTCCTCCGTGCTAGATTGGGGCTGAGAGGCTGGGGGCTCTGGGTTGGCTGGGACAGGACATGGGATTCTTCCTTGTATTGGGGGTTTTGGCTGTTACTCTGTCTCTCCATCAGGTCATCATCCTTTCATCATGGCTCTGTGGATGCATCTCCTCACCGTGCTGGCCCTGCTGGCCCTCTGGGGGCCCAACACTAATCAGGCCTTTGTCAGCCGGCATCTGTGCGGCTCCAACTTAGTGGAGACATTGTATTCAGTGTGTCAGGATGATGGCTTCTTCTATATACCCAAGGACCGTCGGGAGCTAGAGGACCCACAGGGTGAGCCCCTACCTGCCATCCCTGCTGTTTCCGTGCCAGTACCCCAGCTGGCAGGGCATAAGTAAGCAGGAAGCTAATTCCAAGGAGAGTCGATGGGTTTGTTGAAAAGGGAGGCGGCTCTCTTGGTCATTTCGTAAAGTGGTGGTGGCTTCCTATAGCTGCTTTTAAGGGTAAAGGGTAACAGCTGCACCCTTCAGCTGTGGCTTCTGAGCACAACTGGACTCTTCCCTCCACTTGCCTTCGAATGACTGCCCTGGCCTCATGGCAACAGTAGCTCCCTGGTACCAATTTTATTATGCAGATTGCATCTTGGTGTTGATAGCCTTAGGGTAGCCTGGGGGCCATTCATGGGGCGCCCCATCCCTCCTTCCTCCCTGCCTCTGGACAAATGCTCCATGGAGCTCCAAGCTCTGCCACGTGGGAGGTGTGGGTCTCCAGCGCTCTGTGTGCCCAGCATGGCAGCCTCTGTCACCTGGACCAGCTCCCTGGGAGATGCAGTGAGAGGGTGGTAGTGTGGGGCCAGTGCGCAGGCATTCTGCTGCTCCTGACAGCATCTGCCCCTGTCTCTCTCCCCACTGCTGCTGCTCCTGTATTCTGGCACCTCACCCTGCAGTGGAGCAGACAGAACTGGGCATGGGCCTGGGGGCAGGTGGACTACAGCCCTTGGCACTGGAGATGGCACTACAGAAGCGTGGCATTGTGGATCAGTGCTGTACTGGCACCTGCACACGCCACCAGCTGCAGAGCTACTGCAACTAGACACCTGCCTTGAACCTGGCCTCCCACTCTCCCCTGGCAACCAATAAACCCCTTGAATGAGCCCCATTGAATGGTCTGTGTGTCATGGAGGGGGAGGGGCTGACTCAAGGGGGCACATGCATGCCAGCCTATCATCCAGGTTCATTGCAAGACCCCCTCTCTATGCTCTGTGCACCTCTAACACACCC"
    
    # Design primers
    print("Designing primers...")
    result = designer.design_primers(target_sequence, num_return=3)
    
    if result['success']:
        print(f"\nFound {result['num_pairs_found']} primer pairs")
        
        # Get the best primer pair
        best_pair = result['primer_pairs'][0]
        
        # Create a PrimerPair object for specificity checking
        primer_pair = PrimerPair(
            forward=best_pair['left_primer']['sequence'],
            reverse=best_pair['right_primer']['sequence'],
            name="Test_Primers",
            target_region=target_sequence,
            expected_size=best_pair['product_size']
        )
        
        # Initialize specificity checker
        checker = ComprehensiveSpecificityChecker()
        
        # Check primer specificity
        print("\nChecking primer specificity...")
        specificity_result = checker.check_specificity(
            primer_pair,
            organism='human',  # Change this to your target organism
            blast_enabled=True,
            pcr_sim_enabled=True
        )
        
        # Print results
        print("\nPrimer Design Results:")
        print(f"Forward Primer: {primer_pair.forward}")
        print(f"Reverse Primer: {primer_pair.reverse}")
        print(f"Expected Product Size: {primer_pair.expected_size} bp")
        
        print("\nPrimer Analysis:")
        print(f"GC Content (Forward): {best_pair['left_primer']['analysis']['gc_content']}%")
        print(f"GC Content (Reverse): {best_pair['right_primer']['analysis']['gc_content']}%")
        print(f"Tm Difference: {best_pair['tm_difference']}°C")
        
        if best_pair['dimer_analysis']:
            print("\nDimer Analysis:")
            for dimer in best_pair['dimer_analysis']:
                print(f"- {dimer['type']}: {dimer['description']}")
        
        print("\nSpecificity Results:")
        print(f"Overall Specificity Score: {specificity_result.overall_score}")
        if specificity_result.recommendations:
            print("\nRecommendations:")
            for rec in specificity_result.recommendations:
                print(f"- {rec}")
    else:
        print(f"Primer design failed: {result.get('error', 'Unknown error')}")

if __name__ == "__main__":
    main() 