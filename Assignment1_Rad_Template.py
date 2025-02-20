import argparse
import pyfaidx
import pdb
import vcf
from Bio.Restriction.Restriction_Dictionary import rest_dict


"""
Setup: Create a conda environment for this homework and install pyfaidx, vcf, and biopython. Confirm you can import 
each of these modules into python

General Instructions: Below you'll find a code outline for the assignment. Your job is to fill in the missing code
wherever there is a comment starting with "TODO". All of the missing code is contained within functions, which helps
keep individual tasks manageable. Each function already has type-hinted arguments and return variables defined (look up
"type hints in python" if this is confusing), as well as a doc-string with additional instructions. Often, you will
encounter instances where a certain variable is set to "None". These are place-holders, and you will usually want to
replace the "None" with a meaningful definition. 

At the bottom of this script are two sections of code, both commented
out: one called "testing code", the other called "main code". The testing code is meant to be uncommented incrementally,
and is designed to help you check your work as you go. Once you finish filling in all of the missing code segments,
you can uncomment the "main code" section to bring everything together into a working script. 

"""

def my_parse_args() -> argparse.Namespace:
    """
    Task 1:
    Implement a function that uses argparse to parse the command-line options. This function has no arguments itself,
    but should return an argparse.Namespace object named "parsed_args" containing the command-line arguments.
    Set up your parser to take in the following positional arguments (in this order):
        GenomeFile: the path to the .fasta file containing the genome
        VCFFile: the path to the .vcf file with the DNA polymorphisms
        Mode: whether to run SingleRad or ddRad
        RE1: the name (not the motif) of the restriction enzyme
    As well as the following flagged argument:
        RE2: the name (not the motif) of the second restriction enzyme (only used in ddRad mode)

    :return parsed_args: the parsed command-line arguments
    """

    parser = argparse.ArgumentParser(
        description='this script performs either SingleRad or ddRad sequencing on DNA from an input FASTA file and compared'
                    'the results to variants from an input VCF file')
    parser.add_argument("GenomeFile", type=str, help="The FASTA file to analyze")
    parser.add_argument("VCFFile", type=str, help="The FASTA file to analyze")
    parser.add_argument("Mode", type=str, help="The type of RadSeq you want to run. Options: SingleRad, ddRad")
    parser.add_argument("RE1", type=str, help="The name of the restriction enzyme used for RadSeq, and the first enzyme for ddRadseq")
    parser.add_argument("--RE2", type=str, help="The name of the second restriction site for ddRad")

    parsed_args = parser.parse_args()
    return parsed_args


def read_fasta(path_to_fasta: str, chromosome: str) -> str:
    """
    Task 2:
    Implement a function to read the fasta file using pyfaidx. This function should return the dna sequence (as a str)
    associated with the specified chromosome.

    :param path_to_fasta: path to the fasta file, as a str
    :param chromosome: which chromosome to read. The default value of 'NC_036780.1' refers to the first chromosome
    :return dna: the dna associated with the specified chromosome, as a str
    """
    genome = pyfaidx.Fasta(path_to_fasta)
    for nucleotides in genome:
        dna = nucleotides[:].seq

    return dna


def find_motifs(dna: str, motif: str) -> list[int]:
    """
    Task 3:
    Implement the find_motifs function to take DNA sequence (string, as returned by your read_fasta function)
    and a Restriction Enzyme motif (string) and will return a list of locations (integers) where the motif is found

    :param dna: DNA sequence to be analyzed
    :param motif: restriction enzyme motif
    :return positions: list of ints of the motif positions
    """
    
    positions = []
    start = 0
    
    while True:
        start = dna.find(motif, start)
        if start == -1:
            break
        positions.append(start)
        start += len(motif)

    return positions


def run_single_rad(dna: str, re1: str) -> list[tuple[int, int]]:
    """
    Task 4:
    implement a function for performing SingleRad. This function takes in the dna (str) and the name of the restriction
    enzyme (str), and returns a list of tuples where each tuple contains the start index (int) and stop index (int) of a
    site sequenced by SingleRad

    :param dna: DNA sequence to be analyzed
    :param re1: restriction enzyme name
    :return sequenced_sites: a list of tuples, where tuple contains the start and stop index of a sequenced site
    """
    assert re1 in rest_dict, f"Restriction enzyme {re1} not found in rest_dict."

    seq_length = 100
    motif = rest_dict[re1]['site']
    re1_cuts = find_motifs(dna, motif)
    
    sequenced_sites = [
        (max(0, site - seq_length), min(len(dna), site + seq_length))
        for site in re1_cuts
    ]
    

    #TODO: use an assert statement to verify that the RE1 site is in the rest_dict provided by biopython Good
    #TODO: use the rest_dict to look up the motif associated with the enzyme Good
    #TODO: use the find_motifs function to find the cut sites  Good
    #TODO: use list comprehension (preferably) or a loop (acceptable) to create the desired list of sequence sites Good

    return sequenced_sites


def run_ddrad(dna: str, re1: str, re2: str) -> list[tuple[int, int]]:
    """
    Task 5:
    implement a function for performing ddRad. This function takes in the dna (str) and the names of the restriction
    enzymes (str and str), and returns a list of tuples where each tuple contains the start index (int) and stop index
    (int) of a site sequenced by ddRad
    :param dna: DNA sequence to be analyzed
    :param re1: first restriction enzyme name
    :param re2: second restriction enzyme name
    :return sequenced_sites: a list of tuples, where tuple contains the start and stop index of a sequenced site
    """
    min_size = 300 # This variable stores the shortest distance for the DNA piece
    max_size = 700 # This variable stores the longest distance for the DNA piece
    seq_length = 100 # This variable stores how long the sequencing reads are

    assert re1 in rest_dict, f"Restriction enzyme {re1} not found in rest_dict."
    assert re2 in rest_dict, f"Restriction enzyme {re2} not found in rest_dict."
    motif1 = rest_dict[re1]['site']
    motif2 = rest_dict[re2]['site']

    #TODO: Use an assert statement to verify that RE1 and RE2 are both in the rest_dict provided by biopython GOOD
    #TODO: use the rest_dict to look up the motifs associated with the two enzyme  GOOD
    #TODO: use your find_motifs function to find the sequencing sites for RE1 and RE2 GOOD
    re1_sites=find_motifs(dna, motif1)
    re2_sites=find_motifs(dna, motif2)

    #TODO: Complete the below loop, which is currently incomplete but has hints as comments
    sequenced_sites = []
    for r1 in re1_sites: 
        
        # TODO: Look left of the current site
        left_hits_r1 = [r for r in re1_sites if r < r1] # Replace None with a list comprehension to find all re1 sites left of the current site
        left_hits_r2 = [r2 for r2 in re2_sites if r2 < r1] # Replace None with a list comprehension to find all re2 sites left of the current site

        if left_hits_r2:
            nearest_r2 = max(left_hits_r2)
            fragment_size = r1 - nearest_r2

            if min_size < fragment_size < max_size:
                if not left_hits_r1 or max(left_hits_r1) < nearest_r2:
                    sequenced_sites.append((nearest_r2, nearest_r2 + seq_length))
                    sequenced_sites.append((r1 - seq_length, r1))


        #TODO: create a series of nested if statements that check:
        # 1) that left_hits_r2 is not empty good
        # 2) that re2 site is between min and max size good
        # 3) that re2 site is closer re1 site good
        # TODO: if all statements pass, add the two resulting sequence sites to the sequenced_sites list good

        # TODO: Now do the same thing but looking right
        right_hits_r1 = [r for r in re1_sites if r > r1]
        right_hits_r2 = [r2 for r2 in re2_sites if r2 > r1]

        if right_hits_r2:
            nearest_r2 = min(right_hits_r2)
            fragment_size = nearest_r2 - r1

            if min_size < fragment_size < max_size:
                if not right_hits_r1 or min(right_hits_r1) > nearest_r2:
                    sequenced_sites.append((r1, r1 + seq_length))
                    sequenced_sites.append((nearest_r2 - seq_length, nearest_r2))

    return sequenced_sites

def find_variable_sites(vcf_file_path: str, sequenced_sites: list[tuple[int, int]]):
    """
    Task 6:
    Implement this function to take in a path to a VCF file and a list of sequenced sites (as returned by run_single_rad
    or run_ddrad) and return a list of sites that are both sequenced and contain variation.

    :param vcf_file_path: path to the VCF file
    :param sequenced_sites: list of sequenced sites, as returned by run_single_rad or run_ddrad
    :return sequenced_sites_variable: a list of sequenced sites that contain variation. This list is a subset of the
    sequenced_sites list you input to the function
    """

    #TODO: read in the VCF file using the VCFReader class
    with open (vcf_file_path, "r") as vcf_file:
        vcf_obj = vcf.Reader(vcf_file)

        num_called_filtered = 0
        
        #TODO: complete the below loop
        sequenced_sites_variable = [] # This list will keep track of all of the sequenced DNA pieces that also contain variation
        for record in vcf_obj: # Create loop to read through every record
            if record.CHROM != 'NC_036780.1': # These two lines will make sure that you are only looking at records from Chromosome 1
                break
            if record.num_called ==2:
                sample1GT = record.genotype(record.samples[0].sample).gt_nums
                sample2GT = record.genotype(record.samples[1].sample).gt_nums    
                if (sample1GT == "0/0" and sample2GT == "1/1") or (sample1GT == "1/1" and sample2GT == "0/0"):
                    variants_POS = [x for x in sequenced_sites if x[0] < record.POS < x[1]]
                    sequenced_sites_variable.extend(variants_POS)
        
        

        
        # TODO:add an if statement here to ensure that both samples have called genotypes (num_called attribute is your friend) good
        # TODO:add another if statement to ensure that one genotype is 0/0 and the other is 1/1 (use gt_nums)
        # TODO:if both if conditions pass, Use list comprehension to determine if this variants falls within sequenced DNA
        # TODO: if the variant falls within sequenced DNA, add it to the sequenced_sites_variables list
    return sequenced_sites_variable


# ###TESTING CODE###
# # This section of code is provided to help you test your functions as you develop them. Feel free to add additional
# # debugging code to this section, but ensure that all code required to run the script remains in the function
# # definitions above. As you progress, leave the completed sections uncommented, as some variables are reused. Note
# # also that, in much of this testing code, your command-line input (other than file locations) is ignored, but you
# # will likely still need to provide a value for each argument to get your code to run.
#
# Task 1: Parsing Command-line Arguments. Uncomment this section when you have completed the my_parse_args function
# and ensure the args that print match your expectations. This section also checks that the RE1 argument you provide
# # is a valid restriction enzyme name from the rest_dict
# print('testing my_parse_args function')
# args = my_parse_args()
# print(f'the following arguments were received: {args}')
# if args.RE1 not in rest_dict:
#     print(f'no restriction enzyme named {args.RE1} found in rest_dict. Run rest_dict.keys() to see valid names')
# print('\n')
# #
# #
# # # Task2: Reading the Fasta File. Uncomment this section when you have completed the read_fasta function.
# print('testing read_fasta function')
# target_chromosome = 'NC_036780.1'
# dna_string = read_fasta(path_to_fasta=args.GenomeFile, chromosome=target_chromosome)
# print(f'fasta file read. Length of {target_chromosome} is {len(dna_string)}')
# print(f'length matches expected length of 45605991: {len(dna_string)==45605991}')
# print('\n')
# #
# # # Task 3: Finding Motifs. Uncomment this section when you have completed the find_motifs function.
# print(f'testing find_motifs function')
# motif_locations = find_motifs('GGTTAAAGATCGGCGAGCCAATGGATCGACGATCA','GATC')
# print(f'find_motifs returned {motif_locations}')
# print(f'motif locations match expected locations of [7,23,30]: {motif_locations==[7,23,30]}')
# print('\n')
# #
# # Task 4: Implementing SingleRad. Uncomment this section when you have completed the run_single_rad function. Note that
# # this section hard-codes the RE1 argument for consistency, rather than using args.RE1 to get a value from the command
# # line
# print('testing run_single_rad function')
# srad_seq_sites = run_single_rad(dna=dna_string, re1='AanI')
# print(f'located {len(srad_seq_sites)} singlerad sites for re1=AanI')
# print(f'number of sites matches expected number (19280 sites): {len(srad_seq_sites) == 19280}')
# print(f'first sequence site determined to be {srad_seq_sites[0]}')
# print(f'first sequence site matches expected site (44450, 44650): {srad_seq_sites[0] == (44450, 44650)}')
# print('\n')
# #
# #
# # Task 5: Implementing ddRad. Uncomment this section when you have completed the run_ddrad function. Note again that
# # the re1 and re2 arguments are hard-coded here for consistency. This section may  take longer to run than others
# # due to the complexity of ddrad, but should still run in well under a minute. If it doesn't, check your code
# # for infinite loops
# print('testing run_ddrad function')
# ddrad_seq_sites = run_ddrad(dna_string, re1='AanI', re2='MroI')
# print(f'located {len(ddrad_seq_sites)} ddRad sites for re1=AanI, re2=MroI')
# print(f'number of sites matches expected number (848 sites): {len(ddrad_seq_sites) == 848}')
# print(f'first sequence site determined to be {ddrad_seq_sites[0]}')
# print(f'first sequence site matches expected site (91143, 91243): {ddrad_seq_sites[0] == (91143, 91243)}')
# print('\n')

# # # Task 6: Finding Sites with Variation. Uncomment this section when you have completed the find_variable_sites function.
# # # Again, this section will be slightly slower to run, but should still run in less than a minute
# print('testing find_variable_sites function')
# variable_srad_seq_sites = find_variable_sites(vcf_file_path=args.VCFFile, sequenced_sites=srad_seq_sites)
# print(f'found {len(variable_srad_seq_sites)} singleRad sites with variation')
# print(f'number of sites found matches expectation (2483 sites): {len(variable_srad_seq_sites)==2483}')
# print(f'{len(100*variable_srad_seq_sites)/len(srad_seq_sites):.3f}% of sites sequenced with single rad had variation')
# variable_ddrad_seq_sites = find_variable_sites(vcf_file_path=args.VCFFile, sequenced_sites=ddrad_seq_sites)
# print(f'found {len(variable_ddrad_seq_sites)} ddRad sites with variation')
# print(f'number of sites found matches expectation (48 sites): {len(variable_ddrad_seq_sites)==48}')
# print(f'{len(100*variable_ddrad_seq_sites)/len(ddrad_seq_sites):.3f}% of sites sequenced with ddRad had variation')

#
### MAIN CODE ###
# Once you think your code is complete, comment out the above testing section, uncomment the below section,
# and run the script.

args = my_parse_args()
target_chromosome = 'NC_036780.1'
dna_string = read_fasta(path_to_fasta=args.GenomeFile, chromosome=target_chromosome)

if args.Mode == 'SingleRad':
    print(f'running SingleRad for chromosome {target_chromosome} and restriction enzyme {args.RE1}')
    seq_sites = run_single_rad(dna_string, args.RE1)
elif args.Mode == 'ddRad':
    print(f'running ddRad for chromosome {target_chromosome} and restriction enzymes {args.RE1} and {args.RE2}')
    seq_sites = run_ddrad(dna_string, args.RE1, args.RE2)

input_len = len(dna_string)
n_sites = len(seq_sites)
sequencing_length = seq_sites[0][1] - seq_sites[0][0]
sequenced_fraction = sequencing_length*n_sites/input_len
print(f'{args.Mode} sequencing complete.')
print(f'\tlength of input sequence: {input_len}')
print(f'\tnumber of sequencing sites located: {n_sites}')
print(f'\tpercentage of nucleotides sequenced: {100*sequenced_fraction:.3f}%')

print('checking for variation within sequenced sites')
variable_sites = find_variable_sites(args.VCFFile, sequenced_sites=seq_sites)
n_var_sites = len(variable_sites)
var_fraction = n_var_sites / n_sites
print(f'\tnumber of sites with variation: {n_var_sites}')
print(f'\tfraction of sites with variation: {var_fraction}')

print('\nAnalysis Complete')


