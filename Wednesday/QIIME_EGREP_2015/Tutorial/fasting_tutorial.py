from qiime_wrappers import *
import os
from Bio import SeqIO
from subprocess import Popen

output_directory = '2015_fasting_tutorial/'

# Input file names

input_fasta = './Fasting_Example.fna'
input_qual = './Fasting_Example.qual'

fasta_ref = ('./18S_tutorial_files/'+
            '18S_tutorial_sample_seqs.fna')

taxa_ref = ('./18S_tutorial_files/'+
            '18S_tutorial_sample_seqs.txt')

mapping = 'Fasting_Map.txt'

ref_seq_aln = ('./18S_tutorial_files/'+
            '18S_tutorial_sample_seqs_aln.fna')

# a list containing addresses of output files
file_urls = []

def add_urls(to_add, file_urls=file_urls):
    if isinstance(to_add, list):
        file_urls += to_add
    elif isinstance(to_add,str):
        file_urls.append(to_add)
    return file_urls

# Check mapping files

print 'Checking mapping files ...'

outputpath = "validate_mapping"
output = validate_mapping_file(mapping,output_directory+outputpath)
file_urls = add_urls(output)


# Demultiplex the fasta files
print 'Sorting sequence reads by sample tags ...'

barcode=str(12)
outputpath='SplitReads/'
output = split_libraries(mapping, input_fasta, input_qual, barcode, output_directory+outputpath)
file_urls = add_urls(output)

# Picking Operational Taxonomic Units (OTUs) through making OTU table

print 'Clustering reads to OTUs ...'

demultiplexed_fasta = output_directory+'SplitReads/seqs.fna'
output_path = output_directory+'PickedOtus'

pick_otus(demultiplexed_fasta, output_path)


# Pick representative set of sequences
print 'Picking one representatinve sequence read for each OTU ...'

seqs_otus = output_directory+'PickedOtus/seqs_otus.txt'
rep_set_fasta = output_directory+'rep_set.fasta'
pick_rep_set(seqs_otus, demultiplexed_fasta, rep_set_fasta)


# Assign taxonomy to each representative sequence

print 'Matching representative reads with taxonomically identified reference sequences ...'

output_path = output_directory+'AssignedTaxonomy'
rep_set_fasta = output_directory+'rep_set.fasta'
assign_taxonomy(rep_set_fasta, taxa_ref, fasta_ref, output_path, rdp_max_memory=15000)

# Align sequences using mafft

print 'Aligning representative sequences ...'

output_path = output_directory+'AlignedRepSet'
os.mkdir(output_path)

cline = 'mafft %s > %s/%s'%(rep_set_fasta, output_path,'rep_set_aligned.fasta')
p = Popen(cline, shell=True)
out, err = p.communicate()


# Filter with trimal gappyout

print 'Filtering out ambigously aligned positions ...'

aln = output_directory+'AlignedRepSet/rep_set_aligned.fasta'
filtered_aln = output_directory+'AlignedRepSet/rep_set_aligned_filtered.fasta'

filter_aln_with_trimal(aln, filtered_aln)

filtered_corrected_aln = (output_directory+
                          'AlignedRepSet/'+
                          'rep_set_aligned_filtered.fasta.corrected.fasta')


# Make Phylogeny

print 'Reconstructing a phylogenetic tree of the representative sequences ...'

output_path = output_directory+'AlignedRepSet/rep_set.tre'

make_phylogeny(filtered_corrected_aln, output_path)

# Print Summary

print 'Writing OTU summary ...'

assigned_taxa = output_directory+'AssignedTaxonomy/rep_set_tax_assignments.txt'
output_biom = output_directory+'otu_table.biom'
output_path = output_directory+'summary_table.txt'
make_otu_table(seqs_otus, assigned_taxa, output_biom, output_path)
file_urls = add_urls(output_path)

# Make OTU Heatmap

print 'Writing OTU hitmap ...'

output_path = output_directory+'OTU_heatmap'
output = make_otu_heatmap_html(output_biom, output_path)
file_urls = add_urls(output)

# Beta Diversity

print 'Calculating beta diversity (PCA of samples using taxa as variables) ...'

beta_diversity_path = output_directory+'BetaDiversity'
rep_set_tree = output_directory+'AlignedRepSet/rep_set.tre'
output = beta_diversity_through_plots(output_biom, mapping,
                                      beta_diversity_path, rep_set_tree)

file_urls = add_urls(output)


# Summarize Communities by Taxonomic Composition

print 'Writing community taxonomic summary ...'

output_path = output_directory+'CommunitySummary'
output = summarize_taxa_through_plots(output_biom, output_path, mapping)

file_urls = add_urls(output)

# Print html file

print 'Making main HTML file in %s'%output_directory[:-1]+'_main_results_file.html'

def write_html_link(address,text):
    return "<a href=%s >%s</a>"%(address,text)

html = open(output_directory[:-1]+'_main_results_file.html','wt')

for a in file_urls:
    html.write(write_html_link(a,a)+'<br><br>')
html.close()

