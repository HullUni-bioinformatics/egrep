from qiime_wrappers import *
import os
from Bio import SeqIO

output_directory = '2015_EGREP_QIIME_results/'

# Input file names

input_data = {'plate1':{'fasta': 'Plate1_1_subset.fasta.gz', 'qual': 'Plate1_1_subset.qual.gz', 'map': '.MappingFile.txt'},
              'plate2':{'fasta': 'Plate2_1_subset.fasta.gz', 'qual': 'Plate2_1_subset.qual.gz', 'map': '.2MappingFile.txt'},
              'plate3':{'fasta': 'Plate3_1_subset.fasta.gz', 'qual': 'Plate3_1_subset.qual.gz', 'map': '.3MappingFile.txt'},
              'plate4':{'fasta': 'Plate4_1_subset.fasta.gz', 'qual': 'Plate4_1_subset.qual.gz', 'map': '.4MappingFile.txt'},
}
taxa_ref = ('./EUK_MT-CO1_refseq_24SEP15/'+
            '7level_taxonomy.txt')
fasta_ref = ('./EUK_MT-CO1_refseq_24SEP15/'+
            'ref_seq.fasta')
ref_seq_aln = ("./EUK_MT-CO1_refseq_24SEP15/"+
               "ref_seq_aln.fasta")

mapping = '.CombinedMappingFile.txt'

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

for plate in input_data:
    print plate
    outputpath = "%s_validate_mapping"%plate
    output = validate_mapping_file(input_data[plate]['map'],output_directory+outputpath)
    file_urls = add_urls(output)


# Demultiplex the fasta files
print 'Sorting sequence reads by sample tags ...'

barcode_length=str(8)
primer_mismatch = str(8)

for plate in input_data:
    print plate
    outputpath = "%s_demultiplexed_fasta"%plate
    fasta = input_data[plate]['fasta']
    qual = input_data[plate]['qual']
    plate_map = input_data[plate]['map']
    output = split_libraries(plate_map, fasta, qual, barcode_length,
                             output_directory+outputpath,
                             M=primer_mismatch)
    file_urls = add_urls(output)

# Weld fastas

records = []
length = 210
print "Welding demultiplexed fastas and trimming to %i bp"%length
os.mkdir(output_directory+'SplitReads')
demultiplexed_fasta = output_directory+'SplitReads/seqs.fna'
for plate in input_data:
    records += list(SeqIO.parse(output_directory+"%s_demultiplexed_fasta"%plate+'/seqs.fna',
                                'fasta'))
if length:
    for i in range(len(records)):
        records[i].seq = records[i].seq[:210]
SeqIO.write(records, demultiplexed_fasta,'fasta')

# Picking Operational Taxonomic Units (OTUs) through making OTU table

print 'Clustering reads to OTUs ...'

demultiplexed_fasta = output_directory+'SplitReads/seqs.fna'
output_path = output_directory+'PickedOtus'
sim = str(0.96)
pick_otus(demultiplexed_fasta, output_path, similarity=sim)

# Filter out OTUs with less than 3 sequences

seqs_otus = output_directory+'PickedOtus/seqs_otus.txt'
filtered_seqs_otus = output_directory+'PickedOtus/seqs_otus_filtered.txt'
lines = open(seqs_otus,'r').readlines()
count = 0
with open(filtered_seqs_otus,'wt') as hndl:
    for l in lines:
        if l.count('\t') > 2:
            hndl.write(l)
            count += 1
print 'There are %i OTUs with more than 2 sequences\n'%count
# Pick representative set of sequences
print 'Picking one representatinve sequence read for each OTU ...'

seqs_otus = filtered_seqs_otus
rep_set_fasta = output_directory+'rep_set.fasta'
pick_rep_set(seqs_otus, demultiplexed_fasta, rep_set_fasta)


# Assign taxonomy to each representative sequence

print 'Matching representative reads with taxonomically identified reference sequences ...'

output_path = output_directory+'AssignedTaxonomy'
rep_set_fasta = output_directory+'rep_set.fasta'
assign_taxonomy(rep_set_fasta, taxa_ref, fasta_ref, output_path, rdp_max_memory=15000)

# Align sequences using pynast

print 'Aligning representative sequences ...'

output_path = output_directory+'AlignedRepSet'

align_seqs(rep_set_fasta, ref_seq_aln, output_path)

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

# ANOSIM

print 'Computing ANOSIM for the Forest categories and for Species categories ...'

output_path = output_directory+'ANOSIM'
forest = compare_categories(beta_diversity_path, mapping, 'Forest', output_path+'_Forest')
species = compare_categories(beta_diversity_path, mapping, 'Species', 
                              output_path+'_Species')

file_urls = add_urls(forest)
file_urls = add_urls(species)

# Summarize Communities by Taxonomic Composition

print 'Writing community taxonomic summary ...'
output_path = output_directory+'CommunitySummary'
output = summarize_taxa_through_plots(output_biom, output_path, mapping)

file_urls = add_urls(output)

# Compute Alpha Diversity within the Samples and Generate Rarefaction Curves

print 'Computing alpha diversity (rarefaction curves) ...'

rep_set_tree = output_directory+'AlignedRepSet/rep_set.tre'
output_biom = output_directory+'otu_table.biom'
output_path = output_directory+'AlphaDiversity'
output = alpha_rarefaction(output_biom, mapping, output_path, rep_set_tree)
file_urls = add_urls(output)

# Print html file
if False:"""
print 'Repeating all stages for Metazoa only. Starting by extracting Insecta OTUs from the rep set and sequence alignments'

# <codecell>

extract_subset_OTUs(assigned_taxa, assigned_taxa+'.Insecta.txt', 'Insecta')

# <codecell>

extract_subset_seqs_otus(assigned_taxa+'.Insecta.txt', seqs_otus, seqs_otus+'.Insecta.txt')

# <codecell>

aln = output_directory+'AlignedRepSet/rep_set_aligned_filtered.fasta.corrected.fasta'
subset_aln = output_directory+'AlignedRepSet/rep_set_aligned_filtered.Insecta.corrected.fasta'

extract_subset_aln(assigned_taxa, 'Insecta', aln, subset_aln)

# <codecell>


output_path = output_directory+'AlignedRepSet/Insecta_rep_set.tre'
filtered_corrected_aln = output_directory+'AlignedRepSet/rep_set_aligned_filtered.Insecta.corrected.fasta'

make_phylogeny(filtered_corrected_aln, output_path)


assigned_taxa = output_directory+'AssignedTaxonomy/rep_set_tax_assignments.txt.Insecta.txt'
output_biom = output_directory+'otu_table_Insecta.biom'
output_path = output_directory+'summary_table_Insecta.txt'
seqs_otus = output_directory+'PickedOtus/seqs_otus.txt.Insecta.txt'

make_otu_table(seqs_otus, assigned_taxa, output_biom, output_path)

file_urls = add_urls(output_path)

output_path = output_directory+'OTU_heatmap_Insecta'

output = make_otu_heatmap_html(output_biom, output_path)

file_urls = add_urls(output)

beta_diversity_path = output_directory+'BetaDiversityInsecta'
rep_set_tree = output_directory+'AlignedRepSet/Insecta_rep_set.tre'

output = beta_diversity_through_plots(output_biom, mapping,
                                      beta_diversity_path, rep_set_tree)

file_urls = add_urls(output)

output_path = output_directory+'InsectaANOSIM'

forest = compare_categories(beta_diversity_path,mapping, 'Forest', 
                               output_path+'_Forest')
species = compare_categories(beta_diversity_path,mapping, 'Species', 
                              output_path+'_Species')

file_urls = add_urls(forest)
file_urls = add_urls(species)

output_path = output_directory+'CommunitySummaryInsecta'

output = summarize_taxa_through_plots(output_biom, output_path, mapping)

file_urls = add_urls(output)

output_path = output_directory+'AlphaDiversityInsecta'

output = alpha_rarefaction(output_biom, mapping, output_path, rep_set_tree)


file_urls = add_urls(output)
"""
print 'Making main HTML file in %s'%output_directory[:-1]+'_main_results_file.html'

def write_html_link(address,text):
    return "<a href=%s >%s</a>"%(address,text)

html = open(output_directory[:-1]+'_main_results_file.html','wt')

for a in file_urls:
    html.write(write_html_link(a,a)+'<br><br>')
html.close()

