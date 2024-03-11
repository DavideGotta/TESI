from Bio import Entrez, SeqIO
import urllib
# Set your email (required by NCBI)
Entrez.email = "davide.gotta@gmail.com"

# Search for Chroococcidiopsis genomes
search_handle = Entrez.esearch(db="assembly", term="Chroococcidiopsis", retmax=500)
search_results = Entrez.read(search_handle)
search_handle.close()

# Get the list of genome IDs
genome_ids = search_results["IdList"]

# Fetch the summary records of the genome IDs
summary_handle = Entrez.esummary(db="assembly", id=",".join(genome_ids))
summaries = Entrez.read(summary_handle)
summary_handle.close()
print(summaries)
print(type(summaries))
print(len(summaries['DocumentSummarySet']['DocumentSummary']))
print(summaries['DocumentSummarySet']['DocumentSummary'][0]['AssemblyAccession'])
# Print the names associated with the IDs

for i in range(len(summaries['DocumentSummarySet']['DocumentSummary'])):
    assembly_accession = summaries['DocumentSummarySet']['DocumentSummary'][i]['AssemblyAccession']
    organism = summaries['DocumentSummarySet']['DocumentSummary'][i]['Organism'].replace(' ', '_')
    ftp_path = summaries['DocumentSummarySet']['DocumentSummary'][i]['FtpPath_GenBank']
    #downlaod fasta
    filename = ftp_path + "/" + ftp_path.split('/')[-1] + '_genomic.fna.gz'
    localpath = f"/home/davide/Downloads/genomiChro/{organism}_{assembly_accession}.fna.gz"
    urllib.request.urlretrieve(filename, localpath)







