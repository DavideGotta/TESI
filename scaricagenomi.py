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

print(len(summaries['DocumentSummarySet']['DocumentSummary']))
# Print the names associated with the IDs
from pprint import pprint

for i in range(len(summaries['DocumentSummarySet']['DocumentSummary'])):
    assembly_accession = summaries['DocumentSummarySet']['DocumentSummary'][i]['AssemblyAccession']
    organism = summaries['DocumentSummarySet']['DocumentSummary'][i]['Organism'].replace(' ', '_')
    ftp_path = summaries['DocumentSummarySet']['DocumentSummary'][i]['FtpPath_RefSeq']
    print(organism,assembly_accession,ftp_path)
    #downlaod fasta
    filename = ftp_path + "/" + ftp_path.split('/')[-1] + '_protein.faa.gz'
    print(filename)
    #filename = ftp_path + "/" + ftp_path.split('/')[-1] + '_genomic.fna.gz'
    localpath = f"/home/davide/Desktop/genomiChro/{organism}_{assembly_accession}.faa.gz"
    if ftp_path:
        urllib.request.urlretrieve(filename, localpath)







