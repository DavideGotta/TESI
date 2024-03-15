import pandas as pd
from bioservices import UniProt
from Bio import ExPASy
from Bio import SwissProt

class GeneAnalyzer:
    def __init__(self, gene_name):
        self.gene_name = gene_name
        self.uniprot_service = UniProt(verbose=False)

    def get_uniprot_id(self):
        results = self.uniprot_service.search(f"{self.gene_name} AND (taxonomy_id:54298)", frmt="tsv", limit=3)
        if not results:
            return None
        results_lines = results.split("\n")
        if len(results_lines) < 2:
            results = self.uniprot_service.search(self.gene_name, frmt="tsv", limit=3)
            results_lines = results.split("\n")
        try:
            protein_id = results_lines[1].split("\t")[0]
        except IndexError:
            print(f"No Uniprot ID found for {self.gene_name}.")
            protein_id = None
        return protein_id

    def fetch_uniprot_info(self, uniprot_id):
        handle = ExPASy.get_sprot_raw(uniprot_id)
        record = SwissProt.read(handle)
        record_dict = record.__dict__
        return record_dict

    def fetch_go_annotation(self, protein_id):
        if protein_id is None:
            print("No protein ID provided.")
            return
        result = self.uniprot_service.search(protein_id, frmt="tsv", columns="go_p,go_f")
        result_lines = result.split("\n")
        if len(result_lines) < 2:
            print("No GO annotation found.")
            return
        return result_lines[1]

# Uso della classe
import re
df=pd.read_csv("/home/davide/Downloads/geniconMotivoMEMEinGenomiChro.csv")
#rename first column to 'gene'
df.rename(columns={df.columns[0]: 'gene'}, inplace=True)
#add column GO annotation
df['GO annotation'] = ""
for gene in df['gene']:
    if "/" in gene:
        gene = gene.split("/")[0]
    if "_" in gene:
        gene = gene.split("_")[0]
    #create a GeneAnalyzer object
    gene_analyzer = GeneAnalyzer(gene)
    #get the uniprot id
    uniprot_id = gene_analyzer.get_uniprot_id()
    #fetch go annotation
    go_annotation = gene_analyzer.fetch_go_annotation(uniprot_id)
    #add go annotation to the dataframe
    df.loc[df['gene'] == gene, 'GO annotation'] = go_annotation
    print(f"Gene: {gene}")
    print(f"Uniprot ID: {uniprot_id}")
    print(f"GO annotation: {go_annotation}")
    print("\n")

    df.to_csv("/home/davide/Desktop/genomiChro/intergeniche_tutte/motivounicobest/geniconMotivoMEMEinGenomiChroGO.csv", index=False)
    df.to_excel("/home/davide/Desktop/genomiChro/intergeniche_tutte/motivounicobest/geniconMotivoMEMEinGenomiChroGO.xlsx", index=False)