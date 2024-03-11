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
        protein_id = results_lines[1].split("\t")[0]
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
with open("genifimo.txt", "r") as f:
    lines = f.readlines()
    for line in lines:
        #gene = re.search(r'gene=(\w+)', line).group(1)
        gene=line.split()[0]
        print(gene)
        if "/" in gene:
            gene = gene.split("/")[0]
        if "_" in gene:
            gene = gene.split("_")[0]
        analyzer = GeneAnalyzer(gene)
        uniprot_id = analyzer.get_uniprot_id()
        go_annotation = analyzer.fetch_go_annotation(uniprot_id)
        with open("go_annotationsFIMO.txt", "a") as f:
            f.write(f"{gene} {go_annotation}\n")
