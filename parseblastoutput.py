from Bio.Blast import NCBIXML


def parseblastxml(input_file: str, e_value: float = 1e-5):
    """
    Parse the BLAST XML output file and extract the hits.

    Args:
        input_file (str): Path to the BLAST XML output file.

    Returns:
        list: List of dictionaries containing information about the hits.
    """
    # Open the BLAST XML output file
    with open(input_file) as blast_file:
        # Parse the BLAST XML output
        blast_records = NCBIXML.parse(blast_file)
        # Create a list to hold the hits
        hits = []
        queries = []
        # Loop through the BLAST records
        for blast_record in blast_records:

            query = blast_record.query
            queries.append(query)
            # Loop through the alignments in the BLAST record
            for alignment in blast_record.alignments:
                hit_identity = alignment.hsps[0].identities
                # Get the hit description
                hit_description = alignment.hit_def
                # Get the hit length
                hit_length = alignment.length
                # Get the hit score
                hit_score = alignment.hsps[0].score
                # Get the hit e-value
                hit_evalue = alignment.hsps[0].expect
                hit_identity = round(hit_identity / hit_length * 100, 2)
                # Create a dictionary to hold the hit information
                hit_info = f"evalue: {hit_evalue}, id:{hit_identity}%, desc:{hit_description}"

                # Add the hit information to the list
                hits.append(hit_info)

        # Return the list of hit
        diz = dict(zip(queries, hits))
        return diz


import os

dir = "/home/davide/Desktop/genomiChro/blastp_AnaCd"
import pandas as pd

df = pd.DataFrame()
for file in os.listdir(dir):
    if file.endswith(".xml"):
        path = os.path.join(dir, file)
        hits = parseblastxml(path)
        print(hits)
        hits_df = pd.DataFrame.from_dict(hits, orient='index', columns=[file[:-4]])
        # Join this DataFrame with the initial DataFrame using the file name as the column name
        df = df.join(hits_df, how='outer')
        # create dataframe with keys of dict hits as rows and add columns for every new file with column named like file
dir_genomi = "/home/davide/Desktop/genomiChro"
df.to_csv(os.path.join(dir_genomi, "AnaCd.csv"))
