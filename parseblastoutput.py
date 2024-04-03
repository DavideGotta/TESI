from Bio.Blast import NCBIXML
def parseblastxml(input_file:str,e_value:float=1e-5)->list:
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
        # Loop through the BLAST records
        for blast_record in blast_records:
            # Loop through the alignments in the BLAST record
            for alignment in blast_record.alignments:
                print(dir(alignment))

                hit_identity = alignment.hsps[0].identities
                # Get the hit description
                hit_description = alignment.title
                # Get the hit length
                hit_length = alignment.length
                # Get the hit score
                hit_score = alignment.hsps[0].score
                # Get the hit e-value
                hit_evalue = alignment.hsps[0].expect
                # Create a dictionary to hold the hit information
                hit_info = {
                    "description": hit_description.split()[1],
                    "length": hit_length,
                    "score": hit_score,
                    "evalue": hit_evalue,
                    "identity": f"{hit_identity/hit_length*100}.2f%"
                }
                # Add the hit information to the list
                hits.append(hit_info)
        # Return the list of hits
        return hits
hits=parseblastxml("/home/davide/Desktop/genomiChro/blastp_AnaCd/Chroococcidiopsis_sp._CCMEE_29_GCF_023558375.1.xml")
print(hits[:5])