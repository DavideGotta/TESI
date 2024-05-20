from pathlib import Path
dir=Path("/home/davide/Desktop/genomiChro/intergeniche_RefSeq/semplici")
blastp_dir=Path("/home/davide/Desktop/genomiChro/blastp_RefSeq_CCMEE29VSALL")
genomi_dir=Path("/home/davide/Desktop/genomiChro/annotati_Refseq")
from Bio import SeqIO
import pandas as pd
for file in genomi_dir.iterdir():
    locus_to_pid={}
    if file.is_file() and file.suffix==".gbff":
        for record in SeqIO.parse(file, "genbank"):
            for feature in record.features:
                if feature.type=="CDS":
                    if "protein_id" in feature.qualifiers:
                        locus_to_pid[feature.qualifiers["locus_tag"][0]]=feature.qualifiers["protein_id"][0]
    df_locus = pd.DataFrame(locus_to_pid.items(), columns=["locus", "pid"])
    for file2 in blastp_dir.iterdir():
        if file2.is_file():
            if file.stem in file2.stem:
                df_blastp=pd.read_csv(file2,sep="\t",header=None)
                #drop second column
                df_blastp.columns=["pidCCMEE29","pid","info","pident","length","evalue"]
                #drop duplicates of pid
                df_blastp.drop_duplicates(subset="pidCCMEE29",inplace=True,keep="first")
                #drop duplicates of pid after order by evalue
                df_blastp.sort_values(by="evalue",inplace=True)
                df_blastp.drop_duplicates(subset="pid",inplace=True,keep="first")
                df=pd.merge(df_blastp,df_locus,on="pid",how="left")
                print(df_locus)
                df.to_csv(f"/home/davide/Desktop/genomiChro/blastpcsv/{file.stem}.csv",index=False)

from pathlib import Path
blastp_dir=Path("/home/davide/Desktop/genomiChro/blastpcsv/")
intergeniche_dir=Path("/home/davide/Desktop/genomiChro/intergeniche_RefSeq/semplici")
import pandas as pd
from Bio import SeqIO
for file in blastp_dir.iterdir():
    if file.is_file() and file.suffix==".csv":
        seqs=[]
        df=pd.read_csv(file)
        df=df[["locus","pidCCMEE29","info","pident","length","evalue"]]
        #rename info to info_ortologo
        df.rename(columns={"info":"info_ortologo"},inplace=True)
        df.set_index("locus",inplace=True)
        #print udplicated of index
        diz=df.to_dict(index=True,orient="index")
        for file2 in intergeniche_dir.iterdir():
            if file.stem in file2.stem:
                print(file2.stem,file.stem)
                for record in SeqIO.parse(file2, "fasta"):
                    locus=record.id
                    if locus in diz and len(record.seq)>8:
                        record.description=str(diz[locus])[1:-1]+", "+str(record.description)[record.description.find("gene"):]
                        seqs.append(record)
                SeqIO.write(seqs,f"/home/davide/Desktop/genomiChro/intergeniche_RefSeq/ortologhi/{file2.stem}.fasta","fasta")


