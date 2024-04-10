import subprocess
intergen_dir="/home/davide/Desktop/genomiChro/intergeniche_RefSeq"
motivo="/home/davide/Downloads/motivo16.meme"
import os
for file in os.listdir(intergen_dir):
    if file.endswith(".fasta"):
        #rename file to remove (cyanobacteria) from the name
        new_name=file.replace("_(cyanobacteria)","")
        os.rename(os.path.join(intergen_dir,file),os.path.join(intergen_dir,new_name))
        #create directory for each file in intergen_dir/fimo
        dir_name=os.path.join(intergen_dir,"fimo",file[:-6])
        if not os.path.exists(dir_name):
            os.mkdir(dir_name)
        command = ["fimo", "--oc", dir_name, "--verbosity", "1", "--bgfile", "--nrdb", "--thresh", "1.0E-4", "--norc", motivo, os.path.join(intergen_dir,file)]
        command = f"fimo --oc {dir_name} --verbosity 1 --bgfile --nrdb --thresh 1.0E-4 --norc {motivo} {os.path.join(intergen_dir,file)}"
        subprocess.run(command,shell=True)