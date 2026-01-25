#Singularity/docker
docker pull immcantation/suite:4.5.0

docker run -it immcantation/suite:4.5.0 bash

#(1) Converting 10X V(D)J data into the AIRR Community standardized format:

AssignGenes.py igblast -s filtered_contig.fasta -b /usr/local/share/igblast --organism human --loci ig --format blast --outdir results --outname BCR_data_sequences


MakeDb.py igblast -i BCR_data_sequences_igblast.fmt7 -s filtered_contig.fasta -r /usr/local/share/germlines/imgt/human/vdj/ --10x filtered_contig_annotations.csv --extended


#(2) Identifying clones from B cells in AIRR formatted 10X V(D)J data:

 #Splitting into separate light and heavy chain files:

#(a)Heavy chain:


ParseDb.py select -d BCR_data_sequences_igblast_db-pass.tsv -f locus -u "IGH" --logic all --regex --outname BCR_heavy

#(b) Light chain:

ParseDb.py select -d GBCR_data_sequences_igblast_db-pass.tsv -f locus -u "IG[LK]" --logic all --regex --outname BCR_light



#(3) Clustering sequences into clonal groups:

#Heavy chain:

DefineClones.py -d heavy_parse-select.tsv --act set -- model ham --norm len --dist 0.16


#First Analyze Heavy chain then use the following command for both heavy and light chain:
 

light_cluster.py -d BCR_heavy_parse-select_clone-pass.tsv -e BCR_light_parse-select.tsv -o BCR_heavy_light_clone_parse.tsv

#(4) Reconstructing germline sequences:

CreateGermlines.py -d heavy_parse-select_clone-pass.tsv -g dmask --cloned -r IGHV.fasta IGHD.fasta IGHJ.fasta
