categorizes the mapped fragments into genomic categories based on the compatibility with the features defined by gene annotations. Specifically, it categorises the reads into the following categories : CDS, UTR5, UTR5, splice junction, introns, intergenic, deep intergenic, rRNA, MT and multiMapped

To run it simultaneosly for many bam files: 1) make a file with sample name

ls *bam | awk -F "." '{print $1}'>samples.txt

2) while read line ; do echo "~/anaconda/bin/python ~/code2/gprofile/gprofile.py $PWD/../bamPCRDublicatesRemoved2/${line}.bam $PWD/${line}/ $PWD/${line}/genomicFeatures.txt h">run${line}.sh; done<../samples.txt

Every 'run' script contains the comand to run

3) Use qsub to run for all samples ls run*sh | awk '{i+=1;print "qsub -cwd -V -N feature"i" -l h_data=8G,time=03:00:00 "$1}' > all.sh chmod 755 all.sh nohup ./all.sh &
