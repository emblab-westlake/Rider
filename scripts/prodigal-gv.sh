source ~/miniconda3/bin/activate prodigal-gv

for i in /usr/commondata/public/gaoyang/human_genome_project_gut_IBD/assembly_1K/*fa
do
base=$(basename ${i})
prodigal-gv  -i ${i} \
             -a /usr/commondata/public/gaoyang/human_genome_project_gut_IBD/assembly_1K/prodigal_gv/${base}_prodigalgv.faa \
             -p meta
done