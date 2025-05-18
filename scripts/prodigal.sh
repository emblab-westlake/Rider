# source ~/miniconda3/bin/activate prodigal


# for i in /usr/commondata/public/gaoyang/Actic_metatrans/split_prok_3/*fasta
# do
# base=$(basename ${i})

# prodigal -i /usr/commondata/public/gaoyang/Actic_metatrans/split_prok_3/${base} \
#          -o /usr/commondata/public/gaoyang/Actic_metatrans/split_prok_3/prodigal/${base}.gff \
#          -a /usr/commondata/public/gaoyang/Actic_metatrans/split_prok_3/prodigal/${base}.faa \
#          -d /usr/commondata/public/gaoyang/Actic_metatrans/split_prok_3/prodigal/${base}.fna \
#          -f gff -p meta
# done

for i in /usr/commondata/public/gaoyang/yuanlin/test_for_NewVirFinder/*mlen1kb.fa
do
base=$(basename ${i})
prodigal-gv  -i ${i} \
             -a /usr/commondata/public/gaoyang/yuanlin/test_for_NewVirFinder/prodigal-gv/${base}_proteins.faa \
             -p meta
done