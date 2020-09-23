#! /bin/bash

#. $GALAXY_CONDA_PREFIX/etc/profile.d/conda.sh
#conda activate rgi
export PATH=/opt/conda/envs/rgi/bin:$PATH

card_path=$1
contig=$2
prefix=$3
heatmap=$4

/opt/conda/envs/rgi/bin/rgi load --card_json ${card_path}/card.json --local
cp ${card_path}/DB/* localDB/
/opt/conda/envs/rgi/bin/rgi main --input_sequence ${contig} --output_file ${prefix} --input_type contig --local --low_quality
# num=`wc -l out.txt|awk '{print $1}'`
#echo $num
# awk 'BEGIN{FS="\t";OFS="\t"} {if(NR==1) {print "耐药基因/蛋白名称(Best_Hit_ARO)","抗性药物名称(Drug Class)","抗性机制(Resistance Mechanism)","耐药基因家族(AMR Gene Family)"}else{print $9,$15,$16,$17}}END{a=FNR-1;print "#总共检测到"a"个耐药基因"}' out.txt > out.amr.tsv
python $(dirname $0)/amr_format.py --input_file out.txt --data_type contig
#sed -i "1i\#总共检测到${num}个耐药基因" amr.tsv 
#sed -i 's/Best_Hit_ARO/耐药基因\/蛋白名称(Best_Hit_ARO)/g' amr.tsv
#sed -i ''

##plot heatmap
/opt/conda/envs/rgi/bin/rgi heatmap --input ./ --output $heatmap -clus both
mkdir -p png_out
cp ./${heatmap}*png ./png_out/heatmap.png
