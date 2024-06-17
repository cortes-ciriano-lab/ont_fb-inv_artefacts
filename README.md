

```
# A) Parse arguments
sample_id=$1
bam=$2
sample_type=$3
output_dir=$4
echo "sample_type: $sample_type"
echo "sample_id: $sample_id"
echo "bam: $bam"
echo "output_dir: $output_dir"


# B) Run code
for i in {1..22};
do

chromosome="chr"$i
python fb-inv_artefact_rates.py $bam $chromosome $sample_type $sample_id $output_dir

done
```
