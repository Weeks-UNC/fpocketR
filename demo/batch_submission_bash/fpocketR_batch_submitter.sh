source $(conda info --base)/etc/profile.d/conda.sh
conda activate fpocketR
while read line; do echo "$line"; python -m fpocketR $line; done < fpocketR_batch_file.txt