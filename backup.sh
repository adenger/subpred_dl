#!/bin/bash

mkdir backup
conda env export -n subpred_deeplearning | sed -e 's/^.*subpred==5.0.0.*$/# &/' -e 's/^.*pip:*$/# &/' -e '/^prefix/d' > environment_full.yml
conda env export -n dnn_cpu | sed -e 's/^.*subpred==5.0.0.*$/# &/' -e 's/^.*pip:*$/# &/' -e '/^prefix/d' > environment_dnn_cpu_full.yml
tar --exclude=".git" --exclude="data" --exclude="build" --exclude="backup" --exclude="backup_code.tar" -cf backup/backup_code.tar .
tar cf - data/results | pigz > backup/backup_results.tar.gz
tar cf - data/tmp_data | pigz > backup/backup_tmp_data.tar.gz
tar cf - data/datasets | pigz > backup/backup_datasets.tar.gz
tar cf - data/raw | pigz > backup/backup_rawdata.tar.gz
