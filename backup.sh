#!/bin/bash

mkdir backup
tar cf - data/raw | pigz > backup/backup_rawdata.tar.gz
tar cf - data/datasets | pigz > backup/backup_datasets.tar.gz
tar --exclude=".git" --exclude="data" --exclude="build" --exclude="backup_code.tar" -cf backup/backup_code.tar .
