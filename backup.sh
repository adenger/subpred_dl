#!/bin/bash

tar cf - data/raw | pigz > backup_rawdata.tar.gz
tar cf - data/datasets | pigz > backup_datasets.tar.gz
