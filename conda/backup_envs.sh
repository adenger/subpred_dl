#!/bin/bash
conda env export --no-builds -n subpred_deeplearning | sed -e 's/^.*subpred==5.0.0.*$/# &/' -e 's/^.*pip:*$/# &/' -e '/^prefix/d' > conda/environment_full.yml
conda env export  -n subpred_deeplearning | sed -e 's/^.*subpred==5.0.0.*$/# &/' -e 's/^.*pip:*$/# &/' -e '/^prefix/d' > conda/environment_full_with_builds.yml
conda env export --no-builds -n dnn_cpu | sed -e 's/^.*subpred==5.0.0.*$/# &/' -e '/^prefix/d' > conda/environment_dnn_cpu.yml
conda env export  -n dnn_cpu | sed -e 's/^.*subpred==5.0.0.*$/# &/' -e '/^prefix/d' > conda/environment_dnn_cpu_with_builds.yml
