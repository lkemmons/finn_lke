#!/bin/bash

cd /glade/work/emmons/FINN_python/FINNv2.9nrt
echo " $(pwd)"
echo "submitting run_finnpy_nrt_daily.pbs  $(date -Iseconds) "

/opt/pbs/bin/qsub ./run_finnpy_nrt_daily.pbs

exit 0
