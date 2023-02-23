#!/bin/bash

LIST_JOBS=`bjobs | awk '{print $1}' | grep -v JOBID`

for job in $LIST_JOBS;
do
        bkill $job
done
