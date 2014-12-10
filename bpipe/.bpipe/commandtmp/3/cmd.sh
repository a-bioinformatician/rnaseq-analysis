#!/bin/sh
samtools --version
result=$?
echo -n $result > .bpipe/commandtmp/3/cmd.exit
exit $result
