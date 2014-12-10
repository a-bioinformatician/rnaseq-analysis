#!/bin/sh
samtools
result=$?
echo -n $result > .bpipe/commandtmp/4/cmd.exit
exit $result
