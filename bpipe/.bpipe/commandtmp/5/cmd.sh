#!/bin/sh
samtools
result=$?
echo -n $result > .bpipe/commandtmp/5/cmd.exit
exit $result
