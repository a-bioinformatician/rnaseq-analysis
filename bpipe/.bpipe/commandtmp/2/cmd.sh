#!/bin/sh
echo "hello test"
result=$?
echo -n $result > .bpipe/commandtmp/2/cmd.exit
exit $result
