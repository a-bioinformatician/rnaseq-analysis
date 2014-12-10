#!/bin/sh
echo "hello test
result=$?
echo -n $result > .bpipe/commandtmp/1/cmd.exit
exit $result
