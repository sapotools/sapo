#!/usr/bin/env bash

[ -f ./scripts/run-tests.bash ] || cmake .

echo "Running pre-commit hook"
sh ./scripts/run-tests.bash

# $? stores exit value of the last command
if [ $? -ne 0 ]; then
 echo "Tests must pass before commit!"
 exit 1
fi
