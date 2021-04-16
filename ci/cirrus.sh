#!/bin/sh

set -e
set -x

export LC_ALL=C

env >> test_env.log

$CC -v || true
valgrind --version || true

./autogen.sh

./configure --host="$HOST" --enable-benchmark="$BENCH"

# We have set "-j<n>" in MAKEFLAGS.
make

# Print information about binaries so that we can see that the architecture is correct
file test* || true
file bench* || true
file .libs/* || true

if [ -n "$BUILD" ]
then
    make "$BUILD"
fi

if [ -n "$EXEC_CMD" ]; then
    $EXEC_CMD ./test $TESTRUNS
    $EXEC_CMD ./test-verify $TESTRUNS
fi

if [ "$BENCH" = "yes" ]; then
    $EXEC_CMD ./bench
fi
