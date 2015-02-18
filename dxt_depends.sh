#!/bin/sh
SOURCE_PREFIX="2"
DEST_PREFIX="$3"
DXT_PATH=$(dirname "$1")
BUILD_RULE=$(gcc -MM -MG -Isrc -Isrc/linearization -Isrc/LLDLA -Isrc/DLA/ -Isrc/tensors/ "$2/$1")

if [ $DXT_PATH = "." ]
then
BUILD_RULE="$3/${BUILD_RULE}"
else
BUILD_RULE="$3/${DXT_PATH}/${BUILD_RULE}"
fi

echo "${BUILD_RULE}"
echo -e "\t\$(CC) \$(CFLAGS) -c \$< -o \$@"
