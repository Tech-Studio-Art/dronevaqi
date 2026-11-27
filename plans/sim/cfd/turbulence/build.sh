#!/bin/bash
# Build libbuoyantKEpsilonR.so against the installed OpenFOAM v1912 (Ubuntu package openfoam 1912.200626).
# The OpenFOAM headers come from the matching source archive (Ubuntu pool: openfoam_1912.200626.orig.tar.bz2).
# usage: build.sh <OpenFOAM-v1912 src dir> <output lib dir>
set -e
SRC=$1; OUT=$2
HERE=$(cd "$(dirname "$0")" && pwd)
INC=$OUT/inc
mkdir -p "$OUT"
python3 -I "$HERE/mkinc.py" "$SRC" "$INC" > /dev/null
I=""
for d in OpenFOAM OSspecific_POSIX finiteVolume meshTools fileFormats surfMesh TurbulenceModels_turbulenceModels \
         TurbulenceModels_compressible transportModels_compressible thermophysicalModels_basic \
         thermophysicalModels_radiation thermophysicalModels_specie thermophysicalModels_solidThermo \
         thermophysicalModels_solidSpecie thermophysicalModels_thermophysicalProperties; do I="$I -I$INC/$d"; done
FLAGS="-std=c++11 -m64 -pthread -DOPENFOAM=1912 -DWM_DP -DWM_LABEL_SIZE=32 -O3 -DNoRepository -ftemplate-depth-100 -fPIC \
       -Wno-unused-parameter -Wno-invalid-offsetof -Wno-attributes -Wno-unknown-pragmas -Wno-deprecated-copy"
g++ $FLAGS -I"$HERE" $I -c "$HERE/makeBuoyantKEpsilonR.C" -o "$OUT/makeBuoyantKEpsilonR.o"
g++ $FLAGS -shared -Xlinker --add-needed -Xlinker --no-as-needed -o "$OUT/libbuoyantKEpsilonR.so" "$OUT/makeBuoyantKEpsilonR.o" \
    -L/usr/lib -lcompressibleTurbulenceModels -lturbulenceModels -lfluidThermophysicalModels -lspecie \
    -lfiniteVolume -lmeshTools -lOpenFOAM -ldl -lm
echo "built $OUT/libbuoyantKEpsilonR.so"
