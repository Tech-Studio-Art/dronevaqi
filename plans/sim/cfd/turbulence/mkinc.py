"""Make flat include directories (like OpenFOAM lnInclude) from the v1912 source tree, using symlinks."""
import os, sys
src, out = sys.argv[1], sys.argv[2]
libs = ["OpenFOAM", "OSspecific/POSIX", "finiteVolume", "meshTools", "fileFormats", "surfMesh",
        "TurbulenceModels/turbulenceModels", "TurbulenceModels/compressible", "transportModels/compressible",
        "thermophysicalModels/basic", "thermophysicalModels/radiation", "thermophysicalModels/specie",
        "thermophysicalModels/solidThermo", "thermophysicalModels/solidSpecie", "thermophysicalModels/thermophysicalProperties"]
for lib in libs:
    d = os.path.join(out, lib.replace("/", "_"))
    os.makedirs(d, exist_ok=True)
    n = 0
    for root, dirs, files in os.walk(os.path.join(src, lib)):
        dirs[:] = [x for x in dirs if x != "lnInclude" and not x.startswith("Make")]
        for f in files:
            if f.endswith((".H", ".C", ".hxx", ".h", ".cxx")):
                dst = os.path.join(d, f)
                if not os.path.exists(dst):
                    os.symlink(os.path.join(root, f), dst)
                    n += 1
    print(lib, n)
