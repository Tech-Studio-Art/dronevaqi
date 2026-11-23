# Drone V AQI – CFD and pollution-model suite

Everything here is reproducible on a Linux box with `apt install openfoam` (Debian/Ubuntu
package, OpenFOAM v1912), Python 3, `numpy`, `pyvista`, `matplotlib`.

| Script | What it simulates | Method |
|---|---|---|
| `drone3d.py` | 3-D RANS of a 2 kg quad and a 50 kg MTOW (25 kg payload) heavy-lift drone, no filter / ducted filter above / below rotors | actuator disks, k-ω SST, Darcy porous filter, PM tracer with filter sink |
| `finalize_drones.py` | CADR bookkeeping for the drone cases | mass balance (see docstring) |
| `mixing_axi.py` | Solution 1 – 100 drones as inversion mixers over 500 × 500 m | axisymmetric unit-cell URANS, Boussinesq buoyancy, k-ε |
| `cases2d.py` | Solution 2 (2-D shelter), 3a (electrostatic fence), 4a (fog net) | 2-D RANS, Richards & Hoxey ABL inlet, porous screens, Deutsch–Anderson / Langmuir sinks |
| `shelter3d.py` | Solution 2 in 3-D with open shelter ends | 3-D RANS |
| `models0d.py` | 3b vehicle ram-air filters, 4b bipolar-ion coagulation, 4c charged mist | Bernoulli + Davies, Smoluchowski/Fuchs, Slinn + Kraemer–Johnstone |
| `report.py` | figures in `figures/` | |

`REPORT.md` has the results and the literature behind every model choice.

Note: the Debian OpenFOAM v1912 build crashes on any runtime functionObject (SHA1 digest
bug), so PM is carried as the scalar `T` of `buoyantBoussinesq*Foam` with β = 0
(a pure passive tracer, Sc_t = Pr_t = 0.7), and post-processing goes through `foamToVTK`.
Case directories under `runs/` are not committed (large); the JSON summaries are.
