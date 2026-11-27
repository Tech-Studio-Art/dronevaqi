#!/usr/bin/env python3
"""
Build report 2 (three-dimensional drone fleet in wind):
  data_entry.json  - one Research Dex entry {"title": ..., "body": [...]}
  report.html      - the same body, styled for print
  fleet_report.pdf - A4 PDF (<= 15 pages), made with pdf.js
All numbers come from plans/sim/cfd/runs/fleet3d*/ (summary.json, cells.json, calm_rhythm.json, runinfo.json).
Text follows ASD-STE100 rules.
"""
import json
import math
import os
import subprocess

HERE = os.path.dirname(os.path.abspath(__file__))
SIM = os.path.dirname(HERE)
CFD = os.path.join(SIM, "cfd")


def J(*p):
    return json.load(open(os.path.join(*p)))


S = J(CFD, "runs", os.environ.get("FLEET_SUMMARY_DIR", "fleet3d"), "summary.json")
S2 = J(CFD, "runs", os.environ.get("FLEET_SUMMARY_DIR2", "fleet3d_u2"), "summary.json")
CELL = J(CFD, "runs", "fleet3d", "cells.json")
RHY = J(CFD, "runs", "fleet3d", "calm_rhythm.json")["cases"]
RI = J(CFD, "runs", "fleet3d", "runinfo.json")
RI2 = J(CFD, "runs", "fleet3d_u2", "runinfo.json")
refs = J(HERE, "refs_verified.json")
R = {r["key"]: i + 1 for i, r in enumerate(refs)}
C, C2 = S["cases"], S2["cases"]


def c(*keys):
    return "[" + ", ".join(f"<a href='#ref{R[k]}'>{R[k]}</a>" for k in keys) + "]"


def fig(n, src, cap):
    return [f"<img src='assets/{src}'>", f"<a href='#fig{n}' id='fig{n}'>Figure {n}</a>: {cap}"]


def table(rows, head=True):
    out = "<table>"
    for i, r in enumerate(rows):
        tag = head and i == 0
        out += "<tr>" + "".join(f"<td><p>{'<strong>' + str(x) + '</strong>' if tag else x}</p></td>" for x in r) + "</tr>"
    return out + "</table>"


def f0(x):
    return f"{x:.0f}"


def f1(x):
    v = f"{x:.1f}"
    return "0.0" if v == "-0.0" else v


def place(case):
    """Sentence parts for the best position of the target area."""
    x0 = case["best_window_x0"]
    if x0 > 12:
        return (f"starts {f0(x0)} m downwind of the upwind edge of the zone", f"approximately {f0(x0)} m upwind of the "
                f"upwind edge of the area that we want to protect")
    return ("starts at the clusters", "at the upwind edge of the area that we want to protect")


# ---------------------------------------------------------------- numbers
cell = CELL["cell_1col"]["decrease_30min_pct"]              # periodic cell, fleet mesh, calm, original model
cell_r = CELL["cell_1col_rke"]["decrease_30min_pct"]        # same, final model
resolved = CELL["resolved_axi"]["decrease_30min_pct"]       # report 1 best estimate (resolved disk)
big = CELL["cell"]["series"]
big5 = 100 * (1 - big[3]["z<2m"] / big[0]["z<2m"])
pF, pC, pR6, pR60 = (RHY[k] for k in ("pF100", "pC50", "pR600", "pR60"))
K1, K4, F100 = C["K1"], C["K4"], C["F100"]
K1b = C2.get("K1")
F100b = C2["F100"]
U60, U60b = S["U60"], S2["U60"]
NBV = math.sqrt(9.81 / 285.0 * 5.0 / 100.0)
FM1 = 50 * 9.81 / (math.pi * 1.25)


def briggs(n, u):
    """Briggs bent-over momentum rise in stable air: 1.5 (Fm / (U N))^(1/3)"""
    return 1.5 * (n * FM1 / (u * NBV)) ** (1 / 3)


b1, b1b, b25, b100, b100b = briggs(1, U60), briggs(1, U60b), briggs(25, U60), briggs(100, U60), briggs(100, U60b)
pen = lambda case: 60 - case["mix_bottom_m"]
reach = lambda case: case["mix_bottom_m"] <= 1.5
disp_pen = [pen(C[k]) for k in ("F100", "C50", "R60", "R600") if k in C] + [pen(F100b)]
refchk = max(abs(v) for v in S["reference_check"]) if S.get("reference_check") else 0.0
POWER_KW = 6.77                                             # hover power of one heavy drone (report 1)
cpu = sum(v["cpu_core_h"] for v in RI.values() if isinstance(v, dict)) + \
    sum(v["cpu_core_h"] for v in RI2.values() if isinstance(v, dict))


def xpeak(case, src=S):
    pr = case["profile_x"]
    i = max(range(len(pr)), key=lambda j: pr[j])
    return src["x"][i], pr[i]


def at(case, xx, src=S):
    xs = src["x"]
    i = min(range(len(xs)), key=lambda j: abs(xs[j] - xx))
    return case["profile_x"][i]


xk1, rk1 = xpeak(K1)
xk4, rk4 = xpeak(K4)


def plume_end(case, src=S, level=5.0):
    """Largest x (m) where the along-wind street-level decrease is more than `level` %."""
    xs = [xx for xx, v in zip(src["x"], case["profile_x"]) if v > level]
    return max(xs) if xs else 0.0


Lk1, Lk4 = plume_end(K1), plume_end(K4)
zmix_min = min(C[k]["mix_bottom_m"] for k in ("F100", "C50", "R60", "R600"))


def out_share(case):
    """Share (%) of the total street-level benefit that is outside the zone."""
    return 100 * (1 - case["benefit_zone_equiv"]["zone"] / case["benefit_total"]) if case["benefit_total"] > 0 else 0.0

TITLE = ("Drone V AQI report 2: a fleet of mixing drones in wind. Three-dimensional simulation of zone edges, "
         "downwind benefit, drone placement and rhythm")
B = []
A = B.append
A(f"<a href='#title' id='title'>Title</a>: {TITLE}")
A("<a href='#author' id='author'>Author</a>: Noorul Ali, Tech Studio .Art (Drone V AQI programme), India")
A("<a href='#contact' id='contact'>Correspondence</a>: <a href='mailto:noorulali78@gmail.com'>noorulali78@gmail.com</a>")
A("<a href='#abstract'>Abstract</a><br><a href='#introduction'>Introduction</a><br><a href='#methods'>Methods</a><br>"
  "<a href='#results'>Results</a><br><a href='#discussion'>Discussion</a><br><a href='#references'>References</a><br>"
  "<a href='#appendix'>Appendix</a>")

# ---------------------------------------------------------------- abstract
A("<a href='#abstract' id='abstract'>Abstract</a>")
A(f"Report 1 calculated the effect of inversion mixing with one drone in a unit cell. A unit cell represents the "
  f"centre of a very large fleet in calm air. It has no wind and no zone edges. In this report, we simulated a full "
  f"fleet of heavy drones over a 500 m × 500 m zone in three dimensions, with wind.")
A(f"Wind changes the result completely. In a light wind of 1 m/s at 10 m height, the jet of each drone bends over. "
  f"The jets mix only the air from approximately {f0(zmix_min)} m to 70 m. With drones on a 50 m grid, street-level PM did not "
  f"change. This was true for 50 and 100 drones and for a wind of 1 m/s and 2 m/s. The Briggs equation for jets in a stable "
  f"cross-wind gives the same result.")
A(f"When the drones fly in a tight cluster, their jets join and make one strong jet. "
  + (f"One cluster of 100 drones decreased street-level PM by {f1(K1['zone_last'])} % over the zone and by up to "
     f"{f0(K1['max_ground'])} % downwind of the cluster. " if reach(K1) or K1['zone_last'] > 1 else
     f"One cluster of 100 drones decreased street-level PM by only {f1(K1['zone_last'])} % over the zone. ")
  + f"Four clusters of 25 drones gave {f1(K4['zone_last'])} %. The benefit is largest within approximately 200 m "
  f"downwind of each cluster. As a result, the clusters must be at the upwind edge of the area that we want to "
  f"protect. At a wind of 2 m/s, the cluster of 100 drones gave only {f1(K1b['zone_last'])} %.")
A(f"We also tested a rhythm in which the drones alternate between two sets of sites. In calm air, a slow rhythm "
  f"(5 minutes at each set) used the energy of 50 drones. It gave a mean decrease of {f0(pR6['mean_30min_pct'])} % "
  f"over 30 minutes. 100 drones in continuous operation gave {f0(pF['mean_30min_pct'])} %, and 50 drones at fixed sites gave "
  f"{f0(pC['mean_30min_pct'])} %. A fast rhythm (30 s) gave almost no effect.")
A("These are simulation results. The pilot must measure the wind at drone height before each operation. In calm "
  "air, a spread fleet with a slow rhythm is correct. In light wind, the drones must fly in clusters at the upwind "
  "edge of the target area. In a wind of 2 m/s or more, inversion mixing gives little benefit.")

# ---------------------------------------------------------------- introduction
A("<a href='#introduction' id='introduction'>Introduction</a>")
A("A. Background")
A(f"Inversion mixing is one of the three pollution control techniques of the Drone V AQI pilot proposal (report 1). "
  f"On calm winter nights, a temperature inversion keeps PM near the ground in Delhi {c('guttikunda12', 'mahrt')}. "
  f"Heavy-lift drones hover in the inversion and blow air down. Each jet brings cleaner air from aloft to street "
  f"level and mixes the polluted layer.")
A(f"Report 1 calculated that 100 heavy drones at 60 m decrease street-level PM by {f1(resolved)} % in 30 minutes. "
  f"This value comes from a unit cell: one drone in a cell with symmetry walls. A unit cell represents the centre "
  f"of a very large fleet in calm air. It has no wind and no zone edges.")
A("A real zone has edges, and the air moves. Wind moves the mixed air out of the zone and brings unmixed air into "
  "the zone. Wind can also bend the jets. Report 1 estimated that a wind of 1 m/s decreases the effect by "
  "approximately one half. This report replaces that estimate with a calculation.")
A("B. Questions of this study")
A("1. Zone edges: how much of the unit-cell effect remains in a zone of 500 m × 500 m with wind?")
A("2. Downwind benefit: how far does the benefit go downwind of the zone?")
A("3. Placement: where must the drones be, relative to the wind and to the area that we want to protect?")
A("4. Rhythm: do adjacent drones that alternate slowly give more mixing for the same energy?")
A("C. A drone fleet is similar to a wind farm")
A(f"A large group of rotors in the lower atmosphere is a known problem in wind energy. In a wind farm, the first rows "
  f"get the undisturbed wind. The rows behind them operate in the wakes of the first rows {c('stevens', 'porteagel20')}. "
  f"Studies of very large wind farms use a periodic unit cell for the inner part of the farm {c('calaf')}. Stable air "
  f"changes the flow through and above a farm {c('allaerts')}. A drone fleet has the same structure, but the rotors "
  f"blow down, not along the wind. Frost-protection wind machines in orchards and vineyards use the same principle as "
  f"inversion mixing {c('ribeiro', 'battany')}.")
A(f"A drone jet that blows down into stable air is a jet in a stable cross-wind. Such jets bend over and stop at a "
  f"depth that depends on the jet momentum, the wind speed and the stability {c('briggs72', 'wright')}. This physics "
  f"controls the results of this report.")

# ---------------------------------------------------------------- methods
A("<a href='#methods' id='methods'>Methods</a>")
A("A. Flow model")
A(f"We used the open-source CFD code OpenFOAM (version v1912) and the best-estimate model of report 1. The model "
  f"solves the time-dependent Reynolds-averaged flow equations (URANS) with the buoyant k-ε turbulence model "
  f"{c('launder')}. The Boussinesq method gives the buoyancy force. In the buoyant k-ε model, stable air decreases "
  f"the turbulence. Its buoyancy coefficient is 1.5. With this value, the model agrees with the field relations for "
  f"wind and temperature profiles in stable air {c('businger', 'dyer', 'hogstrom')}.")
A("The inversion is the same as in report 1. The potential temperature increases by 5 K from the ground to 100 m. "
  "Above 100 m, the air is neutral. PM moves with the air and mixes like heat. As a result, we calculated PM from the "
  "potential temperature. At the start, PM is 100 % at street level and decreases to 50 % at 100 m and above.")
A("B. Two changes to the turbulence model")
A(f"In the first fleet runs, the buoyant k-ε model of OpenFOAM became unstable after 2 to 15 minutes. We found two "
  f"causes. First, the model has a buoyancy term in the equation for k and in the equation for ε. In the ε equation, "
  f"a factor makes the term zero in horizontal flow. In unstable air with horizontal flow, k then increases without "
  f"limit. This occurred in the wakes of the bent jets, where warm air from aloft stays below cooler air. We used the "
  f"value 1 for this factor in unstable air, as Rodi recommends {c('rodi')}. In stable air, the model did not change. "
  f"As a result, the calibration of report 1 is still correct.")
A(f"Second, where k becomes almost zero, the equation for ε becomes very stiff. This occurred above a drone disk "
  f"when the disk started. We set lower limits of 1e-6 m²/s² for k and 1e-8 m²/s³ for ε. At these limits, the eddy "
  f"viscosity is 9e-6 m²/s. This is less than the viscosity of the air, so the limits do not mix the inversion. In "
  f"the periodic cell, the changed model gave {f1(cell_r)} % and the original model gave {f1(cell)} % after 30 "
  f"minutes. The model code is in the repository.")
A("C. Domain and mesh")
A("The fleet zone is 500 m × 500 m. The zone and the fleet are symmetric about the wind axis. As a result, we "
  "simulated one half of the domain, with a plane of symmetry on the wind axis (Figure 1). The full domain is "
  "1,300 m long, 806 m wide and 200 m high. The zone starts 297 m downwind of the inlet. The outlet is 503 m "
  "downwind of the zone. The sides of the domain are 153 m from the sides of the zone.")
A(f"In the zone and near it, the cells are 6.25 m × 6.25 m. Near the inlet, the outlet and the sides, the cells are "
  f"12.5 m long. The cell height is 1.37 m at the ground, 2.75 m at the drones and 15 m at the top. The mesh has "
  f"{RI['cells']:,} cells. The lowest cell layer (0 m to 1.37 m) gives the street-level PM.")
A("D. Wind and boundary conditions")
A(f"At the inlet, the wind speed, k and ε have the profiles of Richards and Hoxey for a neutral surface layer "
  f"{c('richards', 'hargreaves')}. The roughness length is 0.1 m. The wind speed at 10 m (U10) is 1 m/s in the main "
  f"cases and 2 m/s in two tests. At the drone height, the wind speed was {f1(U60)} m/s and {f1(U60b)} m/s. The inlet "
  f"temperature has the inversion profile. The ground is a rough wall. Its wall functions agree with the inlet "
  f"profiles {c('blockenabl')}. The outlet has the hydrostatic pressure of the inversion. The top is a slip wall, and "
  f"the sides are planes of symmetry {c('franke')}.")
A("First, the flow operated for 30 minutes without drones. During this time, the wind adjusted to the stable air. "
  "Then the drones started and operated for 30 minutes. A reference case continued for the same 30 minutes without "
  "drones. All results are relative to the reference case at the same time.")
A("E. Drone model")
A(f"Each drone is an actuator disk at 60 m. The disk applies a downward force equal to the weight of the drone "
  f"(50 kg, 490 N) {c('porteagel11')}. The heavy drone has 5.9 m² of rotor area. The fleet mesh cannot resolve this "
  f"area. As a result, each disk is one column of cells: 6.25 m × 6.25 m (39 m²) and 5.5 m high. A cluster is a "
  f"block of cell columns with the total force of its drones.")
A(f"A larger disk with the same force makes a slower jet. In calm air, a slower jet goes less deep into stable air "
  f"{c('mtt', 'bloomfield')}. We compared four disk sizes in a periodic cell without wind (Results, section A).")
A("F. Cases")
A(table([["Case", "Drones", "Layout", "Schedule", "Wind U10"],
         ["Reference", "0", "None", "No drones", "1 m/s, 2 m/s"],
         ["F100", "100", "Grid of 50 m, sites A and B", "Continuous", "1 m/s, 2 m/s"],
         ["C50", "50", "Sites A (checkerboard)", "Continuous", "1 m/s"],
         ["R600", "50", "Sites A, then sites B", "5 min at A, 5 min at B", "1 m/s"],
         ["R60", "50", "Sites A, then sites B", "30 s at A, 30 s at B", "1 m/s"],
         ["K4", "100", "4 clusters of 25 at the upwind edge, 125 m apart", "Continuous", "1 m/s"],
         ["K1", "100", "1 cluster of 100 at the upwind edge", "Continuous", "1 m/s, 2 m/s"],
         ["Calm-air cells", "100 or 50 per 25 ha", "Periodic 100 m × 100 m cell, sites A and B",
          "As F100, C50, R600 and R60", "0"]]))
A("Group A and group B make a checkerboard pattern (Figure 1). In the rhythm cases, 50 drones move between the A "
  "sites and the B sites. Two groups of 50 drones that take turns give the same flow. A drone that hovers cannot stop "
  "its jet, so a change of site is the practical form of a rhythm. In a cluster, the drones are 6 m to 7 m apart. "
  "This is possible with satellite navigation that has centimetre accuracy.")
A("We ran C50, R600 and R60 with the original turbulence model for 15 minutes. We stopped them because street-level "
  "PM did not change. F100 at 2 m/s also used the original model for 15 minutes. All other cases used the changed "
  "model for 30 minutes.")
B.extend(fig(1, "aqifleet_fig1_layout.png",
             "Fleet layouts and domain. (a) Drone sites of group A (filled) and group B (open). (b) Four clusters. "
             "(c) One cluster. We simulated the half above the plane of symmetry. (d) Side view with the inversion and "
             "the inlet wind."))

# ---------------------------------------------------------------- results
A("<a href='#results' id='results'>Results</a>")
A("A. Checks of the drone model")
A(f"Figure 2 compares four disk sizes in a periodic cell without wind with the resolved disk of report 1. A disk of "
  f"156 m² delayed the effect by 10 to 15 minutes: at 15 minutes, the decrease was only {f0(big5)} %. The one-column "
  f"disk (39 m²) started at the same time as the resolved disk. After 30 minutes, it gave {f1(cell)} % against "
  f"{f1(resolved)} %. As a result, the fleet mesh calculates approximately {f0(100 * (cell / resolved - 1))} % more "
  f"decrease than the resolved disk. A smaller time step changed the 30-minute result by less than 0.2 percentage "
  f"points. Without drones, the difference of street-level PM between the changed and the original turbulence "
  f"models was less than 0.001 %.")
B.extend(fig(2, "aqifleet_fig2_disk_check.png",
             "Check of the drone disk in a periodic cell without wind. Street-level PM as a percentage of the start "
             "value, for four disk sizes and for the resolved disk of report 1."))

A("B. Rhythm in calm air")
A(f"Figure 3 shows the rhythm cases in calm air. 100 drones in continuous operation decreased street-level PM by "
  f"{f0(pF['decrease_pct'][2])} % in 10 minutes and {f1(pF['decrease_30min_pct'])} % in 30 minutes. 50 drones at the A "
  f"sites gave almost the same value after 30 minutes ({f1(pC['decrease_30min_pct'])} %). But the effect started only "
  f"after 20 minutes, because the drones are farther apart.")
A(f"With the slow rhythm (R600), the same 50 drones gave {f0(pR6['decrease_pct'][2])} % after 10 minutes and "
  f"{f1(pR6['decrease_30min_pct'])} % after 30 minutes. The mean decrease over 30 minutes was "
  f"{f1(pR6['mean_30min_pct'])} %. This is {f0(100 * pR6['mean_30min_pct'] / pF['mean_30min_pct'])} % of the value for 100 "
  f"drones ({f1(pF['mean_30min_pct'])} %). It is {pR6['mean_30min_pct'] / pC['mean_30min_pct']:.1f} times the value for "
  f"50 drones at fixed sites ({f1(pC['mean_30min_pct'])} %).")
A(f"The fast rhythm (R60) gave only {f1(pR60['decrease_30min_pct'])} % after 30 minutes. A jet needs a few minutes to "
  f"reach the ground. When the drones move every 30 s, the jets stop before they reach the ground. As a result, a "
  f"useful rhythm has a period of several minutes. A pulsed jet can mix faster than a steady jet. This occurs "
  f"when the pulse period is near the period of the large eddies of the jet {c('choutapalli')}. For a drone jet, this period is "
  f"approximately 1 s. URANS cannot show this effect, and drones cannot change their thrust this fast.")
B.extend(fig(3, "aqifleet_fig3_rhythm.png",
             "Rhythm in calm air. Decrease of street-level PM in a periodic cell with the checkerboard of sites A and B. "
             "Mean: mean decrease over the 30 minutes."))

A("C. Drones on a 50 m grid in wind")
dmax = max([abs(C[k]['zone_last']) for k in ('F100', 'C50', 'R60', 'R600')] + [abs(F100b['zone_last'])])
lmax = max([C[k]['max_ground'] for k in ('F100', 'C50', 'R60', 'R600')] + [F100b['max_ground']])
A(f"In all cases with drones on a 50 m grid, street-level PM did not change. At a wind of 1 m/s and at 2 m/s, the "
  f"zone-mean decrease was less than {max(0.1, round(dmax + 0.05, 1)):.1f} %. The largest local decrease was "
  f"{f1(lmax)} %. The rhythm cases R600 and R60 gave the same result as C50.")
A(f"Figure 4a shows why. In calm air, the jets mix the full layer from the ground to 75 m. In wind, the jets mix "
  f"only the layer from approximately {f0(min(C[k]['mix_bottom_m'] for k in ('F100', 'C50')))} m to 70 m. Figure 4b "
  f"shows the jets in the wind. Each jet bends over and stops approximately {f0(sum(disp_pen) / len(disp_pen))} m "
  f"below the drone. Below approximately {f0(zmix_min)} m, the air does not change.")
B.extend(fig(4, "aqifleet_fig4_dispersed.png",
             "Drones on a 50 m grid. (a) Zone-mean change of PM with height after 30 minutes (15 minutes for C50 and "
             "for the wind of 2 m/s). (b) Change of PM and flow in a vertical section through one row of drones "
             "(squares), wind 1 m/s, 30 minutes."))
A(f"The Briggs equation for the rise of a momentum jet in a stable cross-wind is Δh = 1.5 (Fm / (U N))^(1/3) "
  f"{c('briggs72', 'hanna')}. Fm is the momentum flux of the jet, and U is the wind speed. N is the "
  f"buoyancy frequency of the inversion ({NBV:.3f} 1/s). For one drone, the equation gives {f0(b1)} m at {f1(U60)} m/s and {f0(b1b)} m at "
  f"{f1(U60b)} m/s. This is the depth of the jet centre. The lower edge of the jet is deeper, so the CFD depth of "
  f"approximately {f0(sum(disp_pen) / len(disp_pen))} m agrees with the equation (Figure 5).")
A(f"The equation also shows the solution. The depth increases with the cube root of the momentum flux. When N drones "
  f"fly together, their jets join, and the depth increases by N^(1/3). For a cluster of 25 drones, the equation gives "
  f"{f0(b25)} m. For a cluster of 100 drones, it gives {f0(b100)} m at {f1(U60)} m/s and {f0(b100b)} m at "
  f"{f1(U60b)} m/s. The drones are at 60 m, so a cluster of 100 drones can reach street level in light wind.")
B.extend(fig(5, "aqifleet_fig5_briggs.png",
             "Depth of mixing below the drones against the Briggs equation for a jet in a stable cross-wind. Open "
             "symbols: the mixing reached the ground."))

A("D. Clusters in wind")
A(f"Figure 6 shows the clusters. The cluster of 100 drones (K1) makes one strong jet. "
  + (f"The jet bends over and reaches the ground. The largest decrease at street level is at x = {f0(xk1)} m. "
     if reach(K1) else f"The jet bends over and mixes the air down to {f0(K1['mix_bottom_m'])} m. ")
  + f"The four clusters of 25 drones (K4) mix the air down to "
  + (f"the ground." if reach(K4) else f"approximately {f0(K4['mix_bottom_m'])} m."))
B.extend(fig(6, "aqifleet_fig6_sections.png",
             "Clusters at a wind of 1 m/s after 30 minutes. Change of PM and flow in a vertical section through a "
             "cluster. Positive values: PM decreased."))
A(f"Figure 7 shows the decrease of street-level PM. With K1, the zone-mean decrease was {f1(K1['zone_last'])} %, and "
  f"the largest local decrease was {f0(K1['max_ground'])} %. The area with a decrease of more than 10 % was "
  f"{f1(K1['area_gt10_ha'])} ha. The area with more than 5 % was {f1(K1['area_gt5_ha'])} ha. With K4, the zone-mean "
  f"decrease was {f1(K4['zone_last'])} %, and the largest local decrease was {f0(K4['max_ground'])} %. The areas with "
  f"more than 10 % and more than 5 % were {f1(K4['area_gt10_ha'])} ha and {f1(K4['area_gt5_ha'])} ha.")
B.extend(fig(7, "aqifleet_fig7_maps.png",
             "Decrease of street-level PM after 30 minutes. The square is the fleet zone. The clusters are at its "
             "upwind edge (left)."))
A(f"Figure 8a shows the effect along the wind. The benefit is small at the upwind edge and largest at "
  f"x = {f0(xk1)} m for K1 and x = {f0(xk4)} m for K4. The decrease is more than 5 % up to x = {f0(Lk1)} m for K1 and "
  f"x = {f0(Lk4)} m for K4. After the jet reaches the ground, the warm jet air rises again, because it is lighter "
  f"than the cold air near the ground (Figure 6). As a result, the strong benefit is near the cluster. A weak "
  f"benefit continues downwind of the zone: in the 500 m downwind of the zone, the mean decrease "
  f"was {f1(K1['downwind_500m'])} % for K1 and {f1(K4['downwind_500m'])} % for K4. The benefit outside the zone was "
  f"{f0(out_share(K1))} % of the total benefit for K1 and {f0(out_share(K4))} % for K4.")
A(f"Figure 8c shows the effect across the wind. The benefit of K1 is largest on the wind axis and decreases to the "
  f"sides. In the zone, the central strip (within 125 m of the axis) had {f1(K1['centre_strip'])} % and the outer "
  f"strips had {f1(K1['edge_strip'])} %. The four clusters of K4 cover the full width of the zone "
  f"({f1(K4['centre_strip'])} % and {f1(K4['edge_strip'])} %). As a result, K4 gave approximately twice the zone-mean "
  f"benefit of K1.")
A(f"Figure 8b shows the best position of a 500 m × 500 m target area. For K1, the best target area "
  f"{place(K1)[0]}, with a mean decrease of {f1(K1['best_window'])} %. For K4, the best target area {place(K4)[0]}, "
  f"with {f1(K4['best_window'])} %. As a result, the clusters must be {place(K1)[1]}.")
if K1b:
    A(f"At a wind of 2 m/s, the jet of K1 mixed the air down to approximately {f0(K1b['mix_bottom_m'])} m. The Briggs "
      f"equation gives a deeper jet ({f0(b100b)} m below the drones). The coarse block of the cluster makes a slower "
      f"jet than real drones (Discussion). "
      f"K1 gave {f1(K1b['zone_last'])} % in the zone and {f1(K1b['downwind_500m'])} % in the 500 m "
      f"downwind of the zone. The largest local decrease was {f0(K1b['max_ground'])} %. The best target area "
      f"{place(K1b)[0]}, with {f1(K1b['best_window'])} %.")
B.extend(fig(8, "aqifleet_fig8_profiles.png",
             "Clusters at the upwind edge of the zone. (a) Decrease of street-level PM along the wind. (b) Mean decrease "
             "in a 500 m × 500 m target area against its position. (c) Decrease across the wind."))
A(f"Figure 9 shows the zone-mean decrease against time. "
  + (f"With K1, the effect started after approximately {next((int(t) for t, v in zip(K1['t_min'], K1['series']) if v > 1.0), 30)} "
     f"minutes. " if K1['series'][-1] > 1.0 else "")
  + f"The zone-mean decrease of K4 changed between {f0(min(K4['series'][1:]))} % and {f0(max(K4['series'][1:]))} %, "
  f"because the flow downwind of the clusters is unsteady. The mean over 30 minutes was {f1(K4['zone_mean'])} % for "
  f"K4 and {f1(K1['zone_mean'])} % for K1. In calm air, the unit cell gave {f1(cell)} % after 30 minutes and "
  f"{f1(pF['mean_30min_pct'])} % as a mean, on the same mesh.")
B.extend(fig(9, "aqifleet_fig9_time.png", "Zone-mean decrease of street-level PM against time."))
A(table([["Case", "Wind U10", "Zone, 30 min", "Mean, 30 min", "Downwind 500 m", "Largest local", "Area > 10 %"],
         ["Unit cell, 100 drones", "Calm", f"{f1(cell)} %", f"{f1(pF['mean_30min_pct'])} %", "-", "-", "Whole zone"],
         ["F100: 100 drones, 50 m grid", "1 m/s", f"{f1(F100['zone_last'])} %", f"{f1(F100['zone_mean'])} %",
          f"{f1(F100['downwind_500m'])} %", f"{f1(F100['max_ground'])} %", f"{f1(F100['area_gt10_ha'])} ha"],
         ["K4: 4 clusters of 25", "1 m/s", f"{f1(K4['zone_last'])} %", f"{f1(K4['zone_mean'])} %",
          f"{f1(K4['downwind_500m'])} %", f"{f1(K4['max_ground'])} %", f"{f1(K4['area_gt10_ha'])} ha"],
         ["K1: 1 cluster of 100", "1 m/s", f"{f1(K1['zone_last'])} %", f"{f1(K1['zone_mean'])} %",
          f"{f1(K1['downwind_500m'])} %", f"{f1(K1['max_ground'])} %", f"{f1(K1['area_gt10_ha'])} ha"]]
        + ([["K1: 1 cluster of 100", "2 m/s", f"{f1(K1b['zone_last'])} %", f"{f1(K1b['zone_mean'])} %",
             f"{f1(K1b['downwind_500m'])} %", f"{f1(K1b['max_ground'])} %", f"{f1(K1b['area_gt10_ha'])} ha"]] if K1b else [])))

# ---------------------------------------------------------------- discussion
A("<a href='#discussion' id='discussion'>Discussion</a>")
A("A. Answers to the four questions")
A(f"1. Zone edges: in wind, drones on a 50 m grid have no effect at street level. As a result, their edges do "
  f"not change the result. With a cluster, the benefit has the shape of a plume that starts at the cluster and goes "
  f"downwind. The sides of the zone between the plumes get little benefit, and four clusters cover the width of the "
  f"zone better than one cluster. We did not simulate a zone with edges in calm air.")
A(f"2. Downwind benefit: with a cluster, a weak benefit of approximately 2 % continues for at least 500 m downwind "
  f"of the zone. After 30 minutes, {f0(out_share(K1))} % (K1) and {f0(out_share(K4))} % (K4) of the total benefit was "
  f"outside the zone.")
A(f"3. Placement: measure the wind at the drone height before each operation. In calm air, use the spread fleet of "
  f"report 1. In light wind, fly the drones in clusters across the wind. Put the clusters {place(K1)[1]}. In our "
  f"cases, four clusters of 25 drones gave twice the zone-mean benefit of one cluster of 100 drones. The four "
  f"clusters covered the full width of the zone. A large cluster reaches the ground in a stronger wind than a small cluster. "
  f"In a wind of 2 m/s or more, even a cluster of 100 drones gives little benefit.")
A(f"4. Rhythm: in calm air, a slow rhythm of approximately 5 minutes gives almost the effect of 100 drones. It "
  f"uses the energy of 50 drones ({f0(50 * POWER_KW)} kW instead of {f0(100 * POWER_KW)} kW). A fast rhythm of 30 s does not "
  f"work. In wind, the rhythm does not help drones on a 50 m grid, because their jets do not reach the ground.")
A("B. Limits of the simulations")
A(f"1. Turbulence model: URANS with the k-ε model calculates the mean flow, not the individual eddies. A study of "
  f"jets in a stable cross-wind with large-eddy simulation (LES) is the next step {c('calaf', 'porteagel11')}. We "
  f"changed the buoyancy term of the model in unstable air. This change is standard {c('rodi')}, but it adds "
  f"uncertainty to the wakes of the jets.")
A(f"2. Drone disk: the fleet mesh calculates approximately {f0(100 * (cell / resolved - 1))} % more decrease than the "
  f"resolved disk in calm air. In wind, the depth of a bent jet depends on its momentum flux, which the coarse disk "
  f"keeps. The real jets of a cluster are faster than the coarse block, so the real clusters can go deeper.")
A("3. Site: the ground is flat and has a uniform roughness. There are no buildings. The ground does not cool the air "
  "during the operation, and there are no new PM emissions. In a real night, the ground continues to cool the air "
  "near it, and traffic and other sources continue to add PM. As a result, the benefit far downwind can be smaller.")
A("4. Wind: the wind is steady and has one direction, along the rows of drones. Real night winds change speed and "
  "direction. A change of the wind direction moves the area of benefit. The cluster must then move.")
A("5. Domain: we simulated one half of the zone, so the flow is symmetric about the wind axis. The 30-minute "
  "operation shows the start of the operation, not a full night.")
A("6. Cluster flight: in a cluster, the drones are 6 m to 7 m apart and fly in the joined jet of their neighbours. "
  "The pilot must test the safe distance between the drones and the effect of the joined jet on their thrust.")
A("C. Conclusion")
A(f"Wind is the most important condition for inversion mixing. In calm air, a spread fleet decreases street-level PM "
  f"by approximately {f0(cell)} % in 30 minutes, and a slow rhythm saves half of the energy. In a light wind of "
  f"1 m/s, the jets of single drones bend over and do not reach street level. A cluster of drones makes one strong "
  f"jet that can reach street level."
  + (f" It decreased PM by up to {f0(K1['max_ground'])} % in a plume downwind of the cluster." if K1['max_ground'] > 2 else "")
  + f" Four clusters of 25 drones decreased the zone-mean PM by {f0(K4['zone_last'])} % after 30 minutes. At 2 m/s, the "
  f"benefit was small. Report 1 estimated that a wind of 1 m/s halves the effect. For a spread fleet, the effect is "
  f"zero, and clusters are necessary. The pilot must test these results with measurements of the wind and of PM at "
  f"different heights.")

# ---------------------------------------------------------------- references
A("<a href='#references' id='references'>References</a>")
for i, r in enumerate(refs, 1):
    A(f"[<a href='#ref{i}' id='ref{i}'>{i}</a>]: {r['citation']} <a href='https://doi.org/{r['doi']}'>Link</a>")

# ---------------------------------------------------------------- appendix
A("<a href='#appendix' id='appendix'>Appendix</a>")
A(table([["Parameter", "Value"],
         ["Inversion", "Potential temperature +5 K from 0 m to 100 m; neutral above 100 m; N = %.3f 1/s" % NBV],
         ["PM at the start", "100 % at street level; 50 % at 100 m and above; PM calculated from the potential temperature"],
         ["Wind at the inlet", f"Richards and Hoxey profiles; U10 = 1 m/s and 2 m/s; roughness length 0.1 m; wind at 60 m: "
                               f"{f1(U60)} m/s and {f1(U60b)} m/s"],
         ["Turbulence model", "Buoyant k-ε, Cg = 1.5, turbulent Prandtl and Schmidt numbers 0.7; C3 = 1 in unstable air; "
                              "k ≥ 1e-6 m²/s², ε ≥ 1e-8 m²/s³"],
         ["Heavy drone", f"50 kg, force 490 N, hover at 60 m; hover power {POWER_KW:.1f} kW (report 1)"],
         ["Drone disk in the fleet mesh", "One cell column, 6.25 m × 6.25 m × 5.5 m (39 m²)"],
         ["Clusters", "K4: 4 blocks of 31 m × 31 m (25 drones each); K1: 1 block of 63 m × 69 m (100 drones)"],
         ["Fleet", "100 sites on a 50 m grid over 500 m × 500 m; groups A and B in a checkerboard pattern"],
         ["Domain", "1,300 m × 806 m × 200 m; one half simulated (plane of symmetry on the wind axis)"],
         ["Mesh", f"{RI['cells']:,} cells; 6.25 m horizontal cells in and near the zone; 1.37 m cell height at the ground"],
         ["Time", "30 min spin-up without drones, then 30 min of operation; time step set by a Courant number of 1.5 or less"],
         ["Computer time", f"{cpu:.0f} core-hours for the fleet cases (4-core cloud computer)"]]))
A("All simulation code and data are in the Drone V AQI repository (github.com/Tech-Studio-Art/dronevaqi), directory "
  "plans/sim. The fleet model is plans/sim/cfd/fleet3d.py, and the changed turbulence model is in "
  "plans/sim/cfd/turbulence. This text follows the rules of ASD-STE100 Simplified Technical English.")

entry = {"title": TITLE, "body": B}
json.dump(entry, open(os.path.join(HERE, "data_entry.json"), "w"), indent=1, ensure_ascii=False)

# ---------------------------------------------------------------- print HTML + PDF
css = """
@page { size: A4; margin: 16mm 15mm 16mm 15mm; }
body { font-family: 'DejaVu Sans', Arial, sans-serif; font-size: 9.6pt; line-height: 1.42; color: #111; }
p { margin: 0 0 5pt 0; }
p.h { font-weight: bold; margin-top: 8pt; }
a { color: #1c5cab; text-decoration: none; }
p.sec { font-size: 13pt; font-weight: bold; border-bottom: 1px solid #bbb; margin-top: 12pt; padding-bottom: 2pt; }
img { display: block; max-width: 100%; max-height: 82mm; margin: 6pt auto 2pt auto; }
p.cap { font-size: 8.4pt; color: #333; margin-bottom: 8pt; }
table { border-collapse: collapse; width: 100%; margin: 4pt 0 8pt 0; font-size: 8.2pt; page-break-inside: avoid; }
td { border: 1px solid #c8c8c8; padding: 3pt 4pt; vertical-align: top; }
td p { margin: 0; }
p.ref { font-size: 8.2pt; margin-bottom: 3pt; }
"""
html = [f"<html><head><meta charset='utf-8'><title>Drone V AQI report 2</title><style>{css}</style></head><body>"]
SEC = {"abstract", "introduction", "methods", "results", "discussion", "references", "appendix"}
for b in B:
    if b.startswith("<img") or b.startswith("<table"):
        html.append(b)
    elif any(f"id='{s}'" in b for s in SEC):
        html.append(f"<p class='sec'>{b}</p>")
    elif b.startswith("<a href='#fig"):
        html.append(f"<p class='cap'>{b}</p>")
    elif b.startswith("[<a href='#ref"):
        html.append(f"<p class='ref'>{b}</p>")
    elif len(b) < 70 and b[:2] in ("A.", "B.", "C.", "D.", "E.", "F.", "G."):
        html.append(f"<p class='h'>{b}</p>")
    else:
        html.append(f"<p>{b}</p>")
html.append("</body></html>")
open(os.path.join(HERE, "report.html"), "w").write("\n".join(html))
subprocess.run(["node", os.path.join(HERE, "pdf.js"), os.path.join(HERE, "report.html"),
                os.path.join(HERE, "fleet_report.pdf")], check=True)
print("entries", len(B))
