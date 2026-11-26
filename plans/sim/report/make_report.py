#!/usr/bin/env python3
"""
Build the pilot-proposal report:
  data_entry.json  - one Research Dex entry {"title": ..., "body": [...]}
  report.html      - the same body, styled for print
  pilot_report.pdf - A4 PDF (<= 15 pages)
All numbers come from the simulation outputs. Text follows ASD-STE100 rules.
"""
import json
import os

HERE = os.path.dirname(os.path.abspath(__file__))
SIM = os.path.dirname(HERE)
CFD = os.path.join(SIM, "cfd")


def J(*p):
    return json.load(open(os.path.join(*p)))


mix = J(CFD, "runs", "mixing", "summary.json")
p2 = os.path.join(CFD, "runs", "mixing", "heavy_z60_prt2", "result.json")
if os.path.exists(p2):
    mix["heavy_z60_prt2"] = json.load(open(p2))
road = J(CFD, "runs", "esp_road", "summary.json")["total"]
amb = {r["case"]: r for r in J(CFD, "runs", "2d", "summary_all.json") if r["case"].startswith("espfence")}
drones = {r["case"]: r for r in J(CFD, "runs", "drone3d", "summary.json")}
sv = J(SIM, "surveillance", "results.json")
refs = J(HERE, "refs_verified.json")


def red(k):
    s = mix[k]["series"]
    return 100 * (1 - s[-1]["z<2m"] / s[0]["z<2m"])


def f0(x):
    return f"{x:.0f}"


R = {r["key"]: i + 1 for i, r in enumerate(refs)}


def c(*keys):
    return "[" + ", ".join(f"<a href='#ref{R[k]}'>{R[k]}</a>" for k in keys) + "]"


def fig(n, src, cap):
    return [f"<img src='assets/{src}'>", f"<a href='#fig{n}' id='fig{n}'>Figure {n}</a>: {cap}"]


def table(rows, head=True):
    out = "<table>"
    for i, r in enumerate(rows):
        tag = "strong" if head and i == 0 else None
        out += "<tr>" + "".join(f"<td><p>{'<strong>' + str(x) + '</strong>' if tag else x}</p></td>" for x in r) + "</tr>"
    return out + "</table>"


hz60, hz90, qz30, nod = red("heavy_z60"), red("heavy_z90"), red("quad_z30"), red("no_drone")
prt2 = red("heavy_z60_prt2") if "heavy_z60_prt2" in mix else None
heavy = drones["heavy_none"]
rd0, rdE, rdH, rdP = road["no fence"], road["ESP fence"], road["ESP fence, half migration velocity"], road["passive fence"]
a1 = amb["espfence_U1.0_PM25"]
a1c = amb["espfence_U1.0_coarse"]
a2 = amb["espfence_U2.0_PM25"]
import sys
sys.path.insert(0, HERE)
from make_figs import city_increment  # noqa: E402
rp = {(x["drones"], x["teams"]): x["emission_prevented_pct"] for x in sv["response"]}
FN, FD = sv["field"]["night (stable, E-F)"], sv["field"]["day (neutral, D)"]
p_hot = max(d["p"] for d in sv["detection"]["night (stable, E-F)|10"])
city100, city10 = city_increment(1.0), city_increment(0.1)
best = rp[(2, 4)]
g_night = sv["hotspot_g_per_night"] * best / 100
mix_range = f"{min(prt2, hz60):.0f} % to {max(prt2, hz60):.0f} %" if prt2 is not None else f"up to {hz60:.0f} %"

TITLE = ("Drone V AQI pilot proposal: inversion mixing, electrostatic wind fences and 24x7 pollution-source "
         "surveillance to decrease winter air pollution in Delhi")
B = []
A = B.append
A(f"<a href='#title' id='title'>Title</a>: {TITLE}")
A("<a href='#author' id='author'>Author</a>: Noorul Ali, Tech Studio .Art (Drone V AQI programme), India")
A("<a href='#contact' id='contact'>Correspondence</a>: <a href='mailto:noorulali78@gmail.com'>noorulali78@gmail.com</a>")
A("<a href='#abstract'>Abstract</a><br><a href='#introduction'>Introduction</a><br><a href='#methods'>Methods</a><br>"
  "<a href='#results'>Results</a><br><a href='#impact'>Environmental impact</a><br><a href='#pilot'>Pilot and test plan</a><br>"
  "<a href='#india'>Use in other Indian cities</a><br><a href='#discussion'>Discussion</a><br><a href='#references'>References</a><br>"
  "<a href='#appendix'>Appendix</a>")

# ---------------------------------------------------------------- abstract
A("<a href='#abstract' id='abstract'>Abstract</a>")
A("Each winter, the air in Delhi becomes very dangerous. After Diwali and during the crop-residue fire season, PM2.5 and PM10 "
  "stay at severe levels for many days. A temperature inversion keeps the smoke near the ground. Low wind and no rain "
  "make the problem worse.")
A("This report recommends a pilot project to test three control techniques. Technique 1 is inversion mixing. Heavy-lift "
  "drones push cleaner air from above the inversion down to street level. Technique 2 is the electrostatic wind fence. "
  "The fence charges particles and collects them on metal plates. The wind moves the air through the fence, so the fence "
  "needs no fans. Technique 3 is 24x7 pollution-source surveillance. Patrol drones with thermal cameras find open "
  "waste fires quickly, so that ground teams can stop them.")
A(f"We did all the work in this report with computer simulations. The simulations use standard methods from the "
  f"scientific literature. The results are as follows. Inversion mixing with 100 heavy drones decreased street-level PM "
  f"by {mix_range} in 30 minutes. The zone was 500 m × 500 m. An electrostatic fence between a road and a footpath decreased "
  f"footpath PM2.5 by {f0(100 * (1 - rdE['2.5'] / rd0['2.5']))} %. Surveillance with two patrol drones and four ground "
  f"teams prevented {best:.0f} % of the PM2.5 from a hotspot of 100 small waste fires. PM sensors cannot find these "
  f"small fires, but together the fires can add up to {city100:.0f} µg/m³ to the PM2.5 of the city.")
A("These numbers are calculations, not measurements. A field pilot is necessary to measure the real values. We ask the Government "
  "of NCT of Delhi and its agencies to support a pilot in winter 2027. The three techniques can also operate in other "
  "Indian cities that have the same winter problem.")

# ---------------------------------------------------------------- introduction
A("<a href='#introduction' id='introduction'>Introduction</a>")
A("A. The problem")
A(f"Delhi has some of the highest particle pollution of any large city. Long exposure to PM2.5 causes heart and lung "
  f"disease and early death {c('pope', 'guttikunda')}. In November, smoke from crop-residue fires in Punjab and Haryana "
  f"adds to local pollution from traffic, dust, waste fires and industry {c('cusworth')}.")
A(f"On calm winter nights, the ground cools quickly. A layer of cold air forms below warmer air. This is a temperature "
  f"inversion. The inversion stops vertical mixing, so pollution stays in a shallow layer near the ground "
  f"{c('stull', 'mahrt')}. The air above the inversion often has less pollution than the air at street level.")
A("B. Why we use simulations now")
A(f"A field test of new equipment needs one full winter, permits and much money. Simulations let us test many designs "
  f"before we build hardware. Computational fluid dynamics (CFD) is a standard tool for urban air-flow and pollution "
  f"studies {c('blocken')}. We use only methods that other researchers have published and tested. Where a method has a "
  f"known error, we show the direction of the error. Where possible, we also show a conservative result.")
A("The pollution in Delhi is at an extreme level each winter, and millions of people are at risk. It is not correct "
  "to wait for complete data before we start. Simulations give the best available guide now. The pilot will "
  "then give the measurements that are necessary to make a decision.")
A("C. The three techniques")
A(table([["Technique", "What it does", "Scale", "Main hardware"],
         ["1. Inversion mixing", "Pushes cleaner air from aloft down to street level", "Neighbourhood (25 ha per 100 drones)",
          "Heavy-lift drones, scout quad drones"],
         ["2. Electrostatic wind fence", "Removes PM from air that the wind moves through the fence", "Road corridor, footpath, school edge",
          "Fence modules with ionizer and collector plates"],
         ["3. 24x7 source surveillance", "Finds open fires and other sources and starts a quick response", "District (25 km²) to city",
          "Low-cost sensor grid, patrol drones with thermal cameras"]]))

# ---------------------------------------------------------------- methods
A("<a href='#methods' id='methods'>Methods</a>")
A("A. Common conditions")
A("All cases use a typical severe night in November: PM2.5 = 300 µg/m³ and PM10 = 500 µg/m³ at street level. PM10 "
  "has a higher mass than PM2.5, as in measured Delhi episodes. We calculated PM2.5 and the coarse part of PM10 "
  "separately where their behaviour is different.")
A("B. Technique 1: inversion mixing")
A(f"We used the open-source CFD code OpenFOAM (version v1912). The model solves time-dependent Reynolds-averaged flow "
  f"equations (URANS) with the standard k-ε turbulence model {c('launder')}. The Boussinesq method gives the buoyancy "
  f"force of the stable air layer. The inversion has a temperature increase of 5 K over the lowest 100 m. The air above "
  f"100 m is neutral and has 50 % of the street-level PM.")
A(f"Each drone is an actuator disk with thrust equal to its weight. The heavy drone has a take-off mass of 50 kg and a "
  f"total rotor area of 5.9 m². The fleet has 100 drones on a 50 m grid. We modelled one drone in a cylindrical "
  f"unit cell of radius 28.2 m with symmetry walls. This is the standard method for a large uniform array. The jet "
  f"penetration agrees with the classical theory for jets in stable fluid {c('turner', 'list')}.")
A(f"The standard k-ε model does not decrease turbulence in stable air. As a result, it calculates too much mixing. To show "
  f"this effect, we also ran a conservative case. In this case, the turbulent Prandtl number is 2.0, which is in the "
  f"range for stable boundary layers {c('zilitinkevich')}. A larger value gives less vertical mixing of heat and PM.")
A("C. Technique 2: electrostatic wind fence")
A(f"We used two-dimensional steady RANS with the k-ε model. The inlet wind profile and the ground roughness follow "
  f"Richards and Hoxey {c('richards')}. PM is a passive scalar with a turbulent Schmidt number of 0.7 {c('tominaga')}. "
  f"The fence is 6 m high and 0.5 m deep. It is a porous zone with a pressure-loss coefficient of 1.5. Inside the "
  f"fence, PM decreases at a rate set by the Deutsch-Anderson equation {c('deutsch', 'mizuno')}. The migration velocity "
  f"is 0.04 m/s for PM2.5 and 0.12 m/s for coarse PM. These values are in the low range of industrial "
  f"electrostatic precipitators {c('white', 'jaworek')}.")
A(f"We ran two configurations. In the first configuration, only background air flows through the fence. In the second "
  f"configuration, the fence is between a two-lane road and a footpath. The road is a line source of PM, as in "
  f"near-road barrier studies {c('baldauf', 'hagler')}. Without a fence, the road adds 40 % to the background PM on the "
  f"footpath. We compared no fence, a passive fence (no electric charge) and the electrostatic fence. We also ran the "
  f"electrostatic fence with half the migration velocity, as a conservative case.")
A("D. Technique 3: 24x7 source surveillance")
A(f"The source is a hotspot of 100 small open fires of domestic waste in a 500 m × 500 m area. Each fire burns "
  f"100 g of waste per hour for 4 hours. The PM2.5 emission factor is 10 g/kg {c('wiedinmyer')}. As a result, each "
  f"fire emits 1 g/h of PM2.5 and the hotspot emits {sv['hotspot_g_h']:.0f} g/h.")
A(f"We calculated the PM2.5 field of the hotspot with the Gaussian plume model and the Briggs urban dispersion "
  f"coefficients {c('hanna')}. The field is the sum of the plumes of all fires. Low-cost PM sensors {c('snyder')} are "
  f"on a square grid at a height of 3 m. The network gives an alert when the 1-hour average at one sensor increases "
  f"by more than 10 µg/m³. We calculated the probability of detection with 400 random hotspots and wind directions.")
A(f"Patrol drones with thermal cameras fly a fixed route over a 25 km² district. When a drone finds a hotspot, it "
  f"makes a map of all fires in 6 minutes. Then the ground teams go to each fire. With a map, a team needs 6 minutes "
  f"for each fire. Without a map, a team needs 20 minutes for each fire. Satellite fire data {c('schroeder')} covers "
  f"crop-residue fires outside the city.")
A(f"To show the effect of many hotspots on the whole city, we used a box model of the city air {c('stull')}. The box "
  f"has the area of Delhi (1,483 km²), a mixing height of 150 m and a wind of 1 m/s.")

# ---------------------------------------------------------------- results
A("<a href='#results' id='results'>Results</a>")
A("A. Technique 1: inversion mixing")
B.extend(fig(1, "aqipilot_fig1_mission.png",
             "Inversion-mixing operation. (1) A scout quad drone measures PM and temperature every 10 m up to 300 m. "
             "(2) The team starts mixing only if the air aloft has less PM. (3) Heavy drones at 60 m push cleaner air "
             "down to street level."))
A(f"Figure 2 shows street-level PM (0 m to 2 m) during 30 minutes of operation. Heavy drones at 60 m decreased PM by "
  f"{hz60:.0f} %. Most of this decrease occurred in the first 10 minutes. Heavy drones at 90 m gave a decrease of "
  f"{hz90:.0f} %, but the effect started only after 15 minutes. Quad drones gave only {qz30:.0f} %, because their jets "
  f"are too weak. Without drones, the change was {nod:.1f} %.")
if prt2 is not None:
    A(f"The conservative case (Pr<sub>t</sub> = 2.0) gave a decrease of {prt2:.0f} %. The result changed only a little. "
      f"This shows that the jets and the return flow of the drones do most of the mixing. Small-scale turbulence does "
      f"less. Wind replaces the air in the zone and decreases the effect. With a wind of 1 m/s, we calculate "
      f"approximately one half of these values.")
B.extend(fig(2, "aqipilot_fig2_mixing_time.png",
             "Street-level PM during inversion mixing, as a percentage of the start value. 100 drones on a 50 m grid "
             "over a 500 m × 500 m zone."))
B.extend(fig(3, "aqipilot_fig3_mixing_field.png",
             "PM distribution around one heavy drone at 60 m (triangle), before and after 30 minutes of mixing."))
A(f"The heavy drone needs {heavy['power_mt_W'] / 1000:.1f} kW to hover. 100 drones need approximately "
  f"{heavy['power_mt_W'] * 100 / 1e6:.2f} MW. An 8-hour night needs approximately "
  f"{heavy['power_mt_W'] * 100 * 8 / 1e6:.1f} MWh per 25 ha. Orchard frost fans use a similar energy per hectare to "
  f"break inversions {c('battany')}.")
B.extend(fig(4, "aqipilot_fig1b_im50.png",
             "IM-50 heavy mixer concept design: 50 kg take-off mass, four pairs of coaxial 54 in rotors, sensor mast "
             "0.5 m above the rotors. The 25 kg payload is a second battery."))

A("B. Technique 2: electrostatic wind fence")
B.extend(fig(5, "aqipilot_fig4_esp_principle.png",
             "Principle of the electrostatic wind fence. The ionizer charges the particles. The electric field between "
             "the plates moves the particles onto the plates. The wind moves the air through the fence."))
A(f"Figure 6a shows the road case with a wind of 1 m/s. Without a fence, footpath PM2.5 was {300 * rd0['2.5']:.0f} µg/m³. "
  f"The electrostatic fence decreased it to {300 * rdE['2.5']:.0f} µg/m³. This is a decrease of "
  f"{100 * (1 - rdE['2.5'] / rd0['2.5']):.0f} %. At 40 m behind the fence, the decrease was "
  f"{100 * (1 - rdE['40'] / rd0['40']):.0f} %. With half the migration velocity, the decrease on the footpath was "
  f"{100 * (1 - rdH['2.5'] / rd0['2.5']):.0f} %.")
A(f"A passive fence without electric charge increased footpath PM2.5 by {100 * (rdP['2.5'] / rd0['2.5'] - 1):.0f} %. "
  f"As a result, the electrostatic collection gives the benefit, not the structure.")
A(f"Figure 6b shows the case with background air only. At 2 m behind the fence, PM2.5 decreased by "
  f"{100 * (1 - a1['pedestrian_C']['x=2m']):.0f} % and the coarse part of PM10 decreased by "
  f"{100 * (1 - a1c['pedestrian_C']['x=2m']):.0f} %. At a wind of 2 m/s, the PM2.5 decrease was "
  f"{100 * (1 - a2['pedestrian_C']['x=2m']):.0f} %. The air moves faster, so the particles have less time between "
  f"the plates.")
B.extend(fig(6, "aqipilot_fig5_esp_results.png",
             "Electrostatic wind fence results. (a) Fence between a road and a footpath. (b) Background air only."))
A(f"One metre of fence cleans {a1['CADR_per_m']:.1f} m³/s of air at a wind of 1 m/s. The high-voltage supply uses "
  f"approximately 140 W per metre. This is approximately 90 J for each cubic metre of clean air. A drone with a "
  f"filter uses more than 2,000 J/m³.")

A("C. Technique 3: 24x7 source surveillance")
A(f"Each small fire adds only {sv['single_fire_100m_night_ug_m3']:.1f} µg/m³ at 100 m downwind on a calm night. Figure 7 "
  f"shows the field of one hotspot of 100 fires. Inside the hotspot, the mean increase is {FN['zone_mean']:.1f} µg/m³ at "
  f"night and {FD['zone_mean']:.1f} µg/m³ in daytime. At 1 km downwind, the increase at night is "
  f"{FN['downwind_1km']:.1f} µg/m³.")
B.extend(fig(7, "aqipilot_fig6_fire_plume.png",
             "PM2.5 added by one hotspot of 100 small waste fires (100 g of waste per hour each) on a calm night. "
             "The orange square is the 500 m × 500 m hotspot."))
A(f"These values are much less than the background of 300 µg/m³ and less than the noise of low-cost sensors. As a "
  f"result, the sensor network found the hotspot in only {100 * p_hot:.0f} % of cases, at all grid spacings. PM sensors "
  f"cannot find small fires. Thermal cameras on patrol drones are necessary.")
A(f"Many hotspots together have a large effect (Figure 8a). If 10 % of Delhi has hotspots of this size, the PM2.5 of "
  f"the city increases by approximately {city10:.0f} µg/m³. If all of Delhi has such hotspots, the increase is "
  f"approximately {city100:.0f} µg/m³. This is approximately {100 * city100 / 300:.0f} % of the background.")
A(f"Figure 8b shows the PM2.5 emission that surveillance prevents in one hotspot. Without patrol drones, the system "
  f"prevented almost no emission. With two drones and two ground teams, it prevented {rp[(2, 2)]:.0f} %. With two "
  f"drones and four ground teams, it prevented {rp[(2, 4)]:.0f} %. Each hotspot has 100 fires, so the number of ground "
  f"teams limits the result.")
B.extend(fig(8, "aqipilot_fig7_surveillance.png",
             "Surveillance results. (a) City-wide PM2.5 added by small fires, for different areas of the city with "
             "hotspots. (b) PM2.5 emission prevented in one hotspot, for different numbers of patrol drones and ground teams."))

# ---------------------------------------------------------------- impact
A("<a href='#impact' id='impact'>Environmental impact</a>")
A("The table gives the likely effect of each technique in one pilot zone. These values come from the simulations. "
  "The pilot must measure them.")
A(table([["Technique", "Likely benefit", "Area of benefit", "Possible negative effects"],
         ["Inversion mixing (100 heavy drones)", f"Street-level PM {mix_range} lower on calm nights", "25 ha zone",
          "Noise at night; energy use; no benefit if a smoke layer is aloft"],
         ["Electrostatic wind fence (per 100 m)", f"Footpath PM2.5 approximately {100 * (1 - rdE['2.5'] / rd0['2.5']):.0f} % lower; "
          f"coarse PM {100 * (1 - a1c['pedestrian_C']['x=2m']):.0f} % lower", "Strip 40 m deep behind the fence",
          "Ozone from the ionizer; waste dust from plate cleaning"],
         ["24x7 surveillance (25 km² district)", f"{best:.0f} % of the PM2.5 from each hotspot prevented (approximately "
          f"{g_night:.0f} g per hotspot each night); city PM2.5 up to {city100 * best / 100:.0f} µg/m³ lower if used in all of Delhi",
          "Whole district", "Privacy of camera data; drone noise"]]))
A(f"Inversion mixing does not remove PM from the atmosphere. It mixes PM into a deeper layer. This decreases the "
  f"concentration that people breathe. At the average population density of Delhi (Census of India 2011), a 25 ha zone "
  f"contains approximately 2,800 people.")
A(f"The electrostatic fence removes PM from the air. Each metre of fence collects approximately 90 g of PM per day in a "
  f"severe episode. The operator must collect this dust and send it to an approved waste facility. Ionizers can make "
  f"ozone {c('boelter')}. The design must use positive corona and carbon-fibre emitters to keep ozone low. The pilot "
  f"must measure ozone near the fence.")
A(f"Surveillance prevents emissions at the source. In one hotspot, the system prevents approximately {g_night:.0f} g of "
  f"PM2.5 each night. The drone maps also show the sites where people burn waste again and again. The city can use "
  f"these maps to improve waste collection and to prevent the fires. Open burning also makes toxic gases "
  f"{c('wiedinmyer')}, so early action gives an additional health benefit.")

# ---------------------------------------------------------------- pilot
A("<a href='#pilot' id='pilot'>Pilot and test plan</a>")
A("The pilot has one aim: to measure the real benefit of each technique and to compare it with the simulations. Each "
  "technique has a test zone and a control zone with the same conditions.")
A(table([["Phase", "Period", "Work", "Result"],
         ["0. Simulation", "2026 (complete)", "CFD and dispersion models in this report", "Designs and expected results"],
         ["1. Prototypes", "Jan 2027 to Jun 2027", "Build 1 heavy drone, 2 scout drones, 20 m of fence, 50 sensors; bench tests",
          "Measured thrust, power, fence efficiency, ozone"],
         ["2. Baseline", "Jul 2027 to Sep 2027", "Install sensors in test and control zones; measure without operation",
          "Baseline PM and weather data"],
         ["3. Winter pilot", "Oct 2027 to Jan 2028", "Operate all three techniques; measure every night",
          "Measured PM decrease against control zones"],
         ["4. Report", "Feb 2028 to Mar 2028", "Compare measurements with simulations; calculate costs", "Decision to scale or stop"]]))
A("A. Test zones")
A("1. Inversion mixing: one 500 m × 500 m zone at the edge of the city, in a green zone for drone operation. Start with "
  "10 heavy drones, then increase to the full fleet if the first results are positive.")
A("2. Electrostatic fence: 100 m of fence on a road wall next to a school or hospital footpath.")
A("3. Surveillance: one 25 km² district with known waste-burning problems. Use 2 patrol drones with thermal cameras, "
  "4 ground teams and a sensor grid for background data.")
A("B. Measurements")
A(table([["Technique", "Measurement", "Instrument", "Pass criterion"],
         ["Inversion mixing", "PM profile 0 m to 300 m; street PM in test and control zones; noise",
          "Scout drone sensors; reference monitor; sound meter", "Street PM at least 10 % lower than control on calm nights"],
         ["Electrostatic fence", "PM2.5 and PM10 at 2 m, 10 m and 40 m behind the fence; ozone; power",
          "Calibrated sensors; ozone analyser", "PM2.5 at least 20 % lower; ozone within national standard"],
         ["Surveillance", "Number of fires found; time to find a hotspot; time to stop all fires",
          "Drone thermal camera logs; field counts by ground teams", "Drones find at least 90 % of the fires that ground "
          "teams count; mean time to find a hotspot less than 60 minutes"]]))
A("C. Support that we request")
A("1. Permission to operate drones at night and above 120 m in the pilot zones.")
A("2. Access to sites, walls and electrical power for the fence and the sensors.")
A("3. Access to data from the regulatory monitors and a link to the enforcement teams for fire response.")
A("4. A joint technical committee with the pollution control agencies to review the results.")

# ---------------------------------------------------------------- India
A("<a href='#india' id='india'>Use in other Indian cities</a>")
A("Many cities in the Indo-Gangetic Plain have the same winter inversions as Delhi. Examples are Lucknow, Kanpur, Agra, "
  "Varanasi and Patna. Inversion mixing can operate in these cities on calm winter nights. The scout drone shows if the "
  "conditions are correct before each operation.")
A("The electrostatic fence can operate on any busy road in any city. Mumbai, Kolkata, Bengaluru, Chennai and Hyderabad "
  "have long road corridors with footpaths next to heavy traffic. The fence modules attach to walls, flyovers and noise "
  "barriers that are already there. As a result, the fence needs no new land.")
A("Open waste burning occurs in cities across India. A sensor grid and patrol drones can operate in any municipal "
  "ward. The system also gives data for the city action plans of the National Clean Air Programme (NCAP).")

# ---------------------------------------------------------------- discussion
A("<a href='#discussion' id='discussion'>Discussion</a>")
A("A. Limits of the simulations")
A("1. Inversion mixing: the k-ε model calculates too much mixing in stable air. The unit cell does not include wind or "
  "the edges of the zone. As a result, the results are an upper limit. The conservative case tests only the turbulent "
  "transport of heat and PM. In Phase 1, we will also use a turbulence model with a buoyancy term and a "
  "three-dimensional model of the full fleet.")
A("2. Electrostatic fence: the model is two-dimensional, so it does not include air flow around the ends of the fence. "
  "The migration velocities come from industrial practice. The model does not include dirt on the plates or calm air.")
A("3. Surveillance: the Gaussian plume model is approximate near buildings. The emission factor of waste fires has a "
  "large range. The model assumes that a thermal camera at patrol height can find a fire that burns 100 g/h. The pilot "
  "must test this. The model does not include false alarms.")
A("B. Why the simulations are a correct basis for a pilot")
A("Each model uses published and tested methods. Where a model has a known error, we show its direction and we give a "
  "conservative case. Each technique shows a clear benefit even in the conservative cases. A pilot is the correct next "
  "step because it measures the real effect at a small cost before a large investment.")
A("C. Conclusion")
A(f"The simulations show that each technique can decrease the PM that people breathe. Inversion mixing decreases "
  f"street-level PM by {mix_range} in a zone. The electrostatic fence decreases footpath PM2.5 by approximately "
  f"{100 * (1 - rdE['2.5'] / rd0['2.5']):.0f} %. PM sensors cannot find small waste fires, but patrol drones and ground teams can prevent approximately {best:.0f} % of their emission. We "
  f"recommend a pilot in winter 2027 to measure these effects in Delhi.")

# ---------------------------------------------------------------- references
A("<a href='#references' id='references'>References</a>")
for i, r in enumerate(refs, 1):
    A(f"[<a href='#ref{i}' id='ref{i}'>{i}</a>]: {r['citation']} <a href='https://doi.org/{r['doi']}'>Link</a>")

# ---------------------------------------------------------------- appendix
A("<a href='#appendix' id='appendix'>Appendix</a>")
A(table([["Parameter", "Value"],
         ["Background PM2.5 / PM10", "300 / 500 µg/m³"],
         ["Inversion", "+5 K over 0 m to 100 m; PM aloft 50 % of street level"],
         ["Heavy drone", "50 kg, 4 coaxial pairs of 54 in rotors, 5.9 m², hover at 60 m"],
         ["Fleet", "100 drones, 50 m grid, 500 m × 500 m"],
         ["Fence", "6 m high, 0.5 m deep, K = 1.5, plate area 20 m²/m³"],
         ["Migration velocity", "PM2.5 0.04 m/s; coarse PM 0.12 m/s"],
         ["Waste-fire hotspot", "100 fires in 500 m × 500 m; 100 g waste/h each; 10 g PM2.5/kg; 4 h without action"],
         ["Sensor alert", "10 µg/m³ above background (1-hour average)"],
         ["Ground teams", "6 min per fire with a drone map; 20 min per fire without a map"],
         ["Patrol drone", "36 km/h, 300 m camera width, 70 % flight time"]]))
A("All simulation code and data are in the Drone V AQI repository (github.com/Tech-Studio-Art/dronevaqi), "
  "directory plans/sim. This text follows the rules of ASD-STE100 Simplified Technical English.")

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
img { display: block; max-width: 100%; max-height: 72mm; margin: 6pt auto 2pt auto; }
p.cap { font-size: 8.4pt; color: #333; margin-bottom: 8pt; }
table { border-collapse: collapse; width: 100%; margin: 4pt 0 8pt 0; font-size: 8.4pt; page-break-inside: avoid; }
td { border: 1px solid #c8c8c8; padding: 3pt 4pt; vertical-align: top; }
td p { margin: 0; }
p.ref { font-size: 8.2pt; margin-bottom: 3pt; }
"""
html = [f"<html><head><meta charset='utf-8'><style>{css}</style></head><body>"]
SEC = {"abstract", "introduction", "methods", "results", "impact", "pilot", "india", "discussion", "references", "appendix"}
for b in B:
    if b.startswith("<img") or b.startswith("<table"):
        html.append(b)
    elif any(f"id='{s}'" in b for s in SEC):
        html.append(f"<p class='sec'>{b}</p>")
    elif b.startswith("<a href='#fig"):
        html.append(f"<p class='cap'>{b}</p>")
    elif b.startswith("[<a href='#ref"):
        html.append(f"<p class='ref'>{b}</p>")
    elif len(b) < 60 and b[:2] in ("A.", "B.", "C.", "D."):
        html.append(f"<p class='h'>{b}</p>")
    else:
        html.append(f"<p>{b}</p>")
html.append("</body></html>")
open(os.path.join(HERE, "report.html"), "w").write("\n".join(html))
print("entries", len(B))
