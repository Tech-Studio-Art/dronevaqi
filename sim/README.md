# Drone V AQI – first-order physics simulation

`python3 sim/drone_aqi_sim.py` (needs only numpy). Output is saved in `results.txt`.

Models: actuator-disk momentum theory for rotor airflow; single-fibre filtration theory
(diffusion, interception, impaction, electret polarisation + Coulomb capture) with Davies
pressure drop; ventilated well-mixed box models for Delhi NCT and a 200 m x 200 m zone;
Slinn (1983) drop–particle collection efficiency for mist scavenging.

Headline results (see results.txt):
- A 25 kg octocopter moves ~19 m³/s through its rotors. Any filter that catches smoke-sized
  PM2.5 (0.1–1 µm) at rotor face velocity (~5 m/s) costs 2–400x hover power.
- Best realistic case (pleated electret, 5x disk area): ~1.8x power, endurance 34 → 18 min,
  PM10 capture ~80 %, but PM2.5 mass capture only ~11 % → ~2 m³/s clean air per drone.
- 100 such drones over Delhi: < 0.01 % PM2.5 reduction. Local 200 m zone: ~1–4 %.
- Mist: PM10 capture plausible locally; PM2.5 needs millimetres of water over the whole city
  (billions of litres) – not viable as a city-scale method.
