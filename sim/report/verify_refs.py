#!/usr/bin/env python3
"""Fetch every reference from Crossref by DOI, build the citation from Crossref's own
metadata, and fail if the registered title does not match the title we cite."""
import difflib, html, json, re, time, urllib.request
REFS = [  # key, DOI, title we intend to cite
 ("pope", "10.1080/10473289.2006.10464485", "Health effects of fine particulate air pollution: lines that connect"),
 ("guttikunda", "10.1016/j.envdev.2012.12.002", "Health impacts of particulate pollution in a megacity—Delhi, India"),
 ("cusworth", "10.1088/1748-9326/aab303", "Quantifying the influence of agricultural fires in northwest India on urban air pollution in Delhi, India"),
 ("stull", "10.1007/978-94-009-3027-8", "An Introduction to Boundary Layer Meteorology"),
 ("mahrt", "10.1146/annurev-fluid-010313-141354", "Stably stratified atmospheric boundary layers"),
 ("blocken", "10.1016/j.buildenv.2015.02.015", "Computational Fluid Dynamics for urban physics: Importance, scales, possibilities, limitations and ten tips and tricks towards accurate and reliable simulations"),
 ("launder", "10.1016/0045-7825(74)90029-2", "The numerical computation of turbulent flows"),
 ("turner", "10.1017/s0022112066001526", "Jets and plumes with negative or reversing buoyancy"),
 ("list", "10.1146/annurev.fl.14.010182.001201", "Turbulent jets and plumes"),
 ("zilitinkevich", "10.1007/s10546-012-9768-8", "A hierarchy of energy- and flux-budget (EFB) turbulence closure models for stably-stratified geophysical flows"),
 ("richards", "10.1016/0167-6105(93)90124-7", "Appropriate boundary conditions for computational wind engineering models using the k-ε turbulence model"),
 ("tominaga", "10.1016/j.atmosenv.2007.06.054", "Turbulent Schmidt numbers for CFD analysis with various types of flowfield"),
 ("mizuno", "10.1109/94.879357", "Electrostatic precipitation"),
 ("white", "10.1109/t-aiee.1951.5060545", "Particle charging in electrostatic precipitation"),
 ("jaworek", "10.1016/j.elstat.2006.07.012", "Modern electrostatic devices and methods for exhaust gas cleaning: A brief review"),
 ("baldauf", "10.1016/j.atmosenv.2008.05.051", "Impacts of noise barriers on near-road air quality"),
 ("hagler", "10.1016/j.atmosenv.2011.02.030", "Model evaluation of roadside barrier impact on near-road air pollution"),
 ("hanna", "10.2172/5591108", "Handbook on atmospheric diffusion"),
 ("wiedinmyer", "10.1021/es502250z", "Global emissions of trace gases, particulate matter, and hazardous air pollutants from open burning of domestic waste"),
 ("snyder", "10.1021/es4022602", "The changing paradigm of air pollution monitoring"),
 ("schroeder", "10.1016/j.rse.2013.12.008", "The New VIIRS 375 m active fire detection data product: Algorithm description and initial assessment"),
 ("battany", "10.1016/j.agrformet.2012.01.009", "Vineyard frost protection with upward-blowing wind machines"),
 ("boelter", "10.1080/02786829708965505", "Ozone generation by indoor, electrostatic air cleaners"),
]
norm = lambda s: re.sub(r"[^a-z0-9]", "", html.unescape(s).lower().replace("ε", "epsilon"))
out, bad = [], []
for key, doi, want in REFS:
    for t in range(4):
        try:
            m = json.load(urllib.request.urlopen(urllib.request.Request(
                f"https://api.crossref.org/works/{doi}", headers={"User-Agent": "dronevaqi-refcheck"}), timeout=30))["message"]
            break
        except Exception as e:
            time.sleep(6 * (t + 1))
    title = html.unescape((m.get("title") or [""])[0])
    sim = difflib.SequenceMatcher(None, norm(title), norm(want)).ratio()
    au = m.get("author", [])
    names = [f"{a.get('family','').title() if a.get('family','').isupper() else a.get('family','')} "
             f"{''.join(p[0] for p in re.split(r'[ .-]+', a.get('given','')) if p)}".strip() for a in au]
    auth = ", ".join(names[:6]) + (", et al." if len(names) > 6 else "")
    if not au:
        auth = {"hanna": "Hanna SR, Briggs GA, Hosker RP", "stull": "Stull RB"}.get(key, "")
    jn = html.unescape((m.get("container-title") or [""])[0])
    yr = (m.get("issued", {}).get("date-parts") or [[None]])[0][0]
    vol, iss, pg = m.get("volume"), m.get("issue"), m.get("page")
    src = jn or m.get("publisher", "")
    cit = f"{auth}{'' if auth.endswith('.') else '.'} {title}. {src}"
    if vol: cit += f" {vol}"
    if iss: cit += f"({iss})"
    if pg: cit += f":{pg}"
    cit += f" ({yr})."
    ok = sim >= 0.9
    if not ok: bad.append((key, title, sim))
    out.append(dict(key=key, doi=doi, crossref_title=title, cited_title=want, match=round(sim, 3), citation=cit))
    print(f"{'OK ' if ok else 'BAD'} {sim:.2f} {key}: {cit[:150]}")
    time.sleep(1.2)
json.dump(out, open("refs_verified.json", "w"), indent=1, ensure_ascii=False)
print("mismatches:", bad)
