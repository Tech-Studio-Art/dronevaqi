#!/usr/bin/env python3
"""Report 2 references: fetch each DOI from Crossref, build the citation from Crossref's own
metadata, and fail if the registered title does not match the title we cite."""
import difflib, html, json, re, sys, time, urllib.request
REFS = [  # key, DOI, title we intend to cite (order = first citation in the report)
 ("guttikunda12", "10.1007/s10661-011-2182-8", "Role of meteorology in seasonality of air pollution in megacity Delhi, India"),
 ("mahrt", "10.1146/annurev-fluid-010313-141354", "Stably stratified atmospheric boundary layers"),
 ("stevens", "10.1146/annurev-fluid-010816-060206", "Flow structure and turbulence in wind farms"),
 ("porteagel20", "10.1007/s10546-019-00473-0", "Wind-turbine and wind-farm flows: a review"),
 ("calaf", "10.1063/1.3291077", "Large eddy simulation study of fully developed wind-turbine array boundary layers"),
 ("allaerts", "10.1007/s10546-017-0307-5", "Gravity waves and wind-farm efficiency in neutral and stable conditions"),
 ("ribeiro", "10.1016/j.agrformet.2006.08.019", "Apple orchard frost protection with wind machine operation"),
 ("battany", "10.1016/j.agrformet.2012.01.009", "Vineyard frost protection with upward-blowing wind machines"),
 ("briggs72", "10.1016/0004-6981(72)90120-5", "Chimney plumes in neutral and stable surroundings"),
 ("wright", "10.1061/(asce)0733-9429(1984)110:5(643)", "Buoyant jets in density-stratified crossflow"),
 ("launder", "10.1016/0045-7825(74)90029-2", "The numerical computation of turbulent flows"),
 ("businger", "10.1175/1520-0469(1971)028<0181:fprita>2.0.co;2", "Flux-profile relationships in the atmospheric surface layer"),
 ("dyer", "10.1007/bf00240838", "A review of flux-profile relationships"),
 ("hogstrom", "10.1007/bf00119875", "Non-dimensional wind and temperature profiles in the atmospheric surface layer: A re-evaluation"),
 ("rodi", "10.1029/jc092ic05p05305", "Examples of calculation methods for flow and mixing in stratified fluids"),
 ("richards", "10.1016/0167-6105(93)90124-7", "Appropriate boundary conditions for computational wind engineering models using the k-ε turbulence model"),
 ("hargreaves", "10.1016/j.jweia.2006.08.002", "On the use of the k–ε model in commercial CFD software to model the neutral atmospheric boundary layer"),
 ("blockenabl", "10.1016/j.atmosenv.2006.08.019", "CFD simulation of the atmospheric boundary layer: wall function problems"),
 ("franke", "10.1504/ijep.2011.038443", "The COST 732 Best Practice Guideline for CFD simulation of flows in the urban environment: a summary"),
 ("porteagel11", "10.1016/j.jweia.2011.01.011", "Large-eddy simulation of atmospheric boundary layer flow through wind turbines and wind farms"),
 ("mtt", "10.1098/rspa.1956.0011", "Turbulent gravitational convection from maintained and instantaneous sources"),
 ("bloomfield", "10.1017/s0022112097008252", "Turbulent fountains in a stratified fluid"),
 ("choutapalli", "10.1017/s0022112009007009", "An experimental study of an axisymmetric turbulent pulsed air jet"),
 ("hanna", "10.2172/5591108", "Handbook on atmospheric diffusion"),
]
norm = lambda s: re.sub(r"[^a-z0-9]", "", html.unescape(re.sub(r"<[^>]+>", "", s)).lower().replace("ε", "epsilon"))
out, bad = [], []
for key, doi, want in REFS:
    m = None
    for t in range(4):
        try:
            m = json.load(urllib.request.urlopen(urllib.request.Request(
                f"https://api.crossref.org/works/{doi}", headers={"User-Agent": "dronevaqi-refcheck"}), timeout=30))["message"]
            break
        except Exception as e:
            print("retry", key, e, file=sys.stderr)
            time.sleep(6 * (t + 1))
    if m is None:
        bad.append((key, "not found", 0)); continue
    title = html.unescape(re.sub(r"<[^>]+>", "", (m.get("title") or [""])[0]))
    sim = difflib.SequenceMatcher(None, norm(title), norm(want)).ratio()
    au = m.get("author", [])
    names = [f"{a.get('family','').title() if a.get('family','').isupper() else a.get('family','')} "
             f"{''.join(p[0] for p in re.split(r'[ .-]+', a.get('given','')) if p)}".strip() for a in au]
    auth = ", ".join(names[:6]) + (", et al." if len(names) > 6 else "")
    if not au:
        auth = {"hanna": "Hanna SR, Briggs GA, Hosker RP Jr", "rodi": "Rodi W"}.get(key, "")
    auth = auth.replace("Hosker RPJ", "Hosker RP Jr")
    jn = html.unescape((m.get("container-title") or [""])[0])
    yr = (m.get("issued", {}).get("date-parts") or [[None]])[0][0]
    vol, iss, pg = m.get("volume"), m.get("issue"), m.get("page")
    cit = f"{auth}{'' if auth.endswith('.') else '.'} {title}. {jn or m.get('publisher', '')}"
    if vol: cit += f" {vol}"
    if iss: cit += f"({iss})"
    if pg: cit += f":{pg}"
    cit += f" ({yr})."
    cit = cit.replace("k– model", "k–ε model").replace("k-ϵ", "k-ε")      # symbol lost in Crossref metadata
    ok = sim >= 0.9
    if not ok: bad.append((key, title, sim))
    out.append(dict(key=key, doi=doi, crossref_title=title, cited_title=want, match=round(sim, 3), citation=cit))
    print(f"{'OK ' if ok else 'BAD'} {sim:.2f} {key}: {cit[:160]}")
    time.sleep(1.2)
json.dump(out, open("refs_verified.json", "w"), indent=1, ensure_ascii=False)
print("mismatches:", bad)
