import json
from pathlib import Path
from obspy import UTCDateTime as U

ROOT = Path(__file__).resolve().parents[1]
R = json.load(open(ROOT / "data/results.json"))
S = R["stations"]
W = R["wadati"]
CF = R["circle_fit"]
PG = R["p_grid"]
EX = json.load(open(ROOT / "data/s_excess.json"))


def dms(x, h):
    d = int(abs(x))
    mf = (abs(x) - d) * 60
    m = int(mf)
    s = round((mf - m) * 60, 1)
    if s >= 60:
        s -= 60
        m += 1
    return f"{d}°{m:02d}′{s:04.1f}″{h}"


_POL = {
    "MAJO": "down (dilatation)",
    "YSS": "up (compression)",
    "KMI": "up (compression)",
    "KAPI": "up (compression)",
    "COLA": "up (compression)",
}
t1 = "".join(
    f"<tr><td>{i}</td><td class='em'>{o['station']}</td><td>{o['net']}</td><td class='num'>{o['lat']:.3f}</td><td class='num'>{o['lon']:.3f}</td><td>{o['p']}</td><td>{o['s']}</td><td class='num'>{o['s_unc']:g} s</td><td>{o['note']}</td></tr>"
    for i, o in enumerate(S, 1)
)
t2 = "".join(
    f"<tr><td class='em'>{o['station']}</td><td class='num'>{o['tp']:.2f}</td><td class='num'>{o['ts']:.2f}</td><td class='num'>{o['sp']:.2f}</td><td class='num'>{o['d_const_km']:,.0f}</td><td class='num em'>{o['d_graph_km']:,.0f}</td><td class='num'>{o['d_graph_deg']:.2f}°</td></tr>"
    for o in S
)
t3 = "".join(
    f"<tr><td class='em'>{o['station']}</td><td class='num'>{o['d_graph_km']:,.0f}</td><td class='num'>{o['d_usgs_km']:,.0f}</td><td class='num'>{o['d_graph_km']-o['d_usgs_km']:+,.0f} km</td><td class='num'>{EX[o['station']]:+.1f} s</td><td class='num'>{o['p_resid']:+.2f} s</td></tr>"
    for o in S
)
import numpy as np

TT = np.loadtxt(ROOT / "data/iasp91_traveltimes.csv", delimiter=",", skiprows=1)[::5]
(ROOT / "assets").mkdir(exist_ok=True)
open(ROOT / "assets/data.js", "w").write(
    "window.TT = "
    + json.dumps({"d": [round(float(x), 2) for x in TT[:, 0]], "sp": [round(float(x), 2) for x in TT[:, 3]]})
    + ";\n"
)
ot = lambda s: (lambda t: t.strftime("%H:%M:") + f"{t.second+t.microsecond/1e6:05.2f}")(
    U(round(U(s).timestamp, 2))
)
rep = dict(
    T1=t1,
    T2=t2,
    T3=t3,
    W_SLOPE=f"{W['slope']:.4f}",
    W_INT=f"{W['intercept']:.3f}",
    W_R2=f"{W['r2']:.4f}",
    W_X0=f"{W['x_intercept_s']:.2f}",
    W_OT=ot(W["origin_time"]),
    W_ERR=f"{abs(W['err_s']):.1f}",
    W_VPVS=f"{W['vp_vs']:.2f}",
    CF_LAT=dms(CF["lat"], "N"),
    CF_LON=dms(CF["lon"], "E"),
    CF_ERR=f"{CF['err_km']:.0f}",
    PG_LAT=dms(PG["lat"], "N"),
    PG_LON=dms(PG["lon"], "E"),
    PG_ERR=f"{PG['err_km']:.0f}",
    PG_OT=ot(PG["origin_time"]),
    PG_TERR=f"{PG['err_t_s']:.1f}",
    PG_RMS=f"{PG['rms_s']:.2f}",
    EX_MED=f"{sorted(EX.values())[len(EX)//2]:.0f}",
)
h = open(ROOT / "site_src/template.html").read()
for k, v in rep.items():
    h = h.replace("{{" + k + "}}", v)
assert "{{" not in h, [l for l in h.splitlines() if "{{" in l][:3]
open(ROOT / "index.html", "w").write(h)
print("index.html written")
