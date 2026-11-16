
import csv, json, numpy as np
from pathlib import Path
from obspy import UTCDateTime as U
from obspy.taup import TauPyModel
from obspy.geodetics import locations2degrees, gps2dist_azimuth

ROOT = Path(__file__).resolve().parents[1]
REF = U("2024-01-01T07:10:00")  
USGS = dict(lat=37.487, lon=137.271, depth=10.0, time="2024-01-01T07:10:09")
rows = list(csv.DictReader(open(ROOT / "data/picks.csv")))
for r in rows:
    r["tp"] = U("2024-01-01T" + r["p_pick_utc"]) - REF
    r["ts"] = U("2024-01-01T" + r["s_pick_utc"]) - REF
    r["lat"] = float(r["lat"])
    r["lon"] = float(r["lon"])
m = TauPyModel("iasp91")
D = np.round(np.arange(0.1, 60.01, 0.05), 2)


def first(arr, names):
    return min(a.time for a in arr if a.name in names)


TP = []
TS = []
for d in D:
    a = m.get_travel_times(10, d, ["p", "P", "Pn", "s", "S", "Sn"])
    TP.append(first(a, ("p", "P", "Pn")))
    TS.append(first(a, ("s", "S", "Sn")))
TP = np.array(TP)
TS = np.array(TS)
np.savetxt(
    ROOT / "data/iasp91_traveltimes.csv",
    np.c_[D, TP, TS, TS - TP],
    delimiter=",",
    header="dist_deg,tP_s,tS_s,tSminusP_s",
    comments="",
    fmt="%.3f",
)
VP, VS = 8.0, 4.3
out = []
for r in rows:
    sp = r["ts"] - r["tp"]
    d_graph = float(np.interp(sp, TS - TP, D))
    d_const = sp / (1 / VS - 1 / VP)
    true = locations2degrees(USGS["lat"], USGS["lon"], r["lat"], r["lon"])
    out.append(
        dict(
            station=r["station"],
            net=r["net"],
            lat=r["lat"],
            lon=r["lon"],
            p=r["p_pick_utc"],
            s=r["s_pick_utc"],
            tp=round(r["tp"], 2),
            ts=round(r["ts"], 2),
            sp=round(sp, 2),
            d_const_km=round(d_const, 1),
            d_graph_deg=round(d_graph, 2),
            d_graph_km=round(d_graph * 111.195, 1),
            d_usgs_deg=round(true, 2),
            d_usgs_km=round(true * 111.195, 1),
            s_unc=float(r["s_unc_s"]),
            note=r["note"],
        )
    )

tp = np.array([o["tp"] for o in out])
sp = np.array([o["sp"] for o in out])
b, a = np.polyfit(tp, sp, 1)
pred = a + b * tp
r2 = 1 - np.sum((sp - pred) ** 2) / np.sum((sp - sp.mean()) ** 2)
x0 = -a / b

lats = np.arange(30, 45, 0.05)
lons = np.arange(128, 147, 0.05)
best = (1e9, None)
P_obs = np.array([o["tp"] for o in out])
sl = np.array([o["lat"] for o in out])
so = np.array([o["lon"] for o in out])
from obspy.geodetics.base import locations2degrees as l2d

G = np.zeros((len(lats), len(lons)))
for i, la in enumerate(lats):
    dd = l2d(la, lons[:, None], sl[None, :], so[None, :])
    tt = np.interp(dd, D, TP)
    t0 = np.mean(P_obs[None, :] - tt, axis=1)
    res = P_obs[None, :] - tt - t0[:, None]
    G[i] = np.sqrt(np.mean(res**2, axis=1))
i, j = np.unravel_index(np.argmin(G), G.shape)

lats2 = np.arange(lats[i] - 0.1, lats[i] + 0.1, 0.005)
lons2 = np.arange(lons[j] - 0.1, lons[j] + 0.1, 0.005)
bestv = 1e9
for la in lats2:
    dd = l2d(la, lons2[:, None], sl[None, :], so[None, :])
    tt = np.interp(dd, D, TP)
    t0 = np.mean(P_obs[None, :] - tt, axis=1)
    rms = np.sqrt(np.mean((P_obs[None, :] - tt - t0[:, None]) ** 2, axis=1))
    k = np.argmin(rms)
    if rms[k] < bestv:
        bestv = rms[k]
        bl = (la, lons2[k], t0[k])
gl_lat, gl_lon, gl_t0 = bl
dd = l2d(gl_lat, gl_lon, sl, so)
ptt = np.interp(dd, D, TP)
pres = P_obs - ptt - gl_t0
for o, rr in zip(out, pres):
    o["p_resid"] = round(float(rr), 2)
np.save(ROOT / "data/grid_rms.npy", G)
np.save(ROOT / "data/grid_axes.npy", np.array([lats[0], lats[-1], lons[0], lons[-1]]))

dg = np.array([o["d_graph_deg"] for o in out])
w = 1 / np.array([o["s_unc"] for o in out])
bestc = (1e9, None)
for la in np.arange(30, 45, 0.02):
    d2 = l2d(la, np.arange(128, 147, 0.02)[:, None], sl[None, :], so[None, :])
    mis = np.sqrt(np.sum((w * (d2 - dg)) ** 2, axis=1) / np.sum(w**2))
    k = np.argmin(mis)
    if mis[k] < bestc[0]:
        bestc = (mis[k], (la, 128 + 0.02 * k))
cl_lat, cl_lon = bestc[1]


def km(a1, o1, a2, o2):
    return gps2dist_azimuth(a1, o1, a2, o2)[0] / 1000


res = dict(
    ref_time=str(REF),
    stations=out,
    wadati=dict(
        slope=round(b, 4),
        intercept=round(a, 3),
        r2=round(r2, 4),
        x_intercept_s=round(x0, 2),
        origin_time=str(t0_wadati),
        vp_vs=round(1 + b, 3),
        err_s=round(t0_wadati - U(USGS["time"]), 2),
    ),
    circle_fit=dict(
        lat=round(cl_lat, 3),
        lon=round(cl_lon, 3),
        misfit_deg=round(bestc[0], 3),
        err_km=round(km(cl_lat, cl_lon, USGS["lat"], USGS["lon"]), 1),
    ),
    p_grid=dict(
        lat=round(gl_lat, 3),
        lon=round(gl_lon, 3),
        origin_time=str(REF + gl_t0),
        rms_s=round(bestv, 2),
        err_km=round(km(gl_lat, gl_lon, USGS["lat"], USGS["lon"]), 1),
        err_t_s=round(REF + gl_t0 - U(USGS["time"]), 2),
    ),
    usgs=USGS,
)
json.dump(res, open(ROOT / "data/results.json", "w"), indent=1)
with open(ROOT / "data/picks_and_distances.csv", "w", newline="") as f:
    w_ = csv.DictWriter(f, fieldnames=list(out[0].keys()))
    w_.writeheader()
    w_.writerows(out)
if __name__ == "__main__":
    for o in out:
        print(o["station"], o["sp"], o["d_const_km"], o["d_graph_km"], o["d_usgs_km"], o["p_resid"])
    print(json.dumps({k: v for k, v in res.items() if k != "stations"}, indent=1))
