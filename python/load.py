import numpy as np
from obspy import read
from obspy.signal.rotate import rotate2zne
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
ORDER = ["MAJO", "ERM", "INCN", "YSS", "TATO", "GUMO", "KMI", "KAPI", "COLA"]
AZ = {
    "KMI": (0, 90),
    "ERM": (359.9, 87.8),
    "KAPI": (144.6, 233.7),
    "COLA": (206, 296),
    "GUMO": (309, 39),
    "INCN": (0, 90),
    "MAJO": (358, 88),
    "TATO": (331, 61),
    "YSS": (0, 90),
}


def stations():
    import csv

    out = {}
    with open(ROOT / "data/stations.txt") as f:
        r = csv.reader(f, delimiter="|")
        next(r)
        for row in r:
            out[row[0]] = dict(net=row[2], lat=float(row[3]), lon=float(row[4]), name=row[8])
    return out


def load():
    st = read(str(ROOT / "data/noto_9stations.mseed"))
    st.merge()
    st.detrend("demean")
    st.detrend("linear")
    res = {}
    for s in ORDER:
        z = st.select(station=s, channel="BHZ")[0]
        h1 = st.select(station=s, channel="BH1")[0]
        h2 = st.select(station=s, channel="BH2")[0]
        n = min(len(z), len(h1), len(h2))
        a1, a2 = AZ[s]
        Z, N, E = rotate2zne(
            z.data[:n].astype(float),
            0,
            -90,
            h1.data[:n].astype(float),
            a1,
            0,
            h2.data[:n].astype(float),
            a2,
            0,
        )
        tr = {}
        for c, d in zip("ZNE", (Z, N, E)):
            t = z.copy()
            t.data = d
            t.stats.channel = "BH" + c
            tr[c] = t
        res[s] = tr
    return res
