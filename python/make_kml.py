# makes the google earth file (stations + circles)
# the answer key folder starts hidden so you can place your own pin first
import json, numpy as np
from pathlib import Path
from obspy.geodetics import kilometers2degrees

ROOT = Path(__file__).resolve().parents[1]
R = json.load(open(ROOT / "data/results.json"))


def circle(lat, lon, r_km, n=360):
    d = np.radians(kilometers2degrees(r_km))
    la = np.radians(lat)
    lo = np.radians(lon)
    az = np.radians(np.linspace(0, 360, n + 1))
    la2 = np.arcsin(np.sin(la) * np.cos(d) + np.cos(la) * np.sin(d) * np.cos(az))
    lo2 = lo + np.arctan2(np.sin(az) * np.sin(d) * np.cos(la), np.cos(d) - np.sin(la) * np.sin(la2))
    return " ".join(f"{(np.degrees(b)+180)%360-180:.5f},{np.degrees(a):.5f},0" for a, b in zip(la2, lo2))


def pm(name, lat, lon, style, desc="", vis=1):
    return f"<Placemark><name>{name}</name><visibility>{vis}</visibility><description><![CDATA[{desc}]]></description><styleUrl>#{style}</styleUrl><Point><coordinates>{lon},{lat},0</coordinates></Point></Placemark>"


def ring(name, lat, lon, r, style, vis=1):
    return f"<Placemark><name>{name}</name><visibility>{vis}</visibility><styleUrl>#{style}</styleUrl><LineString><tessellate>1</tessellate><coordinates>{circle(lat,lon,r)}</coordinates></LineString></Placemark>"


S = R["stations"]
U = R["usgs"]
CF = R["circle_fit"]
PG = R["p_grid"]
k = [
    '<?xml version="1.0" encoding="UTF-8"?><kml xmlns="http://www.opengis.net/kml/2.2"><Document><name>Noto 2024 — triangulation</name>',
    '<Style id="sta"><IconStyle><color>ffffffff</color><scale>1.1</scale><Icon><href>http://maps.google.com/mapfiles/kml/shapes/triangle.png</href></Icon></IconStyle></Style>',
    '<Style id="graph"><LineStyle><color>ff6b6bff</color><width>2.5</width></LineStyle></Style>',
    '<Style id="const"><LineStyle><color>ff4da9ff</color><width>2</width></LineStyle></Style>',
    '<Style id="usgs"><IconStyle><color>ffaef069</color><scale>1.4</scale><Icon><href>http://maps.google.com/mapfiles/kml/shapes/star.png</href></Icon></IconStyle></Style>',
    '<Style id="pg"><IconStyle><color>fff7c34f</color><scale>1.2</scale><Icon><href>http://maps.google.com/mapfiles/kml/shapes/placemark_circle.png</href></Icon></IconStyle></Style>',
    '<Style id="cf"><IconStyle><color>ff6b6bff</color><scale>1.2</scale><Icon><href>http://maps.google.com/mapfiles/kml/shapes/cross-hairs.png</href></Icon></IconStyle></Style>',
]
k.append(
    "<Folder><name>1. Seismic stations</name>"
    + "".join(
        pm(f"{i} {o['station']}", o["lat"], o["lon"], "sta", f"P {o['p']} · S {o['s']} · S−P {o['sp']} s")
        for i, o in enumerate(S, 1)
    )
    + "</Folder>"
)
k.append(
    "<Folder><name>2. Circles — travel-time graph (iasp91) [red]</name>"
    + "".join(
        ring(f"{o['station']} {o['d_graph_km']:.0f} km", o["lat"], o["lon"], o["d_graph_km"], "graph")
        for o in S
    )
    + "</Folder>"
)
k.append(
    "<Folder><name>3. Circles — constant velocity 8 / 4.3 km/s [orange]</name><visibility>0</visibility>"
    + "".join(
        ring(f"{o['station']} {o['d_const_km']:.0f} km", o["lat"], o["lon"], o["d_const_km"], "const", 0)
        for o in S
    )
    + "</Folder>"
)
k.append(
    "<Folder><name>4. ANSWER KEY (turn on after placing your pin)</name><visibility>0</visibility>"
    + pm("USGS epicentre", U["lat"], U["lon"], "usgs", "37.487°N 137.271°E, 07:10:09 UTC", 0)
    + pm(f"S−P circles best fit ({CF['err_km']} km off)", CF["lat"], CF["lon"], "cf", "", 0)
    + pm(f"P-arrival grid search ({PG['err_km']} km off)", PG["lat"], PG["lon"], "pg", "", 0)
    + "</Folder>"
)
k.append("</Document></kml>")
(ROOT / "downloads").mkdir(parents=True, exist_ok=True)
open(ROOT / "downloads/noto_triangulation.kml", "w").write("\n".join(k))
print("ok")
