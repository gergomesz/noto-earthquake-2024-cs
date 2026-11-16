import sys, json, numpy as np
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from load import load, ORDER, ROOT
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.patches import Polygon
from obspy import UTCDateTime as U
from obspy.imaging.beachball import beach
from obspy.geodetics import locations2degrees, kilometers2degrees

THEMES = {
    "light": dict(
        ink="#121316", ink2="#4c4f55", ink3="#8d9096", rule="#dcdcd8", wash="#e9e9e4", paper="#fbfbf9"
    ),
    "dark": dict(
        ink="#ecebe6", ink2="#b4b4ae", ink3="#7c7d7a", rule="#2c2d30", wash="#232428", paper="#121315"
    ),
}
SANS = "Inter" if any("Inter" in f.name for f in font_manager.fontManager.ttflist) else "DejaVu Sans"

R = json.load(open(ROOT / "data/results.json"))
S = R["stations"]
REF = U(R["ref_time"])
POL = json.load(open(ROOT / "data/polarities.json"))
POLV = {"MAJO": -1, "YSS": 1, "KMI": 1, "KAPI": 1, "COLA": 1}
BAZ = {
    "MAJO": 323,
    "ERM": 213,
    "INCN": 83,
    "YSS": 189,
    "TATO": 47,
    "GUMO": 13,
    "KMI": 29,
    "KAPI": 33,
    "COLA": 291,
}
TT = np.loadtxt(ROOT / "data/iasp91_traveltimes.csv", delimiter=",", skiprows=1)
LAND = json.load(open(ROOT / "data/coast50.geojson"))
UO = R["usgs"]
CF = R["circle_fit"]
PG = R["p_grid"]
W = R["wadati"]
NP1 = (47, 37, 100)
NP2 = (215, 54, 82)
D = load()


def style(T):
    plt.rcParams.update(
        {
            "font.family": SANS,
            "mathtext.fontset": "stix",
            "font.size": 10.5,
            "text.color": T["ink"],
            "axes.labelcolor": T["ink"],
            "axes.edgecolor": T["ink"],
            "axes.linewidth": 0.9,
            "xtick.color": T["ink2"],
            "ytick.color": T["ink2"],
            "xtick.labelsize": 9.5,
            "ytick.labelsize": 9.5,
            "axes.facecolor": "none",
            "figure.facecolor": "none",
            "savefig.facecolor": "none",
            "savefig.transparent": True,
            "grid.color": T["rule"],
            "grid.linewidth": 0.8,
            "legend.frameon": False,
            "legend.fontsize": 9.5,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.labelsize": 9.5,
        }
    )


def lab(ax, s, which="x"):
    (ax.set_xlabel if which == "x" else ax.set_ylabel)(
        s.upper(), fontsize=8.8, fontweight="medium", labelpad=8
    )


def save(fig, out, name):
    fig.savefig(out / name, dpi=200, bbox_inches="tight", pad_inches=0.05)
    plt.close(fig)


def trace(s, c, f1, f2, a, b):
    if c != "T":
        tr = D[s][c].copy().filter("bandpass", freqmin=f1, freqmax=f2)
    else:
        bz = np.radians(BAZ[s])
        n = D[s]["N"].copy().filter("bandpass", freqmin=f1, freqmax=f2)
        e = D[s]["E"].copy().filter("bandpass", freqmin=f1, freqmax=f2)
        tr = n.copy()
        tr.data = n.data * np.sin(bz) - e.data * np.cos(bz)
    tr.trim(a, b)
    return tr.times() + (tr.stats.starttime - REF), tr.data


def drawland(ax, ext, T):
    for f in LAND["features"]:
        g = f["geometry"]
        polys = g["coordinates"] if g["type"] == "MultiPolygon" else [g["coordinates"]]
        for p in polys:
            xy = np.array(p[0])
            if (
                xy[:, 0].max() < ext[0] - 5
                or xy[:, 0].min() > ext[1] + 5
                or xy[:, 1].max() < ext[2] - 5
                or xy[:, 1].min() > ext[3] + 5
            ):
                continue
            ax.add_patch(Polygon(xy, closed=True, fc=T["wash"], ec=T["ink3"], lw=0.5, zorder=0))
    ax.set_xlim(ext[:2])
    ax.set_ylim(ext[2:])
    ax.set_aspect(1 / np.cos(np.radians(np.mean(ext[2:]))))
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(True)


def circle_pts(lat, lon, r_km, n=721):
    d = np.radians(kilometers2degrees(r_km))
    la = np.radians(lat)
    lo = np.radians(lon)
    az = np.radians(np.linspace(0, 360, n))
    la2 = np.arcsin(np.sin(la) * np.cos(d) + np.cos(la) * np.sin(d) * np.cos(az))
    lo2 = lo + np.arctan2(np.sin(az) * np.sin(d) * np.cos(la), np.cos(d) - np.sin(la) * np.sin(la2))
    lo2 = np.degrees(lo2) % 360
    la2 = np.degrees(la2)
    j = np.where(np.abs(np.diff(lo2)) > 90)[0] + 1
    return np.insert(la2, j, np.nan), np.insert(lo2.astype(float), j, np.nan)


def marks(ax, T, which=("usgs", "cf", "pg"), legend=True):
    if "usgs" in which:
        ax.plot(
            UO["lon"],
            UO["lat"],
            "*",
            ms=17,
            color=T["ink"],
            mec=T["paper"],
            mew=0.8,
            zorder=6,
            label="USGS epicentre",
        )
    if "cf" in which:
        ax.plot(
            CF["lon"],
            CF["lat"],
            "x",
            ms=11,
            mew=2.2,
            color=T["ink"],
            zorder=6,
            label=f"S−P circles, best fit ({CF['err_km']:.0f} km off)",
        )
    if "pg" in which:
        ax.plot(
            PG["lon"],
            PG["lat"],
            "o",
            ms=10,
            mfc=T["paper"],
            mec=T["ink"],
            mew=1.8,
            zorder=6,
            label=f"P-arrival grid search ({PG['err_km']:.0f} km off)",
        )
    if legend:
        ax.legend(
            loc="lower right",
            handletextpad=0.6,
            borderaxespad=0.8,
            frameon=True,
            facecolor=T["paper"],
            edgecolor="none",
            framealpha=0.94,
        )


def radP(sdr, az, toa):
    s, d, r = np.radians(sdr)
    az = np.radians(az)
    i = np.radians(toa)
    return (
        np.cos(r) * np.sin(d) * np.sin(i) ** 2 * np.sin(2 * (az - s))
        - np.cos(r) * np.cos(d) * np.sin(2 * i) * np.cos(az - s)
        + np.sin(r) * np.sin(2 * d) * (np.cos(i) ** 2 - np.sin(i) ** 2 * np.sin(az - s) ** 2)
        + np.sin(r) * np.cos(2 * d) * np.sin(2 * i) * np.sin(az - s)
    )


def gc(strike, dip):
    t = np.radians(np.linspace(0, 180, 181))
    s = np.radians(strike)
    d = np.radians(dip)
    vs = np.array([np.cos(s), np.sin(s), 0])
    vd = np.array([np.cos(s + np.pi / 2) * np.cos(d), np.sin(s + np.pi / 2) * np.cos(d), np.sin(d)])
    v = np.outer(np.cos(t), vs) + np.outer(np.sin(t), vd)
    pl = np.arcsin(np.clip(v[:, 2], -1, 1))
    tr = np.arctan2(v[:, 1], v[:, 0])
    r = (np.pi / 2 - pl) / (np.pi / 2)
    return r * np.sin(tr), r * np.cos(tr)


def build(theme):
    T = THEMES[theme]
    style(T)
    out = ROOT / "figures" / theme
    out.mkdir(parents=True, exist_ok=True)

    def seis(stns, name, ncol=1):
        rows = int(np.ceil(len(stns) / ncol))
        fig, axs = plt.subplots(rows, ncol, figsize=(11 if ncol == 1 else 15, 2.25 * rows), squeeze=False)
        for ax, s in zip(axs.ravel(), stns):
            o = next(x for x in S if x["station"] == s)
            tp, ts = o["tp"], o["ts"]
            a = REF + tp - 15
            b = REF + ts + 45
            for k, (c, col, lw) in enumerate(
                (
                    ("Z", T["ink2"], 0.55),
                    ("N", T["ink2"], 0.55),
                    ("E", T["ink2"], 0.55),
                    ("T", T["ink"], 0.75),
                )
            ):
                x, y = trace(s, c, 0.05, 2, a, b)
                ax.plot(x, y / np.abs(y).max() * 0.85 - 2 * k, color=col, lw=lw)
                ax.text(
                    x[0] - 0.012 * (x[-1] - x[0]),
                    -2 * k,
                    c,
                    ha="right",
                    va="center",
                    fontsize=9,
                    color=col,
                    fontweight="medium",
                )
            ax.axvline(tp, color=T["ink"], lw=1.3)
            ax.axvline(ts, color=T["ink"], lw=1.3, ls=(0, (4, 3)))
            ax.text(tp, 1.25, " P", fontsize=11, family="STIXGeneral", style="italic")
            ax.text(ts, 1.25, " S", fontsize=11, family="STIXGeneral", style="italic")
            ax.set_yticks([])
            ax.spines["left"].set_visible(False)
            ax.set_ylim(-7, 1.7)
            ax.set_title(
                f"{s}    P {o['p']}    S {o['s']}",
                loc="left",
                fontsize=9.5,
                fontweight="medium",
                color=T["ink"],
            )
            lab(ax, "seconds after 07:10:00 UTC")
        for ax in axs.ravel()[len(stns) :]:
            ax.axis("off")
        fig.tight_layout(h_pad=1.6)
        save(fig, out, name)

    seis(["MAJO", "INCN", "KAPI", "COLA"], "fig1_seismograms.png")
    seis(ORDER, "fig1b_all_stations.png", ncol=3)
    fig, ax = plt.subplots(figsize=(10, 6.2))
    ax.plot(TT[:, 0], TT[:, 1] / 60, color=T["ink"], lw=1.6, label="P  (iasp91)")
    ax.plot(TT[:, 0], TT[:, 2] / 60, color=T["ink"], lw=1.6, ls=(0, (6, 3)), label="S  (iasp91)")
    ax.plot(TT[:, 0], TT[:, 3] / 60, color=T["ink2"], lw=1.4, ls=(0, (1.5, 2)), label="S − P")
    for i, o in enumerate(S, 1):
        x, y = o["d_graph_deg"], o["sp"] / 60
        ax.plot([x, x], [0, y], color=T["ink3"], lw=0.7)
        ax.plot(x, y, "o", ms=7, color=T["ink"], mec=T["paper"], mew=1)
        ax.annotate(
            str(i), (x, y), xytext=(5, -13), textcoords="offset points", fontsize=9.5, fontweight="medium"
        )
    ax.set_xlim(0, 60)
    ax.set_ylim(0, 20)
    ax.set_xticks(range(0, 61, 5))
    ax.set_yticks(range(0, 21, 2))
    ax.grid(True)
    lab(ax, "epicentral distance (degrees)")
    lab(ax, "travel time (minutes)", "y")
    ax.legend(loc="upper left")
    save(fig, out, "fig2_traveltime.png")
    for key, name in (("d_const_km", "fig3a_map_const.png"), ("d_graph_km", "fig3b_map_graph.png")):
        fig, ax = plt.subplots(figsize=(9, 7))
        ext = (95, 215, -12, 72)
        drawland(ax, ext, T)
        for i, o in enumerate(S, 1):
            la, lo = circle_pts(o["lat"], o["lon"], o[key])
            ax.plot(lo, la, color=T["ink"], lw=0.9)
            x = o["lon"] % 360
            ax.plot(x, o["lat"], "^", ms=7, color=T["ink"], mec=T["paper"], mew=0.8, zorder=5)
            ax.annotate(
                f"{i} {o['station']}",
                (x, o["lat"]),
                xytext=(6, -12) if o["station"] == "MAJO" else (5, 4),
                textcoords="offset points",
                fontsize=8.5,
                fontweight="medium",
            )
        marks(ax, T, ("usgs", "cf") if key == "d_graph_km" else ("usgs",))
        ax.set_xticks(range(100, 215, 20))
        ax.set_xticklabels([f"{v}°E" if v <= 180 else f"{360-v}°W" for v in range(100, 215, 20)])
        lab(ax, "longitude")
        lab(ax, "latitude (°N)", "y")
        save(fig, out, name)
    fig, ax = plt.subplots(figsize=(8.6, 7.4))
    ext = (128, 146, 30, 44)
    drawland(ax, ext, T)
    for o in S:
        la, lo = circle_pts(o["lat"], o["lon"], o["d_graph_km"])
        ax.plot(lo, la, color=T["ink"], lw=1.0)
    marks(ax, T)
    lab(ax, "longitude (°E)")
    lab(ax, "latitude (°N)", "y")
    save(fig, out, "fig3c_map_zoom.png")
    tp = np.array([o["tp"] for o in S])
    sp = np.array([o["sp"] for o in S])
    ot = U(round(U(W["origin_time"]).timestamp, 2))
    fig, ax = plt.subplots(figsize=(10, 6.2))
    x = np.linspace(W["x_intercept_s"], tp.max() + 40, 10)
    ax.plot(
        x,
        W["intercept"] + W["slope"] * x,
        color=T["ink"],
        lw=1.5,
        label=f"best fit   $y = {W['slope']:.4f}x + {W['intercept']:.3f}$   $R^2 = {W['r2']:.4f}$",
    )
    ax.errorbar(
        tp,
        sp,
        yerr=[o["s_unc"] for o in S],
        fmt="o",
        ms=7,
        mfc=T["paper"],
        mec=T["ink"],
        mew=1.4,
        ecolor=T["ink2"],
        elinewidth=0.9,
        capsize=3,
        label="stations  (± S-pick uncertainty)",
    )
    for o in S:
        ax.annotate(
            o["station"],
            (o["tp"], o["sp"]),
            xytext=(8, -4),
            textcoords="offset points",
            fontsize=8.8,
            fontweight="medium",
        )
    ax.plot(W["x_intercept_s"], 0, "D", ms=8, color=T["ink"], zorder=5)
    ax.annotate(
        f"x-intercept = {W['x_intercept_s']:.2f} s\norigin time {ot.strftime('%H:%M:%S.%f')[:11]} UTC",
        (W["x_intercept_s"], 0),
        xytext=(200, 25),
        textcoords="offset points",
        fontsize=10,
        fontweight="medium",
        arrowprops=dict(arrowstyle="-", color=T["ink"], lw=0.8),
    )
    ax.set_xlim(W["x_intercept_s"] - 15, tp.max() + 40)
    ax.set_ylim(-10, 500)
    ax.grid(True)
    ax.legend(loc="upper left")
    lab(ax, "P-wave arrival (s after 07:10:00 UTC)")
    lab(ax, "S − P lag (s)", "y")
    save(fig, out, "fig4_wadati.png")
    fig, axs = plt.subplots(1, 2, figsize=(11, 5.6))
    for k, ax in enumerate(axs):
        ax.set_xlim(-1.2, 1.2)
        ax.set_ylim(-1.2, 1.25)
        ax.set_aspect(1)
        ax.axis("off")
        fc = T["ink"] if k == 0 else ("#d2d2cd" if theme == "light" else "#3d3e42")
        ax.add_collection(
            beach(
                NP1,
                xy=(0, 0),
                width=2,
                facecolor=fc,
                bgcolor=T["paper"],
                edgecolor=T["ink"],
                linewidth=1.2,
                zorder=1,
            )
        )
        ax.text(0, 1.1, "N", ha="center", fontsize=10, fontweight="medium")
    for s, p in POL.items():
        toa, az = p["takeoff"], p["az"]
        if toa > 90:
            toa = 180 - toa
            az = (az + 180) % 360
        r = toa / 90
        x = r * np.sin(np.radians(az))
        y = r * np.cos(np.radians(az))
        v = POLV.get(s)
        if v == 1:
            axs[0].plot(x, y, "o", ms=12, mfc=T["paper"], mec=T["paper"], zorder=5)
            axs[0].plot(x, y, "o", ms=5, color=T["ink"], zorder=6)
        elif v == -1:
            axs[0].plot(x, y, "s", ms=11, mfc=T["paper"], mec=T["ink"], mew=1.8, zorder=5)
        else:
            axs[0].plot(x, y, "o", ms=10, mfc="none", mec=T["ink3"], mew=1.6, ls="none", zorder=5)
        axs[0].annotate(
            s,
            (x, y),
            xytext=(8, 5),
            textcoords="offset points",
            fontsize=8.5,
            fontweight="medium",
            color=T["ink"],
            bbox=dict(fc=T["paper"], ec="none", pad=1.2, alpha=0.9),
            zorder=7,
        )
    x, y = gc(*NP1[:2])
    axs[1].plot(x, y, color=T["ink"], lw=2.8, zorder=5)
    x, y = gc(*NP2[:2])
    axs[1].plot(x, y, color=T["ink"], lw=1.6, ls=(0, (4, 3)), zorder=5)
    for stk, col in ((NP1[0], T["ink"]), (NP2[0], T["ink3"])):
        a = np.radians(stk)
        axs[1].plot([np.sin(a), 1.1 * np.sin(a)], [np.cos(a), 1.1 * np.cos(a)], color=col, lw=2)
        axs[1].text(
            1.2 * np.sin(a),
            1.2 * np.cos(a),
            f"{stk}°",
            ha="center",
            va="center",
            fontsize=11,
            family="STIXGeneral",
            color=col,
        )
    axs[0].text(
        -1.2,
        -1.24,
        "●  compression (up)     ■  dilatation (down)     ○  unclear",
        fontsize=8.5,
        color=T["ink2"],
    )
    from matplotlib.lines import Line2D

    axs[1].legend(
        [Line2D([], [], color=T["ink"], lw=2.6), Line2D([], [], color=T["ink"], lw=1.6, ls=(0, (4, 3)))],
        ["plane 1   47° / 37° / 100°   fault (dips SE)", "plane 2   215° / 54° / 82°   auxiliary (dips NW)"],
        loc="upper center",
        bbox_to_anchor=(0.5, 0.0),
        fontsize=8.5,
        labelcolor=T["ink2"],
    )
    save(fig, out, "fig5_focal_mechanism.png")
    G = np.load(ROOT / "data/grid_rms.npy")
    ga = np.load(ROOT / "data/grid_axes.npy")
    fig, ax = plt.subplots(figsize=(8.6, 7.4))
    ext = (128, 146, 30, 44)
    drawland(ax, ext, T)
    cmap = matplotlib.colors.LinearSegmentedColormap.from_list("ink", [T["ink"], T["paper"]])
    gx = np.linspace(ga[2], ga[3], G.shape[1])
    gy = np.linspace(ga[0], ga[1], G.shape[0])
    ax.contourf(gx, gy, G, levels=[0, 3, 5, 10, 15, 20], cmap=cmap, alpha=0.55, zorder=1)
    cs = ax.contour(gx, gy, G, levels=[3, 5, 10, 15, 20], colors=T["ink"], linewidths=0.8, zorder=2)
    ax.clabel(cs, fmt="%d s", fontsize=8.5)
    marks(ax, T)
    lab(ax, "longitude (°E)")
    lab(ax, "latitude (°N)", "y")
    save(fig, out, "fig6_p_gridsearch.png")
    lag = []
    dd = []
    for o in S:
        d = locations2degrees(PG["lat"], PG["lon"], o["lat"], o["lon"])
        lag.append(o["sp"] - np.interp(d, TT[:, 0], TT[:, 3]))
        dd.append(d)
    fig, ax = plt.subplots(figsize=(10, 5))
    ax.axhspan(15, 20, fc="none", ec=T["ink3"], hatch="////", lw=0, zorder=0)
    ax.axhspan(15, 20, fc="none", ec=T["ink2"], lw=0.8)
    ax.text(
        63.5,
        21.2,
        "slow initial rupture, 15–20 s (Ma et al., 2024)",
        ha="right",
        va="bottom",
        fontsize=8.8,
        color=T["ink2"],
    )
    ax.axhline(0, color=T["ink"], lw=0.9)
    ax.errorbar(
        dd,
        lag,
        yerr=[o["s_unc"] for o in S],
        fmt="o",
        ms=7,
        mfc=T["paper"],
        mec=T["ink"],
        mew=1.4,
        ecolor=T["ink2"],
        elinewidth=0.9,
        capsize=3,
    )
    for o, x, y in zip(S, dd, lag):
        ax.annotate(
            o["station"], (x, y), xytext=(8, 3), textcoords="offset points", fontsize=8.8, fontweight="medium"
        )
    ax.set_xlim(-1, 64)
    ax.grid(True, axis="y")
    lab(ax, "distance from epicentre (°)")
    lab(ax, "S − P excess (s)", "y")
    save(fig, out, "fig7_s_delay.png")
    return dict(zip([o["station"] for o in S], [round(float(v), 1) for v in lag]))


if __name__ == "__main__":
    for th in THEMES:
        ex = build(th)
    json.dump(ex, open(ROOT / "data/s_excess.json", "w"))
    print("figures done; S excess", ex)
