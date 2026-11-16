(function () {
  "use strict";
  // dark mode button (same as the puna site)
  (function themeToggle() {
    const btn = document.getElementById("theme-toggle");
    if (!btn) return;
    const root = document.documentElement;
    const mq = window.matchMedia ? window.matchMedia("(prefers-color-scheme: dark)") : null;
    let stored = null;
    try { stored = localStorage.getItem("noto-theme"); } catch (e) { /* storage unavailable */ }
    if (stored === "dark" || stored === "light") root.setAttribute("data-theme", stored);
    const current = () => root.getAttribute("data-theme") || (mq && mq.matches ? "dark" : "light");
    const label = () => {
      const dark = current() === "dark";
      btn.setAttribute("aria-pressed", String(dark));
      btn.querySelector("span").textContent = dark ? "Light mode" : "Dark mode";
    };
    btn.addEventListener("click", () => {
      const next = current() === "dark" ? "light" : "dark";
      root.setAttribute("data-theme", next);
      try { localStorage.setItem("noto-theme", next); } catch (e) {}
      label();
    });
    if (mq && mq.addEventListener) mq.addEventListener("change", label);
    label();
  })();

  const REDUCED = window.matchMedia && window.matchMedia("(prefers-reduced-motion: reduce)").matches;
  function drawOn(root) {
    if (REDUCED || !root) return;
    root.querySelectorAll("[data-draw]").forEach((p) => {
      if (!p.getTotalLength) return;
      const L = p.getTotalLength();
      p.style.transition = "none"; p.style.strokeDasharray = L + " " + L; p.style.strokeDashoffset = L;
      p.getBoundingClientRect();
      requestAnimationFrame(() => {
        p.style.transition = ""; p.style.strokeDashoffset = 0;
        const done = () => { p.style.strokeDasharray = ""; p.style.strokeDashoffset = ""; };
        p.addEventListener("transitionend", done, { once: true }); setTimeout(done, 2200);
      });
    });
  }
  document.querySelectorAll(".eli-btn").forEach((b) => b.addEventListener("click", () => {
    const panel = document.getElementById(b.getAttribute("aria-controls"));
    if (!panel) return;
    const open = b.getAttribute("aria-expanded") !== "true";
    b.setAttribute("aria-expanded", String(open)); panel.hidden = !open;
  }));
  document.querySelectorAll("[data-replay]").forEach((b) => b.addEventListener("click", () => drawOn(document.getElementById(b.dataset.replay))));
  drawOn(document.getElementById("xsec"));

  // slider: S-P lag to distance
  const T = window.TT; const sl = document.getElementById("sp-slider");
  if (T && sl) {
    const VP = 8.0, VS = 4.3;
    const interp = (x, xs, ys) => {
      if (x <= xs[0]) return ys[0];
      for (let i = 1; i < xs.length; i++) if (x <= xs[i]) { const t = (x - xs[i-1]) / (xs[i] - xs[i-1]); return ys[i-1] + t * (ys[i] - ys[i-1]); }
      return ys[ys.length - 1];
    };
    const fmt = (n) => Math.round(n).toLocaleString("en");
    const upd = () => {
      const sp = +sl.value; const deg = interp(sp, T.sp, T.d); const g = deg * 111.195; const c = sp / (1 / VS - 1 / VP);
      document.getElementById("sp-v").textContent = sp + " s";
      document.getElementById("sp-graph").textContent = fmt(g);
      document.getElementById("sp-deg").textContent = deg.toFixed(1) + "°";
      document.getElementById("sp-const").textContent = fmt(c) + " km";
      document.getElementById("sp-diff").textContent = (c - g >= 0 ? "+" : "−") + fmt(Math.abs(c - g)) + " km";
    };
    sl.addEventListener("input", upd); upd();
  }
})();
