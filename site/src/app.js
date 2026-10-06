(function () {
  const D = window.EPIWAVE_DATA, C = window.EPIWAVE_CONTENT, M = window.EpiModel;
  const R1 = D.replicate1;
  const $ = (sel, root = document) => root.querySelector(sel);
  const fmt = (v, d = 3) => (Math.abs(v) >= 1e4 || (Math.abs(v) < 1e-3 && v !== 0)) ? v.toExponential(2) : v.toFixed(d);
  const pct = v => (100 * v).toFixed(1) + "%";

  const store = {
    get(k, fallback) { try { const v = localStorage.getItem(k); return v === null ? fallback : JSON.parse(v); } catch (e) { return fallback; } },
    set(k, v) { try { localStorage.setItem(k, JSON.stringify(v)); } catch (e) { /* storage unavailable */ } }
  };

  // ---------------- theme ----------------
  function applyTheme(t) {
    if (t === "light" || t === "dark") document.documentElement.setAttribute("data-theme", t);
    else document.documentElement.removeAttribute("data-theme");
  }
  let theme = store.get("ew-theme", "system");
  applyTheme(theme);

  // ---------------- study summaries ----------------
  const pred = D.study.predictive;
  const reps = [...new Set(pred.rep)];
  const rmse = { with: {}, without: {} };
  pred.rep.forEach((r, i) => { rmse[pred.model[i]][r] = pred.rmse[i]; });
  const gains = reps.map(r => 1 - rmse.with[r] / rmse.without[r]);
  const wins = gains.filter(g => g > 0).length;
  const medianGain = M.median(gains);
  const recovery = D.study.recovery;
  const recRow = (model, param) => {
    const i = recovery.model.findIndex((m, k) => m === model && recovery.param[k] === param);
    return { bias: recovery.bias[i], coverage: recovery.coverage[i], rhat: recovery.median_rhat[i], ess: recovery.median_ess[i], truth: recovery.true[i] };
  };

  // ---------------- charts ----------------
  function niceTicks(lo, hi, n = 5) {
    const span = hi - lo || 1, step0 = span / n, mag = Math.pow(10, Math.floor(Math.log10(step0)));
    const step = [1, 2, 2.5, 5, 10].map(s => s * mag).find(s => s >= step0);
    const out = [];
    for (let v = Math.ceil(lo / step) * step; v <= hi + step * 1e-9; v += step) out.push(+v.toPrecision(10));
    return out;
  }
  function tickLabel(v) {
    if (v === 0) return "0";
    const a = Math.abs(v);
    if (a >= 1000 || a < 0.01) return v.toExponential(0).replace("e-", "e−");
    return String(+v.toPrecision(3));
  }

  // series: [{x:[], y:[], color:'--with', dash:false, width:2, label}], points: [{x,y,color,r}]
  function lineChart(opt) {
    const W = opt.w || 760, H = opt.h || 300, P = { l: 62, r: 16, t: 14, b: 42 };
    const all = opt.series.flatMap(s => s.y.map((y, i) => [s.x[i], y])).concat((opt.points || []).map(p => [p.x, p.y]));
    const ys = all.map(p => p[1]).filter(Number.isFinite), xs = all.map(p => p[0]);
    const logY = !!opt.logY;
    const tY = v => logY ? Math.log10(Math.max(v, 1e-12)) : v;
    let yLo = opt.yMin !== undefined ? tY(opt.yMin) : Math.min(...ys.map(tY));
    let yHi = opt.yMax !== undefined ? tY(opt.yMax) : Math.max(...ys.map(tY));
    if (!logY && opt.yMin === undefined) yLo = Math.min(yLo, 0);
    if (yHi === yLo) yHi = yLo + 1;
    const pad = (yHi - yLo) * 0.06; if (opt.yMax === undefined) yHi += pad; if (logY && opt.yMin === undefined) yLo -= pad;
    const xLo = opt.xMin !== undefined ? opt.xMin : Math.min(...xs), xHi = opt.xMax !== undefined ? opt.xMax : Math.max(...xs);
    const sx = x => P.l + (x - xLo) / (xHi - xLo || 1) * (W - P.l - P.r);
    const sy = v => H - P.b - (tY(v) - yLo) / (yHi - yLo) * (H - P.t - P.b);
    let g = "";
    const yt = logY ? Array.from({ length: Math.floor(yHi) - Math.ceil(yLo) + 1 }, (_, i) => Math.pow(10, Math.ceil(yLo) + i)) : niceTicks(yLo, yHi);
    yt.forEach(v => { const y = sy(v); g += `<line class="grid" x1="${P.l}" x2="${W - P.r}" y1="${y}" y2="${y}"/><text x="${P.l - 8}" y="${y + 4}" text-anchor="end">${tickLabel(v)}</text>`; });
    niceTicks(xLo, xHi, opt.xTicks || 6).forEach(v => { const x = sx(v); g += `<text x="${x}" y="${H - P.b + 18}" text-anchor="middle">${tickLabel(v)}</text>`; });
    g += `<line class="axis" x1="${P.l}" x2="${W - P.r}" y1="${H - P.b}" y2="${H - P.b}"/>`;
    if (opt.xLabel) g += `<text x="${(P.l + W - P.r) / 2}" y="${H - 6}" text-anchor="middle">${opt.xLabel}</text>`;
    if (opt.yLabel) g += `<text transform="translate(14 ${(P.t + H - P.b) / 2}) rotate(-90)" text-anchor="middle">${opt.yLabel}</text>`;
    (opt.vlines || []).forEach(v => { const x = sx(v.x); g += `<line x1="${x}" x2="${x}" y1="${P.t}" y2="${H - P.b}" style="stroke:var(${v.color || "--muted"});stroke-dasharray:4 4"/><text x="${x + 6}" y="${P.t + 12}">${v.label || ""}</text>`; });
    (opt.hlines || []).forEach(v => { const y = sy(v.y); g += `<line x1="${P.l}" x2="${W - P.r}" y1="${y}" y2="${y}" style="stroke:var(${v.color || "--muted"});stroke-dasharray:4 4"/>`; });
    opt.series.forEach(s => {
      const d = s.y.map((y, i) => (Number.isFinite(y) ? `${i && Number.isFinite(s.y[i - 1]) ? "L" : "M"}${sx(s.x[i]).toFixed(1)},${sy(y).toFixed(1)}` : "")).join("");
      g += `<path d="${d}" style="fill:none;stroke:var(${s.color});stroke-width:${s.width || 2};${s.dash ? "stroke-dasharray:6 5;" : ""}stroke-linejoin:round"/>`;
    });
    (opt.points || []).forEach(p => { g += `<circle cx="${sx(p.x).toFixed(1)}" cy="${sy(p.y).toFixed(1)}" r="${p.r || 3}" style="fill:var(${p.color || "--with"});fill-opacity:${p.opacity || 0.8}"/>`; });
    return `<svg class="chart" viewBox="0 0 ${W} ${H}" role="img" aria-label="${opt.label || "chart"}">${g}</svg>`;
  }

  function cssVar(name) { return getComputedStyle(document.documentElement).getPropertyValue(name).trim(); }
  function hexRgb(h) { h = h.replace("#", ""); if (h.length === 3) h = h.split("").map(c => c + c).join(""); const n = parseInt(h, 16); return [(n >> 16) & 255, (n >> 8) & 255, n & 255]; }
  // diverging colour: negative -> --without, positive -> --with, zero -> surface
  function divColor(v, lim) {
    const s = hexRgb(cssVar("--surface")), pos = hexRgb(cssVar("--with")), neg = hexRgb(cssVar("--without"));
    const t = Math.max(-1, Math.min(1, v / lim)), c = t >= 0 ? pos : neg, w = Math.abs(t);
    return s.map((sv, i) => Math.round(sv + (c[i] - sv) * w));
  }
  function heatmap(canvas, mat, lim) {
    const rows = mat.length, cols = mat[0].length, ctx = canvas.getContext("2d");
    canvas.width = cols; canvas.height = rows;
    const img = ctx.createImageData(cols, rows);
    mat.forEach((row, i) => row.forEach((v, j) => { const c = divColor(v, lim), k = 4 * (i * cols + j); img.data[k] = c[0]; img.data[k + 1] = c[1]; img.data[k + 2] = c[2]; img.data[k + 3] = 255; }));
    ctx.putImageData(img, 0, 0);
  }
  // likelihood surface on a log-compressed scale: the maximum is darkest, and
  // drops of 1, 10, 100, 1000 log-units fade in even steps down to the surface colour
  function surface(canvas, grid, lim) {
    const rows = grid.length, cols = grid[0].length, ctx = canvas.getContext("2d");
    canvas.width = cols; canvas.height = rows;
    const s = hexRgb(cssVar("--surface")), w = hexRgb(cssVar("--with")), img = ctx.createImageData(cols, rows);
    grid.forEach((row, i) => row.forEach((v, j) => {
      const t = Math.max(0, 1 - Math.log1p(-v) / Math.log1p(lim)), k = 4 * ((rows - 1 - i) * cols + j);
      s.forEach((sv, c) => { img.data[k + c] = Math.round(sv + (w[c] - sv) * t); }); img.data[k + 3] = 255;
    }));
    ctx.putImageData(img, 0, 0);
  }

  // ---------------- small builders ----------------
  const legend = items => `<div class="legend">${items.map(([label, color, dash]) => `<span class="${dash ? "dash" : ""}" style="--c:var(${color})">${label}</span>`).join("")}</div>`;
  const stat = (v, l, cls = "") => `<div class="stat ${cls}"><div class="v">${v}</div><div class="l">${l}</div></div>`;
  function slider(id, label, min, max, step, value, fmtFn = v => v) {
    return `<div class="control"><div class="row"><label for="${id}">${label}</label><output id="${id}-out">${fmtFn(value)}</output></div><input type="range" id="${id}" min="${min}" max="${max}" step="${step}" value="${value}"></div>`;
  }
  function bindSliders(root, ids, fmts, onChange) {
    ids.forEach(id => {
      const inp = $("#" + id, root), out = $("#" + id + "-out", root);
      inp.addEventListener("input", () => { out.textContent = (fmts[id] || (v => v))(+inp.value); onChange(); });
    });
  }
  const val = (root, id) => +$("#" + id, root).value;
  const months = R1.times.map((_, i) => i);

  // ---------------- labs ----------------
  const LABS = [
    { id: "fit", title: "The real fit", ch: "study", blurb: "Truth against both fitted models, site by site, from replicate 1 of the study.", mount: labFit },
    { id: "stage1", title: "Stage 1 in your hands", ch: "stage1", blurb: "Move m, a, g and bednet coverage; the Ross-Macdonald model re-solves live.", mount: labStage1 },
    { id: "field", title: "The residual field", ch: "stage2", blurb: "Draw ε from the GP + AR(1) and see what φ, θ and σ do.", mount: labField },
    { id: "prevalence", title: "From incidence to prevalence", ch: "likelihood", blurb: "The detectability kernel, and where the linear and bounded links part.", mount: labPrevalence },
    { id: "identify", title: "Separating α and γ", ch: "likelihood", blurb: "Likelihood surfaces from the real data: cases alone, surveys alone, both.", mount: labIdentify },
    { id: "ridge", title: "The σ²–θ ridge", ch: "stage2", blurb: "Where the chains actually went, and why sampling τ² should help.", mount: labRidge },
    { id: "holdout", title: "Held-out sites", ch: "study", blurb: "Where vector information should pay off: predicting sites with no data.", mount: labHoldout }
  ];

  function labFit(root) {
    root.innerHTML = `
      <div class="controls">
        <div class="control"><label for="fit-site">Site</label><select id="fit-site">${R1.I_true.map((_, i) => `<option value="${i}">Site ${i + 1}</option>`).join("")}</select></div>
        <div class="control"><label for="fit-scale">Scale</label><select id="fit-scale"><option value="log">Log</option><option value="lin">Linear</option></select></div>
        <p class="note">Posterior mean of the latent incidence from the greta fits saved with the 50-replicate study (4 chains × 2,000 samples).</p>
      </div>
      <div class="panels">
        ${legend([["Truth", "--truth"], ["With I*", "--with"], ["I* = 0", "--without", true]])}
        <div id="fit-chart"></div>
        <div class="stats" id="fit-stats"></div>
        <h3>All fifty replicates</h3>
        <p class="note">Each point is one replicate's map error. Points below the diagonal are replicates where the model with I* did better.</p>
        <div id="fit-reps"></div>
      </div>`;
    const draw = () => {
      const s = val(root, "fit-site"), log = $("#fit-scale", root).value === "log";
      $("#fit-chart", root).innerHTML = lineChart({
        label: "Latent incidence at one site", logY: log, xLabel: "Month", yLabel: "Infections per person per day",
        series: [{ x: months, y: R1.I_true[s], color: "--truth", width: 2.5 }, { x: months, y: R1.pred_with[s], color: "--with" }, { x: months, y: R1.pred_without[s], color: "--without", dash: true }]
      });
      const e = (p) => Math.sqrt(M.mean(p[s].map((v, t) => (v - R1.I_true[s][t]) ** 2)));
      $("#fit-stats", root).innerHTML = stat(fmt(e(R1.pred_with), 4), "Map error at this site, with I*", "with") + stat(fmt(e(R1.pred_without), 4), "Map error at this site, I* = 0", "without") +
        stat(fmt(rmse.with[1], 4) + " / " + fmt(rmse.without[1], 4), "Replicate 1, all sites (with / without)");
    };
    $("#fit-site", root).addEventListener("change", draw);
    $("#fit-scale", root).addEventListener("change", draw);
    draw();
    const xs = reps.map(r => rmse.without[r]), lo = Math.min(...xs, ...reps.map(r => rmse.with[r])), hi = Math.max(...xs, ...reps.map(r => rmse.with[r]));
    $("#fit-reps", root).innerHTML = lineChart({
      label: "Map error by replicate", w: 520, h: 340, xLabel: "Map error, I* = 0", yLabel: "Map error, with I*", xMin: lo, xMax: hi, yMin: lo, yMax: hi,
      series: [{ x: [lo, hi], y: [lo, hi], color: "--muted", width: 1, dash: true }],
      points: reps.map(r => ({ x: rmse.without[r], y: rmse.with[r], color: rmse.with[r] < rmse.without[r] ? "--with" : "--without", r: 4 }))
    }) + `<div class="stats">${stat(wins + " of " + reps.length, "Replicates where I* lowers map error", "with")}${stat(pct(medianGain), "Median reduction in map error")}</div>`;
  }

  function labStage1(root) {
    const F = { m: v => v.toFixed(1), amp: v => v.toFixed(2), a: v => v.toFixed(2), g: v => v.toFixed(3), itn: v => Math.round(v * 100) + "%", res: v => Math.round(v * 100) + "%", eip: v => v === 0 ? "off" : v + " days" };
    root.innerHTML = `
      <div class="controls">
        ${slider("s1-m", "Baseline m (mosquitoes per person)", 0.2, 4, 0.1, 2, F.m)}
        ${slider("s1-amp", "Seasonal amplitude", 0, 0.9, 0.05, 0.6, F.amp)}
        ${slider("s1-a", "Biting rate a (per day)", 0.1, 0.5, 0.01, 0.3, F.a)}
        ${slider("s1-g", "Mosquito death rate g (per day)", 0.05, 0.25, 0.005, 0.1, F.g)}
        ${slider("s1-itn", "ITN coverage reached by month 48", 0, 0.95, 0.05, 0.7, F.itn)}
        ${slider("s1-res", "Bioassay susceptibility (IR cube)", 0, 1, 0.05, 0.8, F.res)}
        ${slider("s1-eip", "Extrinsic incubation period", 0, 15, 1, 0, F.eip)}
        <button class="btn" id="s1-reset">Reset to the study's values</button>
      </div>
      <div class="panels">
        <div class="stats" id="s1-stats"></div>
        ${legend([["Human prevalence x", "--truth"], ["Mosquito infectiousness z", "--with", true]])}
        <div id="s1-xz"></div>
        <h3>I* = m · a · b · z</h3>
        <div id="s1-istar"></div>
        <p class="note">Solved with a fourth-order Runge–Kutta step of 6 hours; it matches R's deSolve output to within 0.002% (0.01% with an EIP). The EIP adds four exposed mosquito stages, so a share (1 + g n/4)<sup>−4</sup> of infected mosquitoes live to become infectious.</p>
      </div>`;
    const draw = () => {
      const p = { baselineM: val(root, "s1-m"), amplitude: val(root, "s1-amp"), baselineA: val(root, "s1-a"), baselineG: val(root, "s1-g"), itnMax: val(root, "s1-itn"), susceptibility: val(root, "s1-res"), eipDays: val(root, "s1-eip") || null };
      const s = M.stage1Site(p);
      const R0 = M.R0(p.baselineM, p.baselineA, 0.8, 0.8, p.baselineG, 1 / 7) * (p.eipDays ? M.eipSurvival(p.baselineG, p.eipDays) : 1);
      const end = M.applyItn(p.baselineM, p.baselineA, p.baselineG, p.itnMax, p.susceptibility);
      const eipFactor = p.eipDays ? M.eipSurvival(p.baselineG, p.eipDays) : 1;
      const R0end = M.R0(end.m, end.a, 0.8, 0.8, end.g, 1 / 7) * (p.eipDays ? M.eipSurvival(end.g, p.eipDays) : 1);
      $("#s1-stats", root).innerHTML = stat(R0.toFixed(2), "R₀ before nets") + stat(R0end.toFixed(2), "R₀ at final ITN coverage") +
        stat(pct(eipFactor), "Mosquitoes surviving the EIP") + stat(fmt(s.Istar[48], 4), "I* in month 48");
      $("#s1-xz", root).innerHTML = lineChart({ label: "Prevalence over time", xLabel: "Month", yLabel: "Proportion", yMin: 0, yMax: 1, series: [{ x: months, y: s.x, color: "--truth" }, { x: months, y: s.z, color: "--with", dash: true }] });
      $("#s1-istar", root).innerHTML = lineChart({ label: "I* over time", xLabel: "Month", yLabel: "Infections per person per day", series: [{ x: months, y: s.Istar, color: "--with" }] });
    };
    bindSliders(root, ["s1-m", "s1-amp", "s1-a", "s1-g", "s1-itn", "s1-res", "s1-eip"], { "s1-m": F.m, "s1-amp": F.amp, "s1-a": F.a, "s1-g": F.g, "s1-itn": F.itn, "s1-res": F.res, "s1-eip": F.eip }, draw);
    $("#s1-reset", root).addEventListener("click", () => {
      [["s1-m", 2, F.m], ["s1-amp", 0.6, F.amp], ["s1-a", 0.3, F.a], ["s1-g", 0.1, F.g], ["s1-itn", 0.7, F.itn], ["s1-res", 0.8, F.res], ["s1-eip", 0, F.eip]].forEach(([id, v, f]) => { $("#" + id, root).value = v; $("#" + id + "-out", root).textContent = f(v); });
      draw();
    });
    draw();
  }

  function labField(root) {
    let seed = 7;
    const maxD = Math.max(...R1.coords.flatMap(p => R1.coords.map(q => M.dist(p, q))));
    root.innerHTML = `
      <div class="controls">
        ${slider("f-sigma", "Innovation sd σ", 0.1, 1.5, 0.05, 0.6, v => v.toFixed(2))}
        ${slider("f-phi", "Lengthscale φ", 0.05, 3, 0.05, 3, v => v.toFixed(2))}
        ${slider("f-theta", "AR(1) correlation θ", 0, 0.98, 0.01, 0.75, v => v.toFixed(2))}
        <button class="btn" id="f-redraw">Draw a new field</button>
        <p class="note">Sites are replicate 1's ten locations. The study's true values are σ = 0.6, φ = 3, θ = 0.75.</p>
      </div>
      <div class="panels">
        <div class="stats" id="f-stats"></div>
        <div><canvas class="heat" id="f-heat" aria-label="Residual field, sites by months"></canvas>
        <p class="note">Rows are the ten sites, columns the 49 months. Green: incidence above the mechanistic prediction. Orange: below.</p></div>
        <h3>How far spatial correlation reaches</h3>
        <div id="f-corr"></div>
      </div>`;
    const canvas = $("#f-heat", root);
    canvas.style.aspectRatio = "49 / 18";
    const draw = () => {
      const sigma = val(root, "f-sigma"), phi = val(root, "f-phi"), theta = val(root, "f-theta");
      const eps = M.simulateField(R1.coords, 49, sigma, phi, theta, seed);
      const tau2 = sigma * sigma / (1 - theta * theta);
      heatmap(canvas, eps, 2.5 * Math.sqrt(tau2));
      $("#f-stats", root).innerHTML = stat(tau2.toFixed(2), "Stationary variance τ² = σ²/(1 − θ²)") +
        stat(M.matern52(maxD, phi).toFixed(2), "Correlation between the two farthest sites") + stat(Math.exp(Math.sqrt(tau2)).toFixed(2) + "×", "Incidence one sd above I*");
      const ds = Array.from({ length: 61 }, (_, i) => i * 0.025);
      $("#f-corr", root).innerHTML = lineChart({ label: "Matérn 5/2 correlation by distance", h: 240, xLabel: "Distance (unit square)", yLabel: "Correlation", yMin: 0, yMax: 1, series: [{ x: ds, y: ds.map(d => M.matern52(d, phi)), color: "--with" }], vlines: [{ x: maxD, label: "farthest pair" }] });
    };
    bindSliders(root, ["f-sigma", "f-phi", "f-theta"], { "f-sigma": v => v.toFixed(2), "f-phi": v => v.toFixed(2), "f-theta": v => v.toFixed(2) }, draw);
    $("#f-redraw", root).addEventListener("click", () => { seed += 1; draw(); });
    draw();
  }

  function labPrevalence(root) {
    const w = M.qMonthly(), W = w[0] + w[1];
    const toI = v => Math.pow(10, v);
    root.innerHTML = `
      <div class="controls">
        ${slider("p-log", "Incidence I (per person per day)", -4, -0.5, 0.05, -1.82, v => fmt(toI(v), 4))}
        <p class="note">Incidence is held constant over the two months that contribute. The study's simulated I* spends most of its time between 0.05 and 0.4.</p>
      </div>
      <div class="panels">
        <div class="stats" id="p-stats"></div>
        ${legend([["Linear: p = Λ", "--without", true], ["Bounded: p = 1 − e^−Λ", "--with"]])}
        <div id="p-chart"></div>
        <h3>The detectability kernel q(d)</h3>
        <div id="p-kernel"></div>
        <p class="note">Integrated over 30-day steps the weights are w₀ = ${w[0].toFixed(2)} and w₁ = ${w[1].toFixed(2)} days. Λ = I × (w₀ + w₁) when incidence is constant.</p>
      </div>`;
    const grid = Array.from({ length: 141 }, (_, i) => -4 + i * 0.025);
    const draw = () => {
      const I = toI(val(root, "p-log")), L = I * W, lin = L, bnd = 1 - Math.exp(-L);
      $("#p-stats", root).innerHTML = stat(L.toFixed(3), "Λ, expected detectable infections per person") + stat(lin.toFixed(3), "Linear prevalence", "without") +
        stat(bnd.toFixed(3), "Bounded prevalence", "with") + stat(pct((lin - bnd) / bnd), "Linear above bounded by");
      $("#p-chart", root).innerHTML = lineChart({
        label: "Prevalence against incidence", xLabel: "log₁₀ incidence (per person per day)", yLabel: "Prevalence", yMin: 0, yMax: 1.5,
        series: [{ x: grid, y: grid.map(g => Math.min(toI(g) * W, 1.5)), color: "--without", dash: true }, { x: grid, y: grid.map(g => 1 - Math.exp(-toI(g) * W)), color: "--with" }],
        hlines: [{ y: 1 }], vlines: [{ x: val(root, "p-log"), color: "--truth", label: "" }]
      });
    };
    const days = Array.from({ length: 36 }, (_, d) => d);
    $("#p-kernel", root).innerHTML = lineChart({ label: "Detectability by days since infection", h: 220, xLabel: "Days since infection", yLabel: "P(tests positive)", yMin: 0, yMax: 1, series: [{ x: days, y: days.map(d => M.qDaily(d)), color: "--with" }] });
    bindSliders(root, ["p-log"], { "p-log": v => fmt(toI(v), 4) }, draw);
    draw();
  }

  function labIdentify(root) {
    const ll = M.makeLogLik(R1), nA = 61, nG = 61;
    const aGrid = Array.from({ length: nA }, (_, i) => -1.2 + i * 0.04);
    const gGrid = Array.from({ length: nG }, (_, i) => Math.log(0.02) + i * (Math.log(0.5) - Math.log(0.02)) / (nG - 1));
    const cache = gGrid.map(lg => aGrid.map(a => ll(a, Math.exp(lg))));
    root.innerHTML = `
      <div class="controls">
        <div class="control"><label for="id-mode">Data used</label><select id="id-mode"><option value="cases">Case counts only</option><option value="surveys">Prevalence surveys only</option><option value="both" selected>Both</option></select></div>
        <p class="note">Log-likelihood of replicate 1's real simulated data across α and γ, with the true residual field held fixed. Darker means more likely, on a log scale: each step lighter is ten times further below the maximum, out to 5,000 log-likelihood units. The cross marks the truth (α = 0, γ = 0.1). Holding ε at its true value is an idealisation: in the real fit the surveys inform α + ε together, and the GP's zero-mean prior is what then pins α. The lab shows the geometry, not the posterior.</p>
      </div>
      <div class="panels">
        <div style="position:relative">
          <canvas class="heat" id="id-surf" style="aspect-ratio:1.6/1" aria-label="Likelihood surface over alpha and gamma"></canvas>
          <svg class="chart" viewBox="0 0 160 100" preserveAspectRatio="none" style="position:absolute;inset:0;width:100%;height:100%" aria-hidden="true">
            <line id="id-tx" x1="0" x2="160" style="stroke:var(--truth);stroke-width:0.4;stroke-dasharray:2 2"/>
            <line id="id-ty" y1="0" y2="100" style="stroke:var(--truth);stroke-width:0.4;stroke-dasharray:2 2"/>
          </svg>
        </div>
        <p class="note">Horizontal axis: α from −1.2 to 1.2. Vertical axis: γ from 0.02 to 0.5 on a log scale. With cases only, the high-likelihood region is a diagonal ridge: only γ·e<sup>α</sup> is pinned down. Surveys alone pin α. Together they meet at a single point.</p>
        <h3>Real posterior draws (with I*)</h3>
        <div id="id-draws"></div>
        <h3>The α identity</h3>
        <div id="id-hist"></div>
        <div class="stats" id="id-stats"></div>
      </div>`;
    const tx = $("#id-tx", root), ty = $("#id-ty", root);
    const yTrue = 100 - (Math.log(0.1) - gGrid[0]) / (gGrid[nG - 1] - gGrid[0]) * 100, xTrue = (0 - aGrid[0]) / (aGrid[nA - 1] - aGrid[0]) * 160;
    tx.setAttribute("y1", yTrue); tx.setAttribute("y2", yTrue); ty.setAttribute("x1", xTrue); ty.setAttribute("x2", xTrue);
    const draw = () => {
      const mode = $("#id-mode", root).value;
      const g = cache.map(row => row.map(v => mode === "cases" ? v.cases : mode === "surveys" ? v.surveys : v.cases + v.surveys));
      const mx = Math.max(...g.flat());
      surface($("#id-surf", root), g.map(row => row.map(v => v - mx)), 5000);
    };
    $("#id-mode", root).addEventListener("change", draw);
    draw();
    const dw = R1.draws_with, dn = R1.draws_without;
    $("#id-draws", root).innerHTML = lineChart({
      label: "Posterior draws of alpha and gamma", w: 560, h: 320, xLabel: "α", yLabel: "γ",
      yMin: Math.min(...dw.gamma_rr, R1.truth.gamma) * 0.98, yMax: Math.max(...dw.gamma_rr, R1.truth.gamma) * 1.02,
      xMin: Math.min(...dw.alpha, R1.truth.alpha_effective) - 0.02, xMax: Math.max(...dw.alpha, R1.truth.alpha_effective) + 0.02,
      series: [], points: dw.alpha.map((a, i) => ({ x: a, y: dw.gamma_rr[i], color: "--with", r: 2, opacity: 0.35 })),
      vlines: [{ x: R1.truth.alpha_effective, color: "--truth", label: "true α (centred)" }], hlines: [{ y: R1.truth.gamma, color: "--truth" }]
    });
    const bins = (xs, lo, hi, n) => { const c = new Array(n).fill(0); xs.forEach(x => { const k = Math.floor((x - lo) / (hi - lo) * n); if (k >= 0 && k < n) c[k]++; }); return c; };
    const lo = -2.2, hi = 0.4, n = 104, bx = Array.from({ length: n }, (_, i) => lo + (i + 0.5) * (hi - lo) / n);
    const shift = R1.mean_log_I_star;
    $("#id-hist", root).innerHTML = legend([["α draws, with I*", "--with"], ["α draws, I* = 0", "--without"]]) + lineChart({
      label: "Alpha draws under both models", h: 240, xLabel: "α", yLabel: "Draws", yMin: 0,
      series: [{ x: bx, y: bins(dw.alpha, lo, hi, n), color: "--with" }, { x: bx, y: bins(dn.alpha, lo, hi, n), color: "--without" }],
      vlines: [{ x: R1.truth.alpha_effective, color: "--truth", label: "true α" }, { x: R1.truth.alpha_effective + shift, color: "--muted", label: "true α + mean log I*" }]
    });
    $("#id-stats", root).innerHTML = stat(shift.toFixed(3), "mean(log I*) in replicate 1") + stat(M.mean(dn.alpha).toFixed(3), "Posterior mean α, I* = 0", "without") +
      stat((R1.truth.alpha_effective + shift).toFixed(3), "True α + mean(log I*)") + stat(recRow("without", "alpha").bias.toFixed(3), "\"Bias\" of α across 50 replicates, I* = 0");
  }

  function labRidge(root) {
    const dw = R1.draws_with, n = dw.theta.length, per = n / 4, chainOf = i => Math.floor(i / per);
    const colors = ["--with", "--without", "--chain3", "--chain4"];
    const tauDraws = dw.sigma2.map((s, i) => s / (1 - dw.theta[i] ** 2)), tauMed = M.median(tauDraws);
    const tauTrue = R1.truth.sigma2 / (1 - R1.truth.theta ** 2);
    root.innerHTML = `
      <div class="controls">
        ${slider("r-tau", "Dashed curve: constant τ²", 0.2, 3.5, 0.01, +tauMed.toFixed(2), v => v.toFixed(2))}
        <p class="note">Every point on the dashed curve σ² = τ²(1 − θ²) gives the field the same stationary variance τ², so data that inform only τ² cannot tell those points apart. The grey curve is the truth's τ² = ${tauTrue.toFixed(2)}.</p>
        <p class="note">Draws are from replicate 1's fit with I*, under the earlier priors (θ flat on 0–1), coloured by chain.</p>
      </div>
      <div class="panels">
        <div class="stats">${stat(tauMed.toFixed(2) + " vs " + tauTrue.toFixed(2), "τ²: posterior median vs truth")}${stat(M.median(dw.phi).toFixed(2) + " vs " + R1.truth.phi, "φ: posterior median vs truth")}${stat(recRow("with", "sigma2").rhat.toFixed(2) + " / " + recRow("with", "theta").rhat.toFixed(2), "Median R-hat over 50 replicates: σ² / θ")}</div>
        ${legend([["Chain 1", "--with"], ["Chain 2", "--without"], ["Chain 3", "--chain3"], ["Chain 4", "--chain4"], ["Truth", "--truth"]])}
        <div id="r-chart"></div>
        <p class="note">In this replicate the four chains land close to each other but far from the truth: the fitted field is about three times as variable as the true one and has a much shorter range (φ ≈ 0.7 against 3). Across the fifty replicates the chains for σ² and θ disagree with each other (median R-hat about 2.4), which is the ridge at work. Sampling τ² directly, with a Beta(2, 2) prior on θ, puts the data's information on one parameter. That change has not yet been tested at study scale.</p>
      </div>`;
    const th = Array.from({ length: 99 }, (_, i) => 0.01 + i * 0.01);
    const yMax = Math.max(...dw.sigma2, R1.truth.sigma2) * 1.25;
    const clip = y => (y <= yMax ? y : NaN);
    const draw = () => {
      const tau = val(root, "r-tau");
      $("#r-chart", root).innerHTML = lineChart({
        label: "Posterior draws of theta and sigma squared by chain", h: 340, xLabel: "θ", yLabel: "σ²", xMin: 0, xMax: 1, yMin: 0, yMax,
        series: [{ x: th, y: th.map(t => clip(tau * (1 - t * t))), color: "--truth", dash: true, width: 1.5 },
                 { x: th, y: th.map(t => clip(tauTrue * (1 - t * t))), color: "--muted", width: 1 }],
        points: dw.theta.map((t, i) => ({ x: t, y: dw.sigma2[i], color: colors[chainOf(i)], r: 2, opacity: 0.45 }))
          .concat([{ x: R1.truth.theta, y: R1.truth.sigma2, color: "--truth", r: 6, opacity: 1 }])
      });
    };
    bindSliders(root, ["r-tau"], { "r-tau": v => v.toFixed(2) }, draw);
    draw();
  }

  function labHoldout(root) {
    root.innerHTML = `
      <div class="controls">
        <div class="control"><label for="h-signal">Mechanistic prediction I*</label><select id="h-signal"><option value="flat">Same at every site (as in the study)</option><option value="vary" selected>Varies between sites</option></select></div>
        <div class="control"><label for="h-n">Number of sites</label><select id="h-n"><option>10</option><option selected>30</option><option>60</option></select></div>
        ${slider("h-hold", "Share of sites held out", 0.1, 0.5, 0.05, 0.3, v => Math.round(v * 100) + "%")}
        ${slider("h-phi", "True lengthscale φ", 0.1, 1.5, 0.05, 0.4, v => v.toFixed(2))}
        <button class="btn" id="h-redraw">New sites and field</button>
        <p class="note">A Gaussian approximation, not the greta model: log incidence is observed with small noise and the GP is predicted by kriging with known hyperparameters. It shows the direction of the effect, not the study's numbers.</p>
      </div>
      <div class="panels">
        <div class="stats" id="h-stats"></div>
        ${legend([["Truth", "--truth"], ["With I*", "--with"], ["I* = 0", "--without", true]])}
        <div id="h-chart"></div>
        <p class="note" id="h-caption"></p>
      </div>`;
    let seed = 11;
    const draw = () => {
      const n = +$("#h-n", root).value, hold = val(root, "h-hold"), phi = val(root, "h-phi"), vary = $("#h-signal", root).value === "vary";
      const R = M.rng(seed), coords = Array.from({ length: n }, () => [R.unif(), R.unif()]);
      const base = M.stage1Site({}).Istar.map(v => Math.log(Math.max(v, 1e-6)));
      // spatial signal in vector abundance: a smooth gradient plus a bump, in log units
      const siteEffect = coords.map(([x, y]) => vary ? 1.2 * (x - 0.5) + 0.8 * Math.sin(3 * y) - 0.4 : 0);
      const logIstar = coords.map((_, s) => base.map(b => b + siteEffect[s]));
      const sigma = 0.5, theta = 0.75, tau2 = sigma * sigma / (1 - theta * theta), noise = 0.02;
      const eps = M.simulateField(coords, 49, sigma, phi, theta, seed + 100);
      const truth = logIstar.map((row, s) => row.map((v, t) => v + eps[s][t]));
      const nHold = Math.max(1, Math.round(hold * n)), held = new Set(Array.from({ length: n }, (_, i) => i).sort(() => R.unif() - 0.5).slice(0, nHold));
      const obs = [...Array(n).keys()].filter(i => !held.has(i)), hid = [...held];
      const cObs = obs.map(i => coords[i]), cHid = hid.map(i => coords[i]);
      const predA = hid.map(() => []), predB = hid.map(() => []);
      // variance the I* = 0 model's field must carry: residual plus between-site spread of log I*
      const spread = M.mean(siteEffect.map(e => (e - M.mean(siteEffect)) ** 2));
      for (let t = 0; t < 49; t++) {
        const yObs = obs.map(i => truth[i][t] + Math.sqrt(noise) * R.norm());
        const resA = yObs.map((y, j) => y - logIstar[obs[j]][t]);
        const kA = M.krige(cObs, resA, cHid, tau2, phi, noise);
        const mu = M.mean(yObs), kB = M.krige(cObs, yObs.map(y => y - mu), cHid, tau2 + spread, phi, noise);
        hid.forEach((s, j) => { predA[j].push(logIstar[s][t] + kA[j]); predB[j].push(mu + kB[j]); });
      }
      const err = p => Math.sqrt(M.mean(hid.flatMap((s, j) => p[j].map((v, t) => (Math.exp(v) - Math.exp(truth[s][t])) ** 2))));
      const eA = err(predA), eB = err(predB);
      $("#h-stats", root).innerHTML = stat(fmt(eA, 4), "Map error at held-out sites, with I*", "with") + stat(fmt(eB, 4), "Map error at held-out sites, I* = 0", "without") +
        stat(pct(1 - eA / eB), "Reduction in map error with I*") + stat(nHold + " of " + n, "Sites held out");
      const j = 0, s = hid[j];
      $("#h-chart", root).innerHTML = lineChart({
        label: "Incidence at one held-out site", xLabel: "Month", yLabel: "Infections per person per day", logY: true,
        series: [{ x: months, y: truth[s].map(Math.exp), color: "--truth", width: 2.5 }, { x: months, y: predA[j].map(Math.exp), color: "--with" }, { x: months, y: predB[j].map(Math.exp), color: "--without", dash: true }]
      });
      $("#h-caption", root).textContent = vary
        ? "With I* varying between sites, the model with I* knows each held-out site's mechanistic level; the I* = 0 model can only borrow from neighbours."
        : "With I* the same everywhere, the I* = 0 model recovers the shared seasonal shape from the other sites each month, so the two models predict almost the same. This is the study's current design.";
    };
    ["h-signal", "h-n"].forEach(id => $("#" + id, root).addEventListener("change", draw));
    bindSliders(root, ["h-hold", "h-phi"], { "h-hold": v => Math.round(v * 100) + "%", "h-phi": v => v.toFixed(2) }, draw);
    $("#h-redraw", root).addEventListener("click", () => { seed += 1; draw(); });
    draw();
  }

  // ---------------- views ----------------
  const totals = {
    minutes: C.chapters.reduce((s, c) => s + c.minutes, 0),
    derivations: C.chapters.reduce((s, c) => s + c.derivations.length, 0),
    examples: C.chapters.reduce((s, c) => s + c.examples.length, 0),
    selftest: C.chapters.reduce((s, c) => s + c.selftest.length, 0)
  };
  const footer = `<p class="note">Numbers come from R/epiwave-foi-model.R (commit 64a06b9) and the 50-replicate simulation study run on ${D.study.meta.date}. The labs re-implement the model in JavaScript; their tests check it against R.</p>`;

  function viewHome(main) {
    main.innerHTML = `
      <section class="hero">
        <p class="eyebrow">EpiWave-FOI · vector-informed malaria mapping</p>
        <h1>Ten sites, forty-eight months, one mechanistic prediction.</h1>
        <p class="lede">Every part of this model is a way of turning mosquito numbers into a map of infection, and then asking how much the map should trust them. Read the chapter, work the derivation, move the sliders in the lab, then take the open questions to the next check-in.</p>
      </section>
      <section class="figure">
        <div class="fig-head"><h2>Replicate 1: truth and both fits</h2><a href="#lab-fit">Open the lab</a></div>
        <div class="fit-grid">
          <div style="min-width:0">${legend([["Truth", "--truth"], ["With I*", "--with"], ["I* = 0", "--without", true]])}<div id="home-chart"></div></div>
          <div class="readouts">
            <div class="readout with"><div class="v">${fmt(rmse.with[1], 4)}</div><div class="l">Map error with I*, replicate 1</div></div>
            <div class="readout without"><div class="v">${fmt(rmse.without[1], 4)}</div><div class="l">Map error with I* = 0</div></div>
            <div class="readout"><div class="v">${wins}/${reps.length}</div><div class="l">Replicates where I* lowers map error</div></div>
            <div class="readout"><div class="v">${pct(medianGain)}</div><div class="l">Median reduction in map error</div></div>
          </div>
        </div>
        <p class="note">Site 1, log scale. In this replicate the two fits are almost indistinguishable, and I* = 0 is marginally better. The mechanistic prediction is the same at every simulated site, so it carries no spatial information yet. The study chapter explains why, and the held-out lab shows what changes when it does.</p>
      </section>
      <section class="tiles">
        <div class="tile"><span class="v">${totals.minutes} min</span><span class="l">Reading</span></div>
        <div class="tile"><span class="v">${totals.derivations}</span><span class="l">Derivations, step by step</span></div>
        <div class="tile"><span class="v">${totals.examples}</span><span class="l">Worked examples with real numbers</span></div>
        <div class="tile"><span class="v">${totals.selftest}</span><span class="l">Self-test questions</span></div>
        <div class="tile"><span class="v">${LABS.length}</span><span class="l">Labs that run in the browser</span></div>
      </section>
      <section class="section"><h2>Chapters</h2><div class="cards">${C.chapters.map((c, i) => `<a class="card" href="#ch-${c.id}"><span class="meta">${i + 1} · ${c.minutes} min · ${c.derivations.length} derivations</span><h3>${c.title}</h3><p>${c.lede}</p></a>`).join("")}</div></section>
      <section class="section"><h2>Labs</h2><div class="cards">${LABS.map(l => `<a class="card" href="#lab-${l.id}"><h3>${l.title}</h3><p>${l.blurb}</p></a>`).join("")}</div></section>
      <section class="section"><h2>Before the next check-in</h2><p class="lede">${C.ask.length - doneAsk().size} open questions to raise. <a href="#ask">Open the tracker</a>.</p></section>
      ${footer}`;
    $("#home-chart", main).innerHTML = lineChart({ label: "Replicate 1, site 1", logY: true, xLabel: "Month", yLabel: "Infections per person per day", series: [{ x: months, y: R1.I_true[0], color: "--truth", width: 2.5 }, { x: months, y: R1.pred_with[0], color: "--with" }, { x: months, y: R1.pred_without[0], color: "--without", dash: true }] });
  }

  function viewChapter(main, id) {
    const i = C.chapters.findIndex(c => c.id === id), c = C.chapters[i];
    if (!c) return viewHome(main);
    const prev = C.chapters[i - 1], next = C.chapters[i + 1];
    const labs = LABS.filter(l => l.ch === c.id), asks = C.ask.filter(a => a.ch === c.id);
    main.innerHTML = `
      <section class="hero"><p class="eyebrow">Chapter ${i + 1} · ${c.minutes} min</p><h1>${c.title}</h1><p class="lede">${c.lede}</p></section>
      <section class="reading">${c.reading.map(p => `<p>${p}</p>`).join("")}</section>
      <section class="section reading"><h2>Derivations</h2>${c.derivations.map(d => `<details class="deriv"><summary>${d.title}</summary><ol>${d.steps.map(s => `<li>${s}</li>`).join("")}</ol></details>`).join("")}</section>
      <section class="section reading"><h2>Worked examples</h2>${c.examples.map(e => `<div class="example"><h3>${e.title}</h3><p>${e.body}</p></div>`).join("")}</section>
      <section class="section reading"><h2>Check yourself</h2>${c.selftest.map(q => `<div class="qa"><span class="tag">${q.depth}</span><p>${q.q}</p><button type="button">Show answer</button><p class="ans" hidden>${q.a}</p></div>`).join("")}</section>
      ${labs.length ? `<section class="section"><h2>Labs for this chapter</h2><div class="cards">${labs.map(l => `<a class="card" href="#lab-${l.id}"><h3>${l.title}</h3><p>${l.blurb}</p></a>`).join("")}</div></section>` : ""}
      ${asks.length ? `<section class="section reading"><h2>Open questions</h2><ul>${asks.map(a => `<li>${a.q}</li>`).join("")}</ul><p><a href="#ask">Track them</a></p></section>` : ""}
      <section class="stats">${prev ? `<a href="#ch-${prev.id}">← ${prev.title}</a>` : ""}${next ? `<a href="#ch-${next.id}">${next.title} →</a>` : ""}</section>
      ${footer}`;
    main.querySelectorAll(".qa button").forEach(b => b.addEventListener("click", () => { const a = b.nextElementSibling; a.hidden = !a.hidden; b.textContent = a.hidden ? "Show answer" : "Hide answer"; }));
  }

  function viewLabs(main) {
    main.innerHTML = `<section class="hero"><p class="eyebrow">Labs</p><h1>Seven models you can move</h1><p class="lede">Three use the real numbers from the study; four re-solve the model live as you change it.</p></section>
      <div class="cards">${LABS.map(l => `<a class="card" href="#lab-${l.id}"><span class="meta">${C.chapters.find(c => c.id === l.ch).short}</span><h3>${l.title}</h3><p>${l.blurb}</p></a>`).join("")}</div>${footer}`;
  }

  function viewLab(main, id) {
    const lab = LABS.find(l => l.id === id);
    if (!lab) return viewLabs(main);
    const ch = C.chapters.find(c => c.id === lab.ch);
    main.innerHTML = `<section class="hero"><p class="eyebrow">Lab · <a href="#ch-${ch.id}">${ch.short}</a></p><h1>${lab.title}</h1><p class="lede">${lab.blurb}</p></section><section class="lab" id="lab-root"></section>${footer}`;
    lab.mount($("#lab-root", main));
  }

  const doneAsk = () => new Set(store.get("ew-ask", []));
  function viewAsk(main) {
    let filter = "all";
    const render = () => {
      const done = doneAsk(), items = C.ask.filter(a => filter === "all" || a.ch === filter);
      main.innerHTML = `<section class="hero"><p class="eyebrow">Question tracker</p><h1>Questions for the next check-in</h1><p class="lede">${C.ask.length - done.size} of ${C.ask.length} still open. Ticks are saved in this browser only.</p></section>
        <div class="filters"><button class="btn" data-f="all" aria-pressed="${filter === "all"}">All</button>${C.chapters.filter(c => C.ask.some(a => a.ch === c.id)).map(c => `<button class="btn" data-f="${c.id}" aria-pressed="${filter === c.id}">${c.short}</button>`).join("")}</div>
        <div class="ask-list">${items.map(a => `<label class="ask-item ${done.has(a.id) ? "done" : ""}" for="ask-${a.id}"><input type="checkbox" id="ask-${a.id}" data-id="${a.id}" ${done.has(a.id) ? "checked" : ""}><div><p>${a.q}</p><span class="note">${C.chapters.find(c => c.id === a.ch).short}</span></div></label>`).join("")}</div>${footer}`;
      main.querySelectorAll(".filters button").forEach(b => b.addEventListener("click", () => { filter = b.dataset.f; render(); }));
      main.querySelectorAll(".ask-item input").forEach(inp => inp.addEventListener("change", () => {
        const d = doneAsk(); inp.checked ? d.add(inp.dataset.id) : d.delete(inp.dataset.id); store.set("ew-ask", [...d]); render();
      }));
    };
    render();
  }

  function viewFormulas(main) {
    main.innerHTML = `<section class="hero"><p class="eyebrow">Formula sheet</p><h1>Every equation, and where it lives in the code</h1></section>
      ${C.formulas.map(g => `<section class="formula-group"><h2>${g.group}</h2>${g.items.map(([name, f, fn]) => `<div class="formula-row"><span>${name}</span><span class="f">${f}</span><span class="fn">${fn || ""}</span></div>`).join("")}</section>`).join("")}${footer}`;
  }

  function viewDrill(main) {
    const deck = C.chapters.flatMap(c => c.selftest.map((q, i) => ({ ...q, id: c.id + "-" + i, ch: c.short })));
    let order = deck.map((_, i) => i).sort(() => Math.random() - 0.5), pos = 0, shown = false;
    const render = () => {
      const tally = store.get("ew-drill", {}), q = deck[order[pos]], t = tally[q.id] || { got: 0, again: 0 };
      const known = deck.filter(d => (tally[d.id] || {}).got > (tally[d.id] || {}).again).length;
      main.innerHTML = `<section class="hero"><p class="eyebrow">Drill</p><h1>Say the answer out loud, then check</h1><p class="lede">${known} of ${deck.length} questions answered correctly more often than not. Progress is saved in this browser only.</p></section>
        <div class="flash"><span class="eyebrow">${q.ch} · ${q.depth} · card ${pos + 1} of ${deck.length}</span><p class="q">${q.q}</p>
        ${shown ? `<p>${q.a}</p><div class="stats"><button class="btn primary" id="d-got">I had it</button><button class="btn" id="d-again">Not yet</button></div><p class="note">This card: ${t.got} right, ${t.again} not yet.</p>` : `<button class="btn primary" id="d-show">Show the answer</button>`}</div>${footer}`;
      if (!shown) $("#d-show", main).addEventListener("click", () => { shown = true; render(); });
      else ["got", "again"].forEach(k => $("#d-" + k, main).addEventListener("click", () => {
        const tl = store.get("ew-drill", {}); tl[q.id] = tl[q.id] || { got: 0, again: 0 }; tl[q.id][k]++; store.set("ew-drill", tl);
        pos = (pos + 1) % deck.length; shown = false; render();
      }));
    };
    render();
  }

  // ---------------- shell and routing ----------------
  function buildNav() {
    const nav = $("nav.side");
    nav.innerHTML = `<a class="brand" href="#home">EpiWave-FOI study</a>
      <div class="group"><span class="eyebrow">Chapters</span>${C.chapters.map((c, i) => `<a class="item" href="#ch-${c.id}"><span>${i + 1}. ${c.short}</span><span class="k">${c.minutes}m</span></a>`).join("")}</div>
      <div class="group"><span class="eyebrow">Labs</span>${LABS.map(l => `<a class="item" href="#lab-${l.id}"><span>${l.title}</span></a>`).join("")}</div>
      <div class="group"><span class="eyebrow">Practice</span><a class="item" href="#ask"><span>Question tracker</span></a><a class="item" href="#formulas"><span>Formula sheet</span></a><a class="item" href="#drill"><span>Drill</span></a></div>
      <button class="theme-toggle" id="theme-btn" type="button"></button>`;
    const btn = $("#theme-btn");
    const label = () => { btn.textContent = "Theme: " + theme; };
    btn.addEventListener("click", () => { theme = theme === "system" ? "light" : theme === "light" ? "dark" : "system"; store.set("ew-theme", theme); applyTheme(theme); label(); route(); });
    label();
  }

  function route() {
    const h = (location.hash || "#home").slice(1), main = $("main");
    if (h.startsWith("ch-")) viewChapter(main, h.slice(3));
    else if (h === "labs") viewLabs(main);
    else if (h.startsWith("lab-")) viewLab(main, h.slice(4));
    else if (h === "ask") viewAsk(main);
    else if (h === "formulas") viewFormulas(main);
    else if (h === "drill") viewDrill(main);
    else viewHome(main);
    document.querySelectorAll("nav.side a.item").forEach(a => {
      if (a.getAttribute("href") === "#" + h) a.setAttribute("aria-current", "page");
      else a.removeAttribute("aria-current");
    });
    window.scrollTo(0, 0);
  }

  buildNav();
  window.addEventListener("hashchange", route);
  window.matchMedia("(prefers-color-scheme: dark)").addEventListener("change", () => { if (theme === "system") route(); });
  route();
})();
