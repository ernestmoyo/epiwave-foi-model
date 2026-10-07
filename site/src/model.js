// The model maths from R/epiwave-foi-model.R, ported to JavaScript so the labs
// run in the browser. Tested against R in site/tests/model.test.js.
(function (root) {
  const M = {};

  // ---- random numbers (seeded, so a lab redraws the same field) ----
  M.rng = function (seed) {
    let s = seed >>> 0;
    const unif = () => {
      s = (s + 0x6d2b79f5) >>> 0;
      let t = s;
      t = Math.imul(t ^ (t >>> 15), t | 1);
      t ^= t + Math.imul(t ^ (t >>> 7), t | 61);
      return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
    };
    let spare = null;
    const norm = () => {
      if (spare !== null) { const v = spare; spare = null; return v; }
      let u = 0, v = 0;
      while (u === 0) u = unif();
      v = unif();
      const r = Math.sqrt(-2 * Math.log(u));
      spare = r * Math.sin(2 * Math.PI * v);
      return r * Math.cos(2 * Math.PI * v);
    };
    return { unif, norm };
  };

  // ---- Stage 1 ----
  M.seasonalM = function (tDays, baselineM, amplitude, phase = 0) {
    const year = tDays / 365;
    const s = 1 + amplitude * (Math.sin(4 * Math.PI * year + phase) - 0.3 * Math.cos(2 * Math.PI * year));
    return baselineM * Math.max(s, 0.1);
  };

  // Share of the ITN effect kept under resistance (Symons et al., as in the
  // Vector Atlas IR cube): 1 - 0.46 (1 - susceptibility).
  M.itnRetained = q => 1 - 0.46 * (1 - q);

  M.applyItn = function (m, a, g, cov, susceptibility = 1,
                         kill = 0.5, feed = 0.3, mort = 0.3) {
    const n = cov * M.itnRetained(susceptibility);
    return { m: m * (1 - n * kill), a: a * (1 - n * feed), g: g * (1 + n * mort) };
  };

  M.eipSurvival = (g, n, stages = 4) => Math.pow(1 + g * n / stages, -stages);

  // Linear interpolation held flat beyond the ends (R's approxfun, rule = 2)
  M.interp = function (xs, ys, x) {
    if (x <= xs[0]) return ys[0];
    const n = xs.length;
    if (x >= xs[n - 1]) return ys[n - 1];
    let i = 1;
    while (xs[i] < x) i++;
    const w = (x - xs[i - 1]) / (xs[i] - xs[i - 1]);
    return ys[i - 1] + w * (ys[i] - ys[i - 1]);
  };

  // Ross-Macdonald solved with RK4 on a quarter-day step; output at `times`.
  // eipDays = null: the 2-state model. Otherwise newly infected mosquitoes pass
  // through `stages` exposed classes before becoming infectious.
  // Literature constants, as TRANSMISSION in R/epiwave-foi-model.R
  M.TRANSMISSION = { b: 0.5, c: 0.5, r: 1 / 180, eipDays: 10 };

  // Steady state with parameters held constant (rm_equilibrium() in R)
  M.rmEquilibrium = function (m, a, g, b, c, r, eipDays = null, stages = 4) {
    const S = eipDays === null ? 1 : M.eipSurvival(g, eipDays, stages);
    const zOf = x => S * a * c * x / (a * c * x + g);
    const k = eipDays === null ? 0 : stages;
    if (m * a * a * b * c * S / (g * r) <= 1) return { x: 0, z: 0, y: new Array(k).fill(0) };
    let lo = 1e-12, hi = 1 - 1e-12;
    const f = x => m * a * b * zOf(x) * (1 - x) - r * x;
    for (let i = 0; i < 200; i++) { const mid = (lo + hi) / 2; if (f(mid) > 0) lo = mid; else hi = mid; }
    const x = (lo + hi) / 2, inflow = a * c * x * g / (g + a * c * x), leave = k ? stages / eipDays : 0;
    const y = Array.from({ length: k }, (_, i) => inflow * Math.pow(leave / (leave + g), i) / (leave + g));
    return { x, z: zOf(x), y };
  };

  M.solveRM = function ({ times, m, a, g, b = M.TRANSMISSION.b, c = M.TRANSMISSION.c, r = M.TRANSMISSION.r,
                          x0 = 0.01, z0 = 0.001, h = 0.25, eipDays = null, stages = 4, start = null }) {
    const k = eipDays === null ? 0 : stages;
    const f = (t, s) => {
      const mt = M.interp(times, m, t), at = M.interp(times, a, t), gt = M.interp(times, g, t);
      const x = s[0], z = s[k + 1];
      const dx = mt * at * b * z * (1 - x) - r * x;
      if (k === 0) return [dx, at * c * x * (1 - z) - gt * z];
      const y = s.slice(1, k + 1), leave = k / eipDays, ySum = y.reduce((p, v) => p + v, 0);
      const dy = y.map((v, i) => (i === 0 ? at * c * x * (1 - ySum - z) : leave * y[i - 1]) - (gt + leave) * v);
      return [dx, ...dy, leave * y[k - 1] - gt * z];
    };
    const st = start === "equilibrium" ? M.rmEquilibrium(m[0], a[0], g[0], b, c, r, eipDays, stages) : null;
    const xStart = st ? st.x : x0, zStart = st ? st.z : z0;
    const y0 = !k ? [] : st ? st.y : new Array(k).fill(Number.isFinite(eipDays) ? g[0] * eipDays * z0 / k : 0);
    let s = [xStart, ...y0, zStart], t = times[0];
    const xs = [xStart], zs = [zStart];
    const add = (u, v, w) => u.map((ui, i) => ui + w * v[i]);
    for (let j = 1; j < times.length; j++) {
      while (t < times[j] - 1e-9) {
        const step = Math.min(h, times[j] - t);
        const k1 = f(t, s), k2 = f(t + step / 2, add(s, k1, step / 2));
        const k3 = f(t + step / 2, add(s, k2, step / 2)), k4 = f(t + step, add(s, k3, step));
        s = s.map((v, i) => v + step / 6 * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]));
        t += step;
      }
      xs.push(s[0]); zs.push(s[k + 1]);
    }
    const Istar = zs.map((zz, i) => m[i] * a[i] * b * zz);
    return { x: xs, z: zs, Istar };
  };

  // The demo site of simulate_epiwave_data(): seasonal m, ITN ramp to `itnMax`
  M.stage1Site = function ({ nTimes = 48, baselineM = 0.2, amplitude = 0.6, baselineA = 0.3,
                             baselineG = 0.1, itnMax = 0.7, susceptibility = 0.8, b = M.TRANSMISSION.b,
                             c = M.TRANSMISSION.c, r = M.TRANSMISSION.r, eipDays = M.TRANSMISSION.eipDays }) {
    const times = Array.from({ length: nTimes + 1 }, (_, i) => i * 30);
    const m = [], a = [], g = [];
    times.forEach((t, i) => {
      const cov = itnMax * i / nTimes;
      const adj = M.applyItn(M.seasonalM(t, baselineM, amplitude), baselineA, baselineG, cov, susceptibility);
      m.push(adj.m); a.push(adj.a); g.push(adj.g);
    });
    return { times, m, a, g, ...M.solveRM({ times, m, a, g, b, c, r, eipDays, start: "equilibrium" }) };
  };

  M.R0 = (m, a, b, c, g, r, eipDays = null) => (m * a * a * b * c) / (g * r) * (eipDays === null ? 1 : M.eipSurvival(g, eipDays));
  // Equilibrium human prevalence: x* = (R0 - 1) / (R0 + a c / g)
  M.xStar = (m, a, b, c, g, r) => {
    const R0 = M.R0(m, a, b, c, g, r);
    return R0 <= 1 ? 0 : (R0 - 1) / (R0 + a * c / g);
  };
  M.zStar = (x, a, c, g) => (a * c * x) / (a * c * x + g);

  // ---- detectability kernel ----
  M.qDaily = function (d, maxDays = 30) {
    if (d < 0 || d > maxDays) return 0;
    const up = 1 / (1 + Math.exp(-(d * 2 - 5)));
    const down = 1 - 1 / (1 + Math.exp(-(d / 2 - 10)));
    return up * down;
  };

  // transform_convolution_kernel(): daily kernel -> weight per model step
  M.qMonthly = function (maxDays = 30, stepDays = 30) {
    const tp = Math.round(stepDays);
    const maxT = 1 + Math.ceil(maxDays / tp);
    let total = 0;
    for (let d = 0; d <= maxDays; d++) total += M.qDaily(d, maxDays);
    const sums = {};
    for (let inf = 1; inf <= maxT * tp; inf++) {
      for (let test = 1; test <= tp; test++) {
        const diff = inf - test;
        if (diff < 0) continue;
        const k = Math.floor((inf - 1) / tp);
        sums[k] = sums[k] || {};
        sums[k][test] = (sums[k][test] || 0) + M.qDaily(diff, maxDays);
      }
    }
    const w = [];
    for (let k = 0; k <= maxT; k++) {
      if (!sums[k]) { w.push(0); continue; }
      const parts = Object.values(sums[k]);
      w.push(parts.reduce((s, v) => s + v / total, 0) / parts.length * total);
    }
    return w;
  };

  // ---- Stage 2 pieces ----
  M.matern52 = (d, phi) => { const r = Math.sqrt(5) * d / phi; return (1 + r + r * r / 3) * Math.exp(-r); };

  M.ar1 = function (rho, innov) {           // innov: [sites][times]
    return innov.map(row => {
      const out = [];
      row.forEach((e, t) => out.push(t === 0 ? e : rho * out[t - 1] + e));
      return out;
    });
  };

  M.cholesky = function (A) {
    const n = A.length, L = A.map(() => new Array(n).fill(0));
    for (let i = 0; i < n; i++) {
      for (let j = 0; j <= i; j++) {
        let s = A[i][j];
        for (let k = 0; k < j; k++) s -= L[i][k] * L[j][k];
        L[i][j] = i === j ? Math.sqrt(Math.max(s, 1e-12)) : s / L[j][j];
      }
    }
    return L;
  };

  M.dist = (p, q) => Math.hypot(p[0] - q[0], p[1] - q[1]);

  // Simulate epsilon: Matern 5/2 in space, AR(1) in time (as simulate_gp_residuals)
  M.simulateField = function (coords, nTimes, sigma, phi, rho, seed = 1) {
    const n = coords.length, R = M.rng(seed);
    const K = coords.map(p => coords.map(q => M.matern52(M.dist(p, q), phi)));
    for (let i = 0; i < n; i++) K[i][i] += 1e-6;
    const L = M.cholesky(K);
    const eps = coords.map(() => new Array(nTimes).fill(0));
    for (let t = 0; t < nTimes; t++) {
      const z = Array.from({ length: n }, () => R.norm());
      for (let i = 0; i < n; i++) {
        let v = 0;
        for (let k = 0; k <= i; k++) v += L[i][k] * z[k];
        eps[i][t] = (t === 0 ? 0 : rho * eps[i][t - 1]) + sigma * v;
      }
    }
    return eps;
  };

  // Solve A x = b for symmetric positive-definite A via Cholesky
  M.solveSPD = function (A, b) {
    const L = M.cholesky(A), n = b.length, y = new Array(n), x = new Array(n);
    for (let i = 0; i < n; i++) { let s = b[i]; for (let k = 0; k < i; k++) s -= L[i][k] * y[k]; y[i] = s / L[i][i]; }
    for (let i = n - 1; i >= 0; i--) { let s = y[i]; for (let k = i + 1; k < n; k++) s -= L[k][i] * x[k]; x[i] = s / L[i][i]; }
    return x;
  };

  // Gaussian-process prediction (kriging with a given mean) at held-out sites, one month
  // at a time: y_obs = residual at observed sites, returns predicted residual.
  M.krige = function (coordsObs, yObs, coordsNew, tau2, phi, noise) {
    const K = coordsObs.map((p, i) => coordsObs.map((q, j) => tau2 * M.matern52(M.dist(p, q), phi) + (i === j ? noise : 0)));
    const alpha = M.solveSPD(K, yObs);
    return coordsNew.map(p => coordsObs.reduce((s, q, j) => s + tau2 * M.matern52(M.dist(p, q), phi) * alpha[j], 0));
  };

  // Log-likelihood of replicate 1's data at (alpha, gamma), with the true residual field held fixed.
  M.makeLogLik = function (R1) {
    const S = R1.I_star.length, T = R1.I_star[0].length, w = M.qMonthly();
    const base = R1.I_star.map((row, s) => row.map((v, t) => Math.log(Math.max(v, 1e-6)) + R1.epsilon[s][t]));
    const sv = R1.surveys, idx = sv.survey_indices.map(i => i - 1); // R's column-major, 1-based
    const N = R1.population;
    return function (alpha, gamma) {
      const I = base.map(row => row.map(v => Math.exp(alpha + v)));
      let lc = 0;
      for (let s = 0; s < S; s++) for (let t = 0; t < T; t++) { const mu = gamma * I[s][t] * N, c = R1.cases[s][t]; lc += c * Math.log(mu) - mu; }
      let lp = 0;
      idx.forEach((k, j) => {
        const s = k % S, t = Math.floor(k / S);
        let L = 0; for (let lag = 0; lag < w.length && lag <= t; lag++) L += I[s][t - lag] * w[lag];
        const p = Math.min(1 - Math.exp(-L), 1 - 1e-12), y = sv.n_positive[j], n = sv.n_tested[j];
        lp += y * Math.log(p) + (n - y) * Math.log(1 - p);
      });
      return { cases: lc, surveys: lp };
    };
  };

  M.mean = xs => xs.reduce((s, v) => s + v, 0) / xs.length;
  M.median = xs => { const s = [...xs].sort((a, b) => a - b), n = s.length; return n % 2 ? s[(n - 1) / 2] : (s[n / 2 - 1] + s[n / 2]) / 2; };

  if (typeof module !== "undefined") module.exports = M; else root.EpiModel = M;
})(typeof window !== "undefined" ? window : globalThis);
