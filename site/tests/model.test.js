// Checks the JavaScript port against values exported from R (site/export_data.R).
//   node site/tests/model.test.js
const fs = require("fs");
const path = require("path");
const M = require("../src/model.js");

const src = fs.readFileSync(path.join(__dirname, "../src/data.js"), "utf8");
const window = {};
eval(src);
const { refs, replicate1 } = window.EPIWAVE_DATA;

let failed = 0;
function check(name, got, want, relTol) {
  const g = [].concat(got).flat(), w = [].concat(want).flat();
  let worst = 0;
  w.forEach((v, i) => { worst = Math.max(worst, Math.abs(g[i] - v) / Math.max(Math.abs(v), 1e-12)); });
  const ok = g.length === w.length && worst <= relTol;
  if (!ok) failed++;
  console.log(`${ok ? "PASS" : "FAIL"}  ${name}  (max relative error ${worst.toExponential(2)}, tol ${relTol})`);
}

check("q_daily", refs.q_daily.map((_, d) => M.qDaily(d)), refs.q_daily, 1e-12);
check("q_monthly", M.qMonthly(), refs.q_monthly.slice(0, 3), 1e-12);
check("matern52", refs.matern.d.map(d => M.matern52(d, refs.matern.phi)), refs.matern.k, 1e-12);
check("ar1", M.ar1(refs.ar1.rho, refs.ar1.innovations), refs.ar1.out, 1e-12);

const site = M.stage1Site({});
check("Ross-Macdonald x (vs deSolve lsoda)", site.x, refs.ode_site1.x, 1e-4);
check("I* (vs deSolve lsoda)", site.Istar, refs.ode_site1.I_star, 1e-4);
check("I* without the EIP (vs deSolve lsoda)", M.stage1Site({ eipDays: null }).Istar, refs.ode_no_eip.I_star, 1e-3);
check("equilibrium prevalence before nets is 0.481", [M.rmEquilibrium(0.2, 0.3, 0.1, 0.5, 0.5, 1 / 180, 10).x], [0.4809], 1e-3);
check("ITN effect retained at 0.8 susceptibility", [M.itnRetained(0.8)], [0.908], 1e-12);
check("I* is identical at every site (finding F2)",
      replicate1.I_star.map(row => row[20]), replicate1.I_star.map(() => replicate1.I_star[0][20]), 1e-12);

const meanLog = M.mean(replicate1.I_star.flat().map(v => Math.log(Math.max(v, 1e-6))));
check("mean log I*", [meanLog], [replicate1.mean_log_I_star], 1e-6);

// Likelihood lab: survey prevalence at the truth must equal R's true_prevalence
{
  const S = replicate1.I_star.length, w = M.qMonthly(), sv = replicate1.surveys;
  const I = replicate1.I_star.map((row, s) => row.map((v, t) => Math.exp(Math.log(Math.max(v, 1e-6)) + replicate1.epsilon[s][t])));
  const p = sv.survey_indices.map(k => { k -= 1; const s = k % S, t = Math.floor(k / S); let L = 0; for (let l = 0; l < w.length && l <= t; l++) L += I[s][t - l] * w[l]; return Math.max(1 - Math.exp(-L), 1e-6); });
  check("survey prevalence at the truth (vs R)", p, sv.true_prevalence, 1e-4);
  const ll = M.makeLogLik(replicate1);
  const a = ll(0, 0.1), b = ll(0.5, 0.1 * Math.exp(-0.5));
  check("cases see only gamma * exp(alpha)", [b.cases], [a.cases], 1e-9);
  let best = { v: -Infinity };
  for (let al = -0.6; al <= 0.6; al += 0.02) for (let lg = Math.log(0.05); lg <= Math.log(0.2); lg += 0.02) {
    const r = ll(al, Math.exp(lg)), v = r.cases + r.surveys; if (v > best.v) best = { v, al, g: Math.exp(lg) };
  }
  const ok = Math.abs(best.al) < 0.1 && Math.abs(best.g / 0.1 - 1) < 0.1;
  if (!ok) failed++;
  console.log(`${ok ? "PASS" : "FAIL"}  both likelihoods peak at the truth (alpha ${best.al.toFixed(2)}, gamma ${best.g.toFixed(3)})`);
}

if (failed) { console.log(`\n${failed} check(s) failed`); process.exit(1); }
console.log("\nAll checks passed.");
