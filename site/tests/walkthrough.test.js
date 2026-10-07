// Keeps the code walkthrough honest:
//  1. every code line on the page exists in R/epiwave-foi-model.R, inside the
//     function(s) its section names (comment lines do not count);
//  2. every assignment in fit_epiwave_gp() is explained on the page, so a new
//     model line cannot go unexplained.
//   node site/tests/walkthrough.test.js
const fs = require("fs");
const path = require("path");

const rLines = fs.readFileSync(path.join(__dirname, "../../R/epiwave-foi-model.R"), "utf8").split(/\r?\n/);
const window = {};
eval(fs.readFileSync(path.join(__dirname, "../src/walkthrough.js"), "utf8"));
const sections = window.EPIWAVE_WALKTHROUGH;

// split the file into function bodies; everything else is "(top level)"
const scopes = { "(top level)": [] };
let current = null;
for (const line of rLines) {
  const start = line.match(/^([A-Za-z_][A-Za-z0-9_.]*) <- function\(/);
  if (!current && start) { current = start[1]; scopes[current] = [line.trim()]; continue; }
  if (current) {
    scopes[current].push(line.trim());
    if (line === "}") current = null;
  } else {
    scopes["(top level)"].push(line.trim());
  }
}
const isCode = l => l.length && !l.startsWith("#");

let failed = 0;
for (const s of sections) {
  const pool = new Set(s.fns.flatMap(f => (scopes[f] || []).filter(isCode)));
  if (s.fns.some(f => !scopes[f])) { failed++; console.log(`FAIL  section "${s.id}" names an unknown function: ${s.fns.join(", ")}`); }
  for (const row of s.rows) {
    for (const codeLine of row.c.split("\n").map(l => l.trim())) {
      if (!pool.has(codeLine)) { failed++; console.log(`FAIL  [${s.id}] not found in ${s.fns.join("/")}: ${codeLine}`); }
    }
  }
}

const fitRows = new Set(sections.find(s => s.id === "fit").rows.flatMap(r => r.c.split("\n").map(l => l.trim())));
for (const line of scopes.fit_epiwave_gp.filter(isCode)) {
  if (line.includes("<-") && !line.startsWith("fit_epiwave_gp") && !fitRows.has(line)) {
    failed++; console.log(`FAIL  fit_epiwave_gp line not explained on the page: ${line}`);
  }
}

const nRows = sections.reduce((n, s) => n + s.rows.length, 0);
if (failed) { console.log(`\n${failed} problem(s)`); process.exit(1); }
console.log(`PASS  all ${nRows} walkthrough rows match the R code; every fit_epiwave_gp assignment is explained.`);
