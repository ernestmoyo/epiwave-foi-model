// Assembles site/src/* into one self-contained page: site/dist/epiwave-study.html
//   node site/build.js
const fs = require("fs");
const path = require("path");
const src = f => fs.readFileSync(path.join(__dirname, "src", f), "utf8");

const page = src("shell.html")
  .replace("/*STYLE*/", () => src("style.css"))
  .replace("/*DATA*/", () => src("data.js"))
  .replace("/*CONTENT*/", () => src("content.js"))
  .replace("/*WALK*/", () => src("walkthrough.js"))
  .replace("/*MODEL*/", () => src("model.js"))
  .replace("/*APP*/", () => src("app.js"));

fs.mkdirSync(path.join(__dirname, "dist"), { recursive: true });
const out = path.join(__dirname, "dist", "epiwave-study.html");
fs.writeFileSync(out, page);
// Vercel deploys site/deploy/epiwave-study (the folder name is the project name)
const deployDir = path.join(__dirname, "deploy", "epiwave-study");
fs.mkdirSync(deployDir, { recursive: true });
const head = '<!doctype html>\n<meta charset="utf-8">\n' +
  '<meta name="viewport" content="width=device-width, initial-scale=1, viewport-fit=cover">\n';
fs.writeFileSync(path.join(deployDir, "index.html"), head + page);
console.log(`wrote ${out} (${(page.length / 1024).toFixed(0)} KB) and the Vercel copy in ${deployDir}`);
