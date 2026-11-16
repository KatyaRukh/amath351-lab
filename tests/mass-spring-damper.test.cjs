"use strict";
// tests/mass-spring-damper.test.cjs
// Run with: node --test tests/mass-spring-damper.test.cjs

const { test } = require("node:test");
const assert   = require("node:assert/strict");
const fs = require("node:fs");
const path = require("node:path");
const vm = require("node:vm");
const { solve, classify, exactSolution, ode, rk4Step } =
  require("../src/mass-spring-damper-model.cjs");

const tol = 1e-4; // tolerance for numerical vs analytic comparison

// ----- classify -----
test("classify: undamped when c = 0", () => {
  assert.equal(classify(1, 0, 4), "undamped");
});
test("classify: underdamped disc < 0", () => {
  assert.equal(classify(1, 1, 4), "underdamped");
});
test("classify: critical disc ≈ 0", () => {
  assert.equal(classify(1, 4, 4), "critical");
});
test("classify: overdamped disc > 0", () => {
  assert.equal(classify(1, 10, 4), "overdamped");
});

// ----- ODE residual -----
test("ode: returns correct derivatives", () => {
  // m=1, c=1, k=0.5; state=[0.1, 0]
  const [dy, dv] = ode(0, [0.1, 0], 1, 1, 0.5);
  assert.equal(dy, 0);                       // y' = v = 0
  assert.equal(dv, (-1*0 - 0.5*0.1) / 1);   // = -0.05
});

// ----- initial conditions satisfied -----
test("solve: initial values match y0 and yp0", () => {
  const { y, yp } = solve({ m: 1, c: 1, k: 0.5, y0: 0.1, yp0: 0 });
  assert.ok(Math.abs(y[0]  - 0.1) < 1e-14);
  assert.ok(Math.abs(yp[0] - 0.0) < 1e-14);
});

// ----- ODE residual throughout the trajectory -----
// At every step, the derivative y[i+1]-y[i] should be approximately h * ode(...)[0]
test("solve: trajectory satisfies ODE residual", () => {
  const { t, y, yp } = solve({ m: 1, c: 1, k: 0.5, y0: 0.1, yp0: 0, tEnd: 50, steps: 500 });
  for (let i = 0; i < t.length - 1; i++) {
    const h = t[i + 1] - t[i];
    const [dydt, dvdt] = ode(t[i], [y[i], yp[i]], 1, 1, 0.5);
    // Forward difference should be close to the derivative (first-order check)
    const fwdY  = (y[i + 1]  - y[i])  / h;
    const fwdYp = (yp[i + 1] - yp[i]) / h;
    assert.ok(Math.abs(fwdY  - dydt) < 0.05, `displacement residual at i=${i}: ${fwdY} vs ${dydt}`);
    assert.ok(Math.abs(fwdYp - dvdt) < 0.05, `velocity residual at i=${i}: ${fwdYp} vs ${dvdt}`);
  }
});

// ----- numerical vs analytic (underdamped) -----
test("solve: matches analytic solution for underdamped case (MATLAB defaults)", () => {
  const params = { m: 1, c: 1, k: 0.5, y0: 0.1, yp0: 0, tEnd: 50, steps: 2000 };
  const { t, y } = solve(params);
  const exact = exactSolution(1, 1, 0.5, 0.1, 0, t);
  for (let i = 0; i < t.length; i++) {
    assert.ok(Math.abs(y[i] - exact[i]) < tol,
      `t=${t[i].toFixed(2)}: numerical=${y[i].toFixed(6)}, exact=${exact[i].toFixed(6)}`);
  }
});

// ----- numerical vs analytic (undamped) -----
test("solve: matches analytic for undamped oscillator", () => {
  const { t, y } = solve({ m: 1, c: 0, k: 4, y0: 1, yp0: 0, tEnd: 20, steps: 2000 });
  const exact = exactSolution(1, 0, 4, 1, 0, t);
  for (let i = 0; i < t.length; i++) {
    assert.ok(Math.abs(y[i] - exact[i]) < tol,
      `t=${t[i].toFixed(2)}: numerical=${y[i].toFixed(6)}, exact=${exact[i].toFixed(6)}`);
  }
});

// ----- energy is conserved when c = 0 -----
test("solve: energy conserved for undamped oscillator", () => {
  const m = 2, k = 3;
  const { t, y, yp } = solve({ m, c: 0, k, y0: 1, yp0: 0, tEnd: 30, steps: 3000 });
  const E0 = 0.5 * m * yp[0] ** 2 + 0.5 * k * y[0] ** 2;
  for (let i = 0; i < t.length; i++) {
    const E = 0.5 * m * yp[i] ** 2 + 0.5 * k * y[i] ** 2;
    assert.ok(Math.abs(E - E0) < 1e-3, `Energy drift at t=${t[i].toFixed(2)}: ${E} vs ${E0}`);
  }
});

// ----- energy dissipates when c > 0 -----
test("solve: energy decreases monotonically for overdamped case", () => {
  const m = 1, c = 10, k = 4;
  const { y, yp } = solve({ m, c, k, y0: 1, yp0: 0, tEnd: 10, steps: 1000 });
  const energy = (i) => 0.5 * m * yp[i] ** 2 + 0.5 * k * y[i] ** 2;
  for (let i = 1; i < y.length; i++) {
    assert.ok(energy(i) <= energy(i - 1) + 1e-6,
      `Energy increased at step ${i}`);
  }
});

// ----- parameter validation -----
test("solve: throws for m <= 0", () => {
  assert.throws(() => solve({ m: 0, c: 1, k: 1, y0: 0, yp0: 0 }), /m must be positive/);
});
test("solve: throws for k <= 0", () => {
  assert.throws(() => solve({ m: 1, c: 1, k: -1, y0: 0, yp0: 0 }), /k must be positive/);
});
test("solve: throws for c < 0", () => {
  assert.throws(() => solve({ m: 1, c: -1, k: 1, y0: 0, yp0: 0 }), /c must be non-negative/);
});

test("the shared model loads in a browser without CommonJS globals", () => {
  const context = vm.createContext({});
  const source = fs.readFileSync(path.join(__dirname, "../src/mass-spring-damper-model.cjs"), "utf8");
  vm.runInContext(source, context);
  const params = { m: 1, c: 4, k: 4, y0: 1, yp0: 0 };
  const browser = context.Amath351MassSpring;
  assert.equal(browser.classify(1, 4, 4), "critical");
  assert.deepEqual(Array.from(browser.solve(params).y), solve(params).y);
});

test("the page embeds the shared model once and generates valid browser JavaScript", () => {
  const source = fs.readFileSync(path.join(__dirname, "../src/mass-spring-damper.html"), "utf8");
  const model = fs.readFileSync(path.join(__dirname, "../src/mass-spring-damper-model.cjs"), "utf8");
  assert.equal(source.split("{{MODEL}}").length - 1, 1);
  const scripts = [...source.replaceAll("{{MODEL}}", model).matchAll(/<script>([\s\S]*?)<\/script>/g)];
  assert.equal(scripts.length, 1);
  assert.doesNotThrow(() => new vm.Script(scripts[0][1]));
  assert.match(source, /const \{ solve, classify \} = globalThis\.Amath351MassSpring/);
});

test("the damping control includes the overdamped preset", () => {
  const source = fs.readFileSync(path.join(__dirname, "../src/mass-spring-damper.html"), "utf8");
  const input = source.match(/<input\b[^>]*id="ctrl-c"[^>]*>/)[0];
  const max = Number(input.match(/\bmax="([\d.]+)"/)[1]);
  const preset = Number(source.match(/overdamped:\s*\{[^}]*\bc:\s*([\d.]+)/)[1]);
  assert.ok(preset <= max);
  assert.equal(Number(input.match(/aria-valuemax="([\d.]+)"/)[1]), max);
});
