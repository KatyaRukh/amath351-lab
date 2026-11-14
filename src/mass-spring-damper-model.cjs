"use strict";
// mass-spring-damper-model.cjs
// Numerical model for the second-order mass-spring-damper ODE:
//   m y'' + c y' + k y = 0
//
// Reduced to a first-order system (as in MATLAB's ode45 workflow):
//   Let y1 = y,  y2 = y'
//   y1' = y2
//   y2' = (-c*y2 - k*y1) / m

// ----- ODE right-hand side -----
function ode(t, y, m, c, k) {
  // y = [displacement, velocity]
  return [
    y[1],
    (-c * y[1] - k * y[0]) / m
  ];
}

// ----- Single RK4 step -----
function rk4Step(t, y, h, m, c, k) {
  const f = (tt, yy) => ode(tt, yy, m, c, k);
  const add = (a, b, s) => [a[0] + s * b[0], a[1] + s * b[1]];

  const k1 = f(t,         y);
  const k2 = f(t + h / 2, add(y, k1, h / 2));
  const k3 = f(t + h / 2, add(y, k2, h / 2));
  const k4 = f(t + h,     add(y, k3, h));

  return [
    y[0] + (h / 6) * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0]),
    y[1] + (h / 6) * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1])
  ];
}

// ----- Full solver -----
// Returns arrays t, y (displacement), yp (velocity)
function solve({ m, c, k, y0, yp0, tEnd = 50, steps = 1000 }) {
  if (m <= 0) throw new Error("m must be positive");
  if (k <= 0) throw new Error("k must be positive");
  if (c < 0)  throw new Error("c must be non-negative");

  const h = tEnd / steps;
  const t  = new Array(steps + 1);
  const y  = new Array(steps + 1);
  const yp = new Array(steps + 1);

  let state = [y0, yp0];
  for (let i = 0; i <= steps; i++) {
    const ti = i * h;
    t[i]  = ti;
    y[i]  = state[0];
    yp[i] = state[1];
    if (i < steps) {
      state = rk4Step(ti, state, h, m, c, k);
    }
  }
  return { t, y, yp };
}

// ----- Damping classification -----
// Returns one of: 'undamped' | 'underdamped' | 'critical' | 'overdamped'
function classify(m, c, k) {
  if (c === 0) return "undamped";
  const disc = c * c - 4 * m * k;
  if (Math.abs(disc) < 1e-9) return "critical";
  return disc > 0 ? "overdamped" : "underdamped";
}

// ----- Exact analytic solution (for comparison / tests) -----
// Returns y(t) evaluated at array of times; null if no closed form is handy.
function exactSolution(m, c, k, y0, yp0, tArr) {
  const type = classify(m, c, k);
  const omega0 = Math.sqrt(k / m);

  if (type === "undamped") {
    // y = A cos(omega0 t) + B sin(omega0 t)
    const A = y0;
    const B = yp0 / omega0;
    return tArr.map(t => A * Math.cos(omega0 * t) + B * Math.sin(omega0 * t));
  }

  const alpha = c / (2 * m);

  if (type === "underdamped") {
    const omegaD = Math.sqrt(omega0 * omega0 - alpha * alpha);
    // y = e^(-alpha t) [A cos(omegaD t) + B sin(omegaD t)]
    const A = y0;
    const B = (yp0 + alpha * y0) / omegaD;
    return tArr.map(t =>
      Math.exp(-alpha * t) * (A * Math.cos(omegaD * t) + B * Math.sin(omegaD * t))
    );
  }

  if (type === "critical") {
    // y = e^(-alpha t) (A + B t)
    const A = y0;
    const B = yp0 + alpha * y0;
    return tArr.map(t => Math.exp(-alpha * t) * (A + B * t));
  }

  // overdamped
  const beta = Math.sqrt(alpha * alpha - omega0 * omega0);
  const r1 = -alpha + beta;
  const r2 = -alpha - beta;
  // y = A e^(r1 t) + B e^(r2 t)
  // y(0) = A + B = y0
  // y'(0) = A r1 + B r2 = yp0
  const A = (yp0 - r2 * y0) / (r1 - r2);
  const B = y0 - A;
  return tArr.map(t => A * Math.exp(r1 * t) + B * Math.exp(r2 * t));
}

module.exports = { solve, classify, exactSolution, ode, rk4Step };
