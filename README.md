# [WIP] ode-rk45-cash-karp [![Build Status](https://travis-ci.org/scijs/ode-rk45-cash-karp.svg)](https://travis-ci.org/scijs/ode-rk45-cash-karp) [![npm version](https://badge.fury.io/js/ode-rk45-cash-karp.svg)](http://badge.fury.io/js/ode-rk45-cash-karp) [![Dependency Status](https://david-dm.org/scijs/ode-rk45-cash-karp.svg)](https://david-dm.org/scijs/ode-rk45-cash-karp)

> Integrate a system of ODEs using the Fifth Order Adaptive Cash-Karp (RKCK) method


## Introduction

This module integrates a system of ordinary differential equations of the form

<p align="center"><img alt="&bsol;begin&lcub;eqnarray&midast;&rcub; y&apos;&lpar;t&rpar; &amp;&equals;&amp; f&lpar;t&comma; y&lpar;t&rpar;&rpar;&comma; &bsol;&bsol; y&lpar;t&lowbar;0&rpar; &amp;&equals;&amp; y&lowbar;0 &bsol;end&lcub;eqnarray&midast;&rcub;" valign="middle" src="docs/images/begineqnarray-yt-ft-yt-yt_0-y_0-endeqnarray-0298eae3db.png" width="187" height="61"></p>

where <img alt="y" valign="middle" src="docs/images/y-720f311276.png" width="14.5" height="20"> is a vector of length <img alt="n" valign="middle" src="docs/images/n-9baedbc330.png" width="16" height="16">. Given time step <img alt="&bsol;Delta t" valign="middle" src="docs/images/delta-t-a20a5fe4f2.png" width="28" height="16">, the Cash-Karp method uses a fifth order Runge-Kutta scheme with a fourth order embedded estimator in order to compute an error estimate. In other words, the same coefficients used in calculating the fifth order update can be used to calculate a fourth order estimate. The difference yields an error estimate, and the error estimate controls the timestep <img alt="&bsol;Delta t" valign="middle" src="docs/images/delta-t-a20a5fe4f2.png" width="28" height="16">.

## Install

```bash
$ npm install ode-rk45-cash-karp
```

## Example

```javascript
var rk45 = require('ode-rk45-cash-karp')

var mu = 4.0

var deriv = function(dydt, y, t) {
  dydt[0] = y[1]
  dydt[1] = mu * (1-y[0]*y[0])*y[1] - y[0]
}

var y0 = [2,0]
var n = 1000
var t0 = 0
var dt = 1

var integrator = rk45( y0, deriv, t0, dt )

// Store the output in simple arrays:
var t = []
var x = []
var y = []

while( integrator.step( tmax ) ) {
  integrator.steps(n)
  t.push( integrator.t )
  x.push( integrator.y[0] )
  y.push( integrator.y[1] )
}
```



## API

### `require('ode-rk45-cash-karp')( y0, deriv, t0, dt0 [, options] )`
#### Arguments:
- `y0`: an array or typed array containing initial conditions. This vector is updated in-place with each integrator step.
- `deriv`: a function that calculates the derivative. Format is `function( dydt, y, t )`. Inputs are current state `y` and current time `t`, output is the calculated derivative `dydt`.
- `t0`: initial time <img alt="t" valign="middle" src="docs/images/t-fc93da6f4d.png" width="11.5" height="16">.
- `dt0`: initial time step <img alt="&bsol;Delta t" valign="middle" src="docs/images/delta-t-a20a5fe4f2.png" width="28" height="16">.
- `options`: an optional associative array of options. Valid parameters are:
  - `tol`: The target error level to be acheived. Default is: <img alt="1 &bsol;cdot 10&Hat;&lcub;-8&rcub;" valign="middle" src="docs/images/1-cdot-10-8-9ceae48083.png" width="70" height="20">.
  - `maxIncreaseFactor`: The maximum factor by which to increase the timestep if the error tolerance is met. Default value is 10. This limit is applied at the end of a successful timestep.
  - `maxDecreaseFactor`: The maximum factor by which to decrease the timestep if the error tolerance is not met. Default value is 10. This limit is applied multiple times on each adaptive timestep until the error tolerance is acheived.
  - `dtMinMag`: The minimum allowed magnitude of <img alt="&bsol;Delta t" valign="middle" src="docs/images/delta-t-a20a5fe4f2.png" width="28" height="16">. If this limit is exceeded during adaptation, a warning will be printed to the console and the timestep will complete with <img alt="&bsol;Delta t" valign="middle" src="docs/images/delta-t-a20a5fe4f2.png" width="28" height="16"> clipped to the the prescribed magnitude. If `undefined`, this limit is ignored. Default value is `undefined`.
  - `dtMaxMag`: The maximum allowed magnitude of <img alt="&bsol;Delta t" valign="middle" src="docs/images/delta-t-a20a5fe4f2.png" width="28" height="16">. This limit is applied at the beginning of each step. If a timestep larger than this magnitude is requested, the timestep will commence with <img alt="&bsol;Delta t" valign="middle" src="docs/images/delta-t-a20a5fe4f2.png" width="28" height="16"> clipped to the prescribed magnitude. If `undefined`, this limit is ignored. Default value is `undefined`.
  - `errorScaleFunction`: The function used to compute a normalizing factor for the error in a given dimension. See below for details.
  - `errorReduceFunction`: The reduce operation by which errors in each dimension are combined into a single error metric. See below for details.
  - `errorPostFunction`: An operation applied to the total error. For example, if using the <img alt="L&lowbar;2" valign="middle" src="docs/images/l_2-23fd536b11.png" width="26.5" height="19"> norm this would be a square root. See below for details.


#### Returns:
Initialized integrator object.

**Properties:**
- `n`: dimension of `y0`.
- `y`: current state. Initialized as a shallow copy of input `y0`.
- `deriv`: function that calculates the derivative. Initialized from input. May be changed.
- `t`: current time, incremented by `dt` with each time step.
- `dt`: current time step <img alt="&bsol;Delta t" valign="middle" src="docs/images/delta-t-a20a5fe4f2.png" width="28" height="16">. Initialized from input `dt0`. May be changed, but will be overwritten with each adaptive step in order to acheive the desired error bound.
- `maxDecreaseFactor`: set from options or defaults. May be changed.
- `maxIncreaseFactor`: set from options or defaults. May be changed.
- `dtMinMag`: set from options or defaults. May be changed.
- `dtMaxMag`: set from options or defaults. May be changed.

#### Methods:
- `.step( [tLimit] )`: takes a single step of the integrator and stores the result in-place in the `y` property. Returns true if `tLimit` was not provided or if `t` has not reached the limit, otherwise returns false, meaning `t` has reached `tLimit`.
- `.steps( n, [tLimit] )`: takes `n` steps of the integrator, storing the result in-place in the `y` property. Exits early if `tLimit` is reached. Returns true if `tLimit` was not provided or if `t` has not reached the limit, otherwise returns false, meaning `t` has reached `tLimit`. Note that `.steps( Infinity, 10 )` is valid and will take whatever number of step is required to reach <img alt="t&equals;10" valign="middle" src="docs/images/t10-b3a7f6176e.png" width="58" height="16">.

### Error Estimation
Ideally, there would be no choices in error computation since this library would implement the best possible choices, but I've left this configurable.

#### `errorScaleFunction: function( i, dt, y, dydt )`
This function receives the dimension number `i`, the current timestep `dt`, the current state `y`, and the derivative calculated at the beginning of the step, `dydt`. It must return a normalization factor by which the error in the given dimension is normalized. By default, this is:

```
function errorScaleFunction( i, dt, y, dydt ) {
  return Math.abs(y) + Math.abs(dt * dydt) + 1e-32
}
```

#### `errorReduceFunction: function( i, accumulatedError, errorEstimate )`
This function performs a reduce operation on the per-dimension error. `accumulatedError` is initially zero. The function must add the error estimate in dimension `i` and return a new error estimate. By default, the error reduce function simply returns the maximum error:

```
function errorReduceFunction( i, accumulatedError, errorEstimate ) {
  return Math.max( accumulatedError, Math.abs(errorEstimate))
}
```

#### `errorPostFunction: function( i, accumulatedError, errorEstimate )`
This function applies a mapping to the total reduced error resulting from `errorReduceFunction`. For the <img alt="L&lowbar;2" valign="middle" src="docs/images/l_2-23fd536b11.png" width="26.5" height="19"> norm, this would just be `Math.sqrt`; for the <img alt="L&lowbar;&bsol;infty" valign="middle" src="docs/images/l_infty-c904452e37.png" width="34.5" height="19"> norm, this is simply a no-op:

```
function errorPostFunction( accumulatedError ) {
  return accumulatedError
}
```



## Credits

(c) 2015 Ricky Reusser. MIT License