'use strict'

module.exports = IntegratorFactory

function defaultErrorScaleFunction( i, dt, y, dydt ) {
  return Math.abs(y) + Math.abs(dt * dydt) + 1e-32
}

function defaultErrorReduceFunction( i, accumulatedError, errorEstimate ) {
  return Math.max( accumulatedError, Math.abs(errorEstimate))
}

function defaultErrorPostFunction( accumulatedError ) {
  return accumulatedError
}

function minMag (a, b) {
  return (a > 0 ? Math.min : Math.max)(a, b);
}

function maxMag (a, b) {
  return (a > 0 ? Math.max : Math.min)(a, b);
}

var Integrator = function Integrator( y0, deriv, t0, dt0, options ) {
  var opts = options || {}
  this.tol = opts.tol===undefined ? 1e-8 : opts.tol
  this.maxIncreaseFactor = opts.maxIncreaseFactor===undefined ? 10 : opts.maxIncreaseFactor
  this.maxDecreaseFactor = opts.maxDecreaseFactor===undefined ? 10 : opts.maxDecreaseFactor
  this.dtMinMag = opts.dtMinMag===undefined ? 0 : Math.abs(opts.dtMinMag)
  this.dtMaxMag = opts.dtMaxMag===undefined ? undefined : Math.abs(opts.dtMaxMag)
  this.verbose = opts.verbose===undefined ? true : !!opts.verbose;

  var logCnt = 0
  var maxLogs = 10
  var maxLogWarningIssued = false
  this.__log = function (method, msg) {
    if (!this.verbose) return;
    if (logCnt < maxLogs) {
      console.log('ode45-cash-karp::' + method + '(): ' + msg)
      logCnt++
    } else {
      if (!maxLogWarningIssued) {
        console.log('ode45-cash-karp: too many warnings. Silencing further output')
        maxLogWarningIssued = true
      }
    }
  }.bind(this)

  this.errorScaleFunction = opts.errorScaleFunction === undefined ? defaultErrorScaleFunction : opts.errorScaleFunction
  this.errorReduceFunction = opts.errorReduceFunction === undefined ? defaultErrorReduceFunction : opts.errorReduceFunction
  this.errorPostFunction = opts.errorPostFunction === undefined ? defaultErrorPostFunction : opts.errorPostFunction

  // This is technically a parameter, but I think the value of leaving this undocumented exceeds the
  // value of documenting this and only adding confusion. I can't imagine this will even need to be
  // modified.
  this.safetyFactor = opts.safetyFactor===undefined ? 0.9 : opts.safetyFactor

  // Bind variables to this:
  this.deriv = deriv
  this.y = y0
  this.n = this.y.length
  this.dt = dt0
  this.t = t0

  // Create a scratch array into which we compute the derivative:
  this._ctor = this.y.constructor

  this._errorScale = new this._ctor( this.n )
  this._w = new this._ctor( this.n )
  this._k1 = new this._ctor( this.n )
  this._k2 = new this._ctor( this.n )
  this._k3 = new this._ctor( this.n )
  this._k4 = new this._ctor( this.n )
  this._k5 = new this._ctor( this.n )
  this._k6 = new this._ctor( this.n )
}

Integrator.prototype._calculateK1 = function() {
  this.deriv( this._k1, this.y, this.t )

  return this
}

Integrator.prototype._calculateKs = function(dt) {
  var i

  //var a21 =  0.200000000000000000 // 1/5
  //var a31 =  0.075000000000000000 // 3/40
  //var a32 =  0.225000000000000000 // 9/40
  //var a41 =  0.300000000000000000 // 3/10
  //var a42 = -0.900000000000000000 // -9/10
  //var a43 =  1.200000000000000000 // 6/5
  //var a51 = -0.203703703703703703 // -11/54
  //var a52 =  2.500000000000000000 // 5/2
  //var a53 = -2.592592592592592592 // -70/27
  //var a54 =  1.296296296296296296 // 35/27
  //var a61 =  0.029495804398148148 // 1631/55296
  //var a62 =  0.341796875000000000 // 175/512
  //var a63 =  0.041594328703703703 // 575/13824
  //var a64 =  0.400345413773148148 // 44275/110592
  //var a65 =  0.061767578125000000 // 253/4096

  //var b1  =  0.000000000000000000 // 0
  //var b2  =  0.200000000000000000 // 1/5
  //var b3  =  0.300000000000000000 // 3/10
  //var b4  =  0.600000000000000000 // 3/5
  //var b5  =  1.000000000000000000 // 1
  //var b6  =  0.875000000000000000 // 7/8

  // Same for every step, so don't repeat:
  //this.deriv( this._k1, this.y, this.t )

  for(i=0; i<this.n; i++) {
    this._w[i] = this.y[i] + dt * (
      0.2 * this._k1[i] )
  }

  this.deriv( this._k2, this._w, this.t + dt * 0.2 )

  for(i=0; i<this.n; i++) {
    this._w[i] = this.y[i] + dt * (
      0.075 * this._k1[i] +
      0.225 * this._k2[i] )
  }

  this.deriv( this._k3, this._w, this.t + dt * 0.3 )

  for(i=0; i<this.n; i++) {
    this._w[i] = this.y[i] + dt * (
       0.3 * this._k1[i] +
      -0.9 * this._k2[i] +
       1.2 * this._k3[i] )
  }

  this.deriv( this._k4, this._w, this.t + dt * 0.6 )

  for(i=0; i<this.n; i++) {
    this._w[i] = this.y[i] + dt * (
      -0.203703703703703703 * this._k1[i] +
       2.5                  * this._k2[i] +
      -2.592592592592592592 * this._k3[i] +
       1.296296296296296296 * this._k4[i] )
  }

  this.deriv( this._k5, this._w, this.t + dt /* * b5 */ )

  for(i=0; i<this.n; i++) {
    this._w[i] = this.y[i] + dt * (
      0.029495804398148148 * this._k1[i] +
      0.341796875          * this._k2[i] +
      0.041594328703703703 * this._k3[i] +
      0.400345413773148148 * this._k4[i] +
      0.061767578125       * this._k5[i] )
  }

  this.deriv( this._k6, this._w, this.t + dt * 0.875 )

  return this
}

Integrator.prototype._calculateError = function(dt) {
  //var cs1 =  0.102177372685185185 // 2825/27648
  //var cs2 =  0.000000000000000000 // 0
  //var cs3 =  0.383907903439153439 // 18575/48384
  //var cs4 =  0.244592737268518518 // 13525/55296
  //var cs5 =  0.019321986607142857 // 277/14336
  //var cs6 =  0.250000000000000000 // 1/4

  //var dc1 =  0.004293774801587301 // cs1 - c1
  //var dc2 =  0.000000000000000000 // cs2 - c2
  //var dc3 = -0.018668586093857832 // cs3 - c3
  //var dc4 =  0.034155026830808080 // cs4 - c4
  //var dc5 =  0.019321986607142857 // cs5 - c5
  //var dc6 = -0.039102202145680406 // cs6 - c6

  var error = 0
  for(var i=0; i<this.n; i++) {
    error =  this.errorReduceFunction( i, error,
      dt * (
         0.004293774801587301 * this._k1[i] +
        -0.018668586093857832 * this._k3[i] +
         0.034155026830808080 * this._k4[i] +
         0.019321986607142857 * this._k5[i] +
        -0.039102202145680406 * this._k6[i]
      ) / this._errorScale[i]
    )
  }

  return this.errorPostFunction(error)
}

Integrator.prototype._update = function(dt) {
  //var c1  =  0.097883597883597883 // 37/378
  //var c2  =  0.000000000000000000 // 0
  //var c3  =  0.402576489533011272 // 250/621
  //var c4  =  0.210437710437710437 // 125/594
  //var c5  =  0.000000000000000000 // 0
  //var c6  =  0.289102202145680406 // 512/1771

  for(var i=0; i<this.n; i++) {
    this.y[i] += dt * (
      0.097883597883597883 * this._k1[i] +
      0.402576489533011272 * this._k3[i] +
      0.210437710437710437 * this._k4[i] +
      0.289102202145680406 * this._k6[i]
    )
  }
  this.t += dt
  return this
}

Integrator.prototype._calculateErrorScale = function(dt) {
  for(var i=0; i<this.n; i++) {
    this._errorScale[i] = this.errorScaleFunction(i, dt, this.y[i], this._k1[i])
  }
  return this
}

Integrator.prototype.step = function( tLimit ) {
  // Bail out early if we're *at* the limit:
  if (Math.abs(this.t - tLimit) < this.dt * 1e-10) {
    return false;
  }

  var thisDt = this.dt;

  // Don't integrate past a tLimit, if provided:
  if( tLimit !== undefined ) {
    thisDt = thisDt > 0 ? Math.min( tLimit - this.t, thisDt ) : Math.max( tLimit - this.t, thisDt )
  }

  // Limit the magnitude of dt to dtMaxMag
  if( this.dtMaxMag !== undefined && Math.abs( thisDt ) > this.dtMaxMag ) {
    this.__log('step', 'step greater than maximum stepsize requested. dt magnitude has been limited.')
    thisDt = thisDt > 0 ? this.dtMaxMag : -this.dtMaxMag
  }

  // Limit the magnitude of dt to dtMinMag
  if( this.dtMinMag !== undefined && Math.abs( thisDt ) < this.dtMinMag ) {
    this.__log('step', 'step smaller than minimum stepsize requested. dt magnitude has been limited.')
    thisDt = thisDt > 0 ? this.dtMinMag : -this.dtMinMag
  }

  // The first derivative doesn't change even if dt does, so only calculate this once:
  this._calculateK1()

  // The scale factor per-dimension probably doesn't need to change either across a single adaptive step:
  this._calculateErrorScale(thisDt)

  var error = Infinity
  var maxError = 0
  var nextDt
  var lowerDtLimitReached = false

  while(true) {

    // Calculate intermediate k's for the proposed step:
    this._calculateKs(thisDt)

    // Calculate the max error of the proposed step:
    error = this._calculateError(thisDt)

    if( error < this.tol || lowerDtLimitReached ) {
      // Success! Exit:
      break
    }

    if( ! Number.isFinite(error) ) {
      throw new Error('ode45-cash-karp::step() NaN encountered while integrating.')
    }

    // Failure. Adapt the timestep:
    nextDt = this.safetyFactor * thisDt * Math.pow( this.tol / error, 0.2 )

    // Cut the timestep, but not by more than maxDecreaseFactor
    thisDt = maxMag( thisDt / this.maxDecreaseFactor, nextDt )

    // If stepsize too small, finish off by taking the currently proposed step and logging a warning:
    if( this.dtMinMag !== undefined && Math.abs(thisDt) < this.dtMinMag ) {
      thisDt = this.dtMinMag * (thisDt > 0 ? 1 : -1);
      this.__log('step', 'minimum stepsize reached.')
      lowerDtLimitReached = true
    }
  }

  // Apply this update:
  this._update(thisDt)

  // Calculate the next timestep size:
  nextDt = this.safetyFactor * thisDt * Math.pow( this.tol / error, 0.25 )

  // Increase the timestep for the next time around, but not by more than the maxIncreaseFactor:
  this.dt = maxMag(this.dt / this.maxDecreaseFactor, minMag( this.dt * this.maxIncreaseFactor, nextDt ));

  if( tLimit !== undefined ) {
    return Math.abs(this.t - tLimit) > this.dt * 1e-8;
  } else {
    return true
  }
}

Integrator.prototype.steps = function( n, tLimit ) {
  for(var step=0; step<n; step++) {
    if( ! this.step(tLimit) ) return false;
  }
}

function IntegratorFactory( y0, deriv, t, dt, options ) {
  return new Integrator( y0, deriv, t, dt, options )
}
