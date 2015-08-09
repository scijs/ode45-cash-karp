'use strict'

module.exports = IntegratorFactory

var Integrator = function Integrator( y0, deriv, t, dt, options ) {

  var opts = options || {}
  this.tol = opts.tol===undefined ? 1e-8 : opts.tol
  this.safety = opts.safety===undefined ? 0.9 : opts.safety
  this.maxIncreaseFactor = opts.maxIncreaseFactor===undefined ? 10 : opts.maxIncreaseFactor
  this.maxDecreaseFactor = opts.maxDecreaseFactor===undefined ? 10 : opts.maxDecreaseFactor
  this.dtMin = opts.dtMin===undefined ? 0 : Math.abs(opts.dtMin)
  //this.dtClip = opts.dtClip===undefined ? undefined : Math.abs(opts.dtClip)

  // Bind variables to this:
  this.deriv = deriv
  this.y = y0
  this.n = this.y.length
  this.dt = dt
  this.t = t

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

Integrator.prototype._step = function() {
  var i

  var a21 =  0.200000000000000000 // 1/5
  var a31 =  0.075000000000000000 // 3/40
  var a32 =  0.225000000000000000 // 9/40
  var a41 =  0.300000000000000000 // 3/10
  var a42 = -0.900000000000000000 // -9/10
  var a43 =  1.200000000000000000 // 6/5
  var a51 = -0.203703703703703703 // -11/54
  var a52 =  2.500000000000000000 // 5/2
  var a53 = -2.592592592592592592 // -70/27
  var a54 =  1.296296296296296296 // 35/27
  var a61 =  0.029495804398148148 // 1631/55296
  var a62 =  0.341796875000000000 // 175/512
  var a63 =  0.041594328703703703 // 575/13824
  var a64 =  0.400345413773148148 // 44275/110592
  var a65 =  0.061767578125000000 // 253/4096

//var b1  =  0.000000000000000000 // 0
  var b2  =  0.200000000000000000 // 1/5
  var b3  =  0.300000000000000000 // 3/10
  var b4  =  0.600000000000000000 // 3/5
//var b5  =  1.000000000000000000 // 1
  var b6  =  0.875000000000000000 // 7/8

  // Same for every step, so don't repeat:
  //this.deriv( this._k1, this.y, this.t )

  for(i=0; i<this.n; i++) {
    this._w[i] = this.y[i] + this.dt * (
      a21 * this._k1[i] )
  }

  this.deriv( this._k2, this._w, this.t + this.dt * b2 )

  for(i=0; i<this.n; i++) {
    this._w[i] = this.y[i] + this.dt * (
      a31 * this._k1[i] +
      a32 * this._k2[i] )
  }

  this.deriv( this._k3, this._w, this.t + this.dt * b3 )

  for(i=0; i<this.n; i++) {
    this._w[i] = this.y[i] + this.dt * (
      a41 * this._k1[i] +
      a42 * this._k2[i] +
      a43 * this._k3[i] )
  }

  this.deriv( this._k4, this._w, this.t + this.dt * b4 )

  for(i=0; i<this.n; i++) {
    this._w[i] = this.y[i] + this.dt * (
      a51 * this._k1[i] +
      a52 * this._k2[i] +
      a53 * this._k3[i] +
      a54 * this._k4[i] )
  }

  this.deriv( this._k5, this._w, this.t + this.dt /* * b5 */ )

  for(i=0; i<this.n; i++) {
    this._w[i] = this.y[i] + this.dt * (
      a61 * this._k1[i] +
      a62 * this._k2[i] +
      a63 * this._k3[i] +
      a64 * this._k4[i] +
      a65 * this._k5[i] )
  }

  this.deriv( this._k6, this._w, this.t + this.dt * b6 )

  return this
}

Integrator.prototype._error = function() {
//var cs1 =  0.102177372685185185 // 2825/27648
//var cs2 =  0.000000000000000000 // 0
//var cs3 =  0.383907903439153439 // 18575/48384
//var cs4 =  0.244592737268518518 // 13525/55296
//var cs5 =  0.019321986607142857 // 277/14336
//var cs6 =  0.250000000000000000 // 1/4

  var dc1 =  0.004293774801587301 // cs1 - c1
//var dc2 =  0.000000000000000000 // cs2 - c2
  var dc3 = -0.018668586093857832 // cs3 - c3
  var dc4 =  0.034155026830808080 // cs4 - c4
  var dc5 =  0.019321986607142857 // cs5 - c5
  var dc6 = -0.039102202145680406 // cs6 - c6

  var maxError = 0
  for(var i=0; i<this.n; i++) {
    maxError =  Math.max( maxError,
      Math.abs(
        this.dt * (
          dc1*this._k1[i] +
          //dc2*this._k2[i] +
          dc3*this._k3[i] +
          dc4*this._k4[i] +
          dc5*this._k5[i] +
          dc6*this._k6[i]
        ) / this._errorScale[i]
      )
    )
  }

  return maxError
}

Integrator.prototype._update = function() {
  var c1  =  0.097883597883597883 // 37/378
//var c2  =  0.000000000000000000 // 0
  var c3  =  0.402576489533011272 // 250/621
  var c4  =  0.210437710437710437 // 125/594
//var c5  =  0.000000000000000000 // 0
  var c6  =  0.289102202145680406 // 512/1771

  for(var i=0; i<this.n; i++) {
    this.y[i] += this.dt * (
      c1*this._k1[i] +
    //c2*this._k2[i] +
      c3*this._k3[i] +
      c4*this._k4[i] +
    //c5*this._k5[i] +
      c6*this._k6[i]
    )
  }
  this.t += this.dt
  return this
}

Integrator.prototype._calculateErrorScale = function() {
  for(var i=0; i<this.n; i++) {
    this._errorScale[i] = Math.max( Math.abs(this.y[i]), Math.abs(this.dt * this._k1[i]) ) + 1e-50;
  }
  return this
}

Integrator.prototype.step = function() {
  //console.log('\nTIMESTEP, dt =',this.dt)

  this._calculateK1()
  this._calculateErrorScale()

  var error = Infinity
  var maxError = 0
  var prevDt = this.dt
  var nextDt

  while(true) {

    //console.log('ATTEMPT, dt =',this.dt)

    // Calculate intermediate k's for the proposed step:
    this._step()

    // Calculate the max error of the proposed step:
    error = this._error()
    //console.log('  error =',error)

    if( error < this.tol ) {
      // Success! Exit:
      break;
    }

    // Failure. Adapt the timestep:
    nextDt = this.safety * this.dt * Math.pow( this.tol / error, 0.2 );

    // Cut the timestep, but not by more than maxDecreaseFactor
    this.dt = (this.dt > 0 ? Math.max : Math.min)( this.dt/this.maxDecreaseFactor, nextDt )

    // If stepsize too small, finish off by taking the currently proposed step and throwing an error
    if( Math.abs(this.dt) < this.dtMin ) {
      this._update()
      throw new Error('rk45(): minimum stepsize exceeded.')
    }
  }

  // Apply this update;
  this._update()

  // Calculate a new timestep:
  nextDt = this.safety * this.dt * Math.pow( this.tol / error, 0.25 );

  this.dt = (this.dt > 0 ? Math.min : Math.max)( this.dt * this.maxIncreaseFactor, nextDt )

  return this
}

Integrator.prototype.steps = function( n ) {
  for(var step=0; step<n; step++) {
    this.step()
  }
  return this
}

function IntegratorFactory( y0, deriv, t, dt, options ) {
  return new Integrator( y0, deriv, t, dt, options )
}
