'use strict'

var fs = require('fs')
  , rkckOut = fs.openSync('rkck.dat','w')
  , rk4Out = fs.openSync('rk4.dat','w')
  , rk4 = require('ode-rk4')
  , rkck = require('../lib')

var deriv = function(dydt, y, t) {
  dydt[0] = Math.exp(-t*t * 100) * 2/Math.sqrt(Math.PI/100)
}

var t0 = -10
var tmax = 10

var i1 = rkck( [-1], deriv, t0, 1e-4, {
  tol: 1e-8,
  dtMaxMag: 0.5,
  maxIncreaseFactor: 2
})
var i2 = rk4( [-1], deriv, t0, 1e-2 )

fs.writeSync( rkckOut, i1.t + '\t' + i1.y[0] + '\t' + i1.dt + '\n' )
while( i1.step( tmax ) ) {
  fs.writeSync( rkckOut, i1.t + '\t' + i1.y[0] + '\t' + i1.dt + '\n' )
}

fs.writeSync( rk4Out, i2.t + '\t' + i2.y[0] + '\t' + i2.dt + '\n' )
for(var i=0; i<(tmax-t0)/i2.dt; i++, i2.step()) {
  fs.writeSync( rk4Out, i2.t + '\t' + i2.y[0] + '\t' + i2.dt + '\n' )
}

fs.closeSync( rkckOut )
fs.closeSync( rk4Out )
