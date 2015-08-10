'use strict'

var rkck = require('../lib')

var evaluations = 0
var deriv = function(dydt, y, t) {
  evaluations ++
  dydt[0] = 1/t * Math.cos(1/t)
}

var ta = 0.01
var tb = 1

var i = rkck( [-1], deriv, ta, 1e-8, {
  tol: 5e-8,
  maxIncreaseFactor: 2
})

var rkck_y0 = i.y[0]
i.steps( Infinity, tb )
var rkck_y1 = i.y[0]

var actual = -0.34255274804359265
console.log('Absolute error:',Math.abs( actual - (rkck_y1-rkck_y0) ))
console.log('Derivative evaluations:', evaluations )
