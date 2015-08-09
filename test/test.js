'use strict'

var rk45 = require('../lib')
  , assert = require('chai').assert
  , richardson = require('richardson-extrapolation')

var ctors = {
  //'float32': Float32Array,
  'float64': Float64Array,
  'array': function(){ return arguments[0] }
}


Object.keys(ctors).forEach(function(dtype) {
  var ctor = ctors[dtype]

  describe('rk45 integration (' + dtype + ')', function() {

    describe('setup', function() {
      var integrator, f, y0, t0, n

      beforeEach(function() {
        f = function(dydt, y) { dydt[0] = -y[0] }
        t0 = 1.5
        y0 = new ctor([1])
        n = 10

        integrator = rk45( new ctor([1]), function(){}, 1, 1)
      })

      it('creates work arrays of the same type as the input',function() {
        assert.equal( integrator._w.constructor,  y0.constructor )
        assert.equal( integrator._k1.constructor, y0.constructor )
        assert.equal( integrator._k2.constructor, y0.constructor )
        assert.equal( integrator._k3.constructor, y0.constructor )
        assert.equal( integrator._k4.constructor, y0.constructor )
        assert.equal( integrator._k5.constructor, y0.constructor )
        assert.equal( integrator._k6.constructor, y0.constructor )
      })

      it('creates work arrays of the same size as the input',function() {
        assert.equal( integrator._w.length,  y0.length )
        assert.equal( integrator._k1.length, y0.length )
        assert.equal( integrator._k2.length, y0.length )
        assert.equal( integrator._k3.length, y0.length )
        assert.equal( integrator._k4.length, y0.length )
        assert.equal( integrator._k5.length, y0.length )
        assert.equal( integrator._k6.length, y0.length )
      })
    })

    describe('adaptive timestepping', function() {

      it('is sign-independent in the independent variable', function() {
        var f = function(dydt, y) {
          dydt[0] = -y[1]
          dydt[1] =  y[0]
        }
        // Integrate around a circle and confirm that it doesn't matter
        // whether we integrate one way or the other:
        var i1 = rk45( new ctor([1,0]), f, 0, 1e10).step()
        var i2 = rk45( new ctor([1,0]), f, 0, -1e10).step()

        assert.closeTo( i1.y[0], i2.y[0], 1e-8, 'x-coordinates are equal' )
        assert.closeTo( i1.y[1], -i2.y[1], 1e-8, 'y-coordinates are opposite' )
        assert.closeTo( i1.dt, -i2.dt, 1e-8, 'dt has been adapted identically' )
      })

      it('is scale-independent in calculating the error', function() {
        var f = function(dydt, y) {
          dydt[0] = -y[1]
          dydt[1] =  y[0]
        }
        // Integration around a circle is scale-independent in dt, so
        // ensure that the adaptation is the same no matter the radius:
        var i1 = rk45( new ctor([1e30,0]), f, 0, 1e10).step()
        var i2 = rk45( new ctor([1e-30,0]), f, 0, 1e10).step()
        assert.closeTo( i1.dt, i2.dt, 1e-3 )
      })

      it('updates dt according to the timestep taken', function() {
        // Integrate dy/dt = constant
        var c = 5.2
        var f = function(dydt, y) { dydt[0] = c }
        var i = rk45( new ctor([0]), f, 0, 1).step()
        assert.closeTo( i.y[0], i.t * c, 1e-3, 'answer is correct' )
      })

      it('integrates the 0 without incident', function() {
        // Just to make sure there aren't any divide-by-zero issues
        var f = function(dydt, y) { dydt[0] = 0 }
        var i = rk45( new ctor([0]), f, 0, 1).step()
        assert.closeTo( i.y[0], 0, 1e-3, 'answer is correct' )
        assert.closeTo( i.dt, 10, 1e-3, 'dt has been increased for the next step' )
      })

      it('increases the timestep by no more than maxIncreaseFactor if tolerance met', function() {
        // Integrating a straight line should increase the timestep by maxIncreaseFactor:
        var dt0 = 15
        var factor = 11
        var f = function(dydt, y) { dydt[0] = 1 }
        var i = rk45( new ctor([0]), f, 0, dt0, {maxIncreaseFactor: factor}).step()
        assert.closeTo( i.dt, dt0 * factor, 1e-3, 'increased dt by maxIncreaseFactor' )
      })

      it('throws an error on step size hitting dtMin', function() {
        var i
        var f = function(dydt, y) { dydt[0] = Math.cos(1e10*y[0]) }
        assert.throws(function() {
          i = rk45( new ctor([200]), f, 0, 1, {dtMin: 1e-4})
          i.step()
        },Error,/minimum stepsize exceeded/)

        assert( Number.isFinite(i.y[0]), 'y is finite' )
        assert( i.y[0] !== 200, 'y has been timestepped' )
      })

    })

    describe('convergence', function() {

      it('total accumulated error of high order scheme is order O(h^5)', function() {

        var result = richardson(function(h) {
          // Integrate around a circle at an accelerating rate
          var f = function(dydt, y, t) {
            var s =  Math.sin(t * Math.PI) * Math.PI / 2
            dydt[0] = -y[1]* 2 * Math.PI * s
            dydt[1] =  y[0]* 2 * Math.PI * s
          }
          var i = rk45( new ctor([1,0,0]), f, 0, h )

          var n = Math.floor(1/h+0.5)
          for(var j=0; j<n; j++) {
            i._calculateK1()._step()._update()
          }

          // Return the distance from the expected endpoint:
          return Math.sqrt( Math.pow(i.y[0]-1,2) + Math.pow(i.y[1],2) )

          // Step size chosen to be the minimum before we lose too much precision in the
          // float32 scheme to calculate anything meaningful:
        }, 1/18, { f: 0 } )

        assert.closeTo( result.n, 5, 0.15, 'n ~ 5' )
      })

    })

  })

})
