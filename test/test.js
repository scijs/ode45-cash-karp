'use strict'

var ode45 = require('../lib')
  , assert = require('chai').assert
  , richardson = require('richardson-extrapolation')

var ctors = {
  'float32': Float32Array,
  'float64': Float64Array,
  'array': function(){ return arguments[0] }
}


Object.keys(ctors).forEach(function(dtype) {
  var ctor = ctors[dtype]

  describe('ode45 integration (' + dtype + ')', function() {

    describe('setup', function() {
      var integrator, f, y0, t0, n

      beforeEach(function() {
        f = function(dydt, y) { dydt[0] = -y[0] }
        t0 = 1.5
        y0 = new ctor([1])
        n = 10

        integrator = ode45( new ctor([1]), function(){}, 1, 1)
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
        var i1 = ode45( new ctor([1,0]), f, 0, 1e4)
        i1.step()
        var i2 = ode45( new ctor([1,0]), f, 0, -1e4)
        i2.step()

        assert.closeTo( i1.y[0], i2.y[0], 1e-6, 'x-coordinates are equal' )
        assert.closeTo( i1.y[1], -i2.y[1], 1e-6, 'y-coordinates are opposite' )
        assert.closeTo( i1.dt, -i2.dt, 1e-6, 'dt has been adapted identically' )
      })

      it('is scale-independent in calculating the error', function() {
        var f = function(dydt, y) {
          dydt[0] = -y[1]
          dydt[1] =  y[0]
        }
        // Integration around a circle is scale-independent in dt, so
        // ensure that the adaptation is the same no matter the radius:
        var i1 = ode45( new ctor([1e5,0]), f, 0, 1e4)
        i1.step()
        var i2 = ode45( new ctor([1e-5,0]), f, 0, 1e4)
        i2.step()
        assert.closeTo( i1.dt, i2.dt, 1e-2 )
      })

      it('throws an error if NaN encountered', function() {
        var f = function(dydt, y) { dydt[0] = Math.pow(y[0],4) }
        assert.throws(function() {
          var i = ode45( new ctor([100,0]), f, 0, 1)
          i.steps(10)
        },Error,/NaN encountered/)
      })

      it('updates dt according to the timestep taken', function() {
        // Integrate dy/dt = constant
        var c = 5.2
        var f = function(dydt, y) { dydt[0] = c }
        var i = ode45( new ctor([0]), f, 0, 1)
        i.step()
        assert.closeTo( i.y[0], i.t * c, 1e-3, 'answer is correct' )
      })

      it('integrates with dt>0 and a limit', function() {
        var c = 5.2
        var t0 = 1
        var t1 = 3
        var dt = 1
        var f = function(dydt, y) { dydt[0] = c }
        var i = ode45( new ctor([0]), f, t0, dt)
        i.step( t1 )
        assert.closeTo( i.y[0], (2-t0) * c, 1e-3, 'answer is correct' )
        assert.closeTo( i.t, 2, 1e-3, 'dt has been clipped')
        i.step( t1 )
        assert.closeTo( i.y[0], (t1-t0) * c, 1e-3, 'answer is correct' )
        assert.closeTo( i.t, t1, 1e-3, 'dt has been clipped')
      })

      it('integrates with dt<0 and a limit', function() {
        var c = 5.2
        var t0 = 3
        var t1 = 1
        var dt = -1
        var f = function(dydt, y) { dydt[0] = c }
        var i = ode45( new ctor([0]), f, t0, dt)

        assert.isTrue( i.step( t1 ) )
        assert.closeTo( i.y[0], (2-t0) * c, 1e-3, 'answer is correct' )
        assert.closeTo( i.t, t0+dt, 1e-3, 'dt is updated')
        assert.isTrue( i.step( t1 ) )
        assert.closeTo( i.y[0], (t1-t0) * c, 1e-3, 'answer is correct' )
        assert.closeTo( i.t, t1, 1e-3, 'dt has been clipped')
      })

      it('integrates the 0 without incident', function() {
        // Just to make sure there aren't any divide-by-zero issues
        var f = function(dydt, y) { dydt[0] = 0 }
        var i = ode45( new ctor([0]), f, 0, 1)
        i.step()
        assert.closeTo( i.y[0], 0, 1e-3, 'answer is correct' )
        assert.closeTo( i.dt, 10, 1e-3, 'dt has been increased for the next step' )
      })

      it('increases the timestep by no more than maxIncreaseFactor if tolerance met', function() {
        // Integrating a straight line should increase the timestep by maxIncreaseFactor:
        var dt0 = 15
        var factor = 11
        var f = function(dydt, y) { dydt[0] = 1 }
        var i = ode45( new ctor([0]), f, 0, dt0, {maxIncreaseFactor: factor})
        i.step()
        assert.closeTo( i.dt, dt0 * factor, 1e-3, 'increased dt by maxIncreaseFactor' )
      })
    })

    describe('stepsize limiting', function() {
      it('doesn\'t decrease the step size past dtMinMag (dt > 0)', function() {
        var dtMinMag = 1e-4
        var i
        var f = function(dydt, y) { dydt[0] = Math.cos(1e5*y[0]) }
        i = ode45( new ctor([0]), f, 0, 1, {dtMinMag: dtMinMag})
        i.step()
        assert( Number.isFinite(i.y[0]), 'y is finite' )
        assert( i.y[0] !== 0, 'y has been timestepped' )
        assert.closeTo( i.t, dtMinMag, 1e-8, 'dt is at the lower limit' )
        i.step()
        assert.closeTo( i.t, 2 * dtMinMag, 1e-8, 'dt is at the lower limit' )
      })

      it('doesn\'t decrease the step size past dtMinMag (dt < 0)', function() {
        var dtMinMag = 1e-4
        var i
        var f = function(dydt, y) { dydt[0] = Math.cos(1e5*y[0]) }
        i = ode45( new ctor([0]), f, 0, -1, {dtMinMag: dtMinMag})
        i.step()
        assert( Number.isFinite(i.y[0]), 'y is finite' )
        assert( i.y[0] !== 0, 'y has been timestepped' )
        assert.closeTo( i.t, -dtMinMag, 1e-8, 'dt is at the lower limit' )
        i.step()
        assert.closeTo( i.t, -2 * dtMinMag, 1e-8, 'dt is at the lower limit' )
      })

      it('doesn\'t increase the step size past dtMax (dt > 0)', function() {
        var dtMaxMag = 2
        var i, dt0 = 1e4
        var f = function(dydt, y) { dydt[0] = 1 }
        i = ode45( new ctor([0]), f, 0, dt0, {dtMaxMag: dtMaxMag})
        i.step()
        assert( Number.isFinite(i.y[0]), 'y is finite' )
        assert( i.y[0] !== 0, 'y has been timestepped' )
        assert.closeTo( i.t, dtMaxMag, 1e-8, 'dt is at the upper limit' )
        i.step()
        assert.closeTo( i.t, dtMaxMag * 2, 1e-8, 'dt is at the upper limit' )
      })

      it('doesn\'t increase the step size past dtMax (dt < 0)', function() {
        var dtMaxMag = 2
        var i, dt0 = -1e4
        var f = function(dydt, y) { dydt[0] = 1 }
        i = ode45( new ctor([0]), f, 0, dt0, {dtMaxMag: dtMaxMag})
        i.step()
        assert( Number.isFinite(i.y[0]), 'y is finite' )
        assert( i.y[0] !== 0, 'y has been timestepped' )
        assert.closeTo( i.t, -dtMaxMag, 1e-8, 'dt is at the upper limit' )
        i.step()
        assert.closeTo( i.t, -dtMaxMag * 2, 1e-8, 'dt is at the upper limit' )
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
          var i = ode45( new ctor([1,0,0]), f, 0, h )

          var n = Math.floor(1/h+0.5)
          for(var j=0; j<n; j++) {
            i._calculateK1()
            i._calculateKs()
            i._update()
          }

          // Return the distance from the expected endpoint:
          return Math.sqrt( Math.pow(i.y[0]-1,2) + Math.pow(i.y[1],2) )

          // Step size chosen to be the minimum before we lose too much precision in the
          // float32 scheme to calculate anything meaningful. There's probably a better
          // equation for this, but long story short, fifth order convergence is difficult
          // to test in 32 bit floating point.
        }, 1/18, { f: 0 } )

        assert.closeTo( result.n, 5, 0.15, 'n ~ 5' )
      })

    })

    describe

  })

})
