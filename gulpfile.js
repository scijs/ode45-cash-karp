var gulp = require('gulp')
  , mdEqs = require('gulp-markdown-equations')
  , tap = require('gulp-tap')
  , filter = require('gulp-filter')
  , latex = require('gulp-latex')
  , pdftocairo = require('gulp-pdftocairo')


gulp.task('mdtex',function() {
  var texFilter = filter('*.tex')
  var mdFilter = filter('*.md')

  var transform = mdEqs({
    defaults: {
      display: { margin: '1pt 1pt 1pt -10pt' },
      inline:  { margin: '1pt 1pt 1pt -5pt' }
    }
  })

  return gulp.src('*.mdtex')
    .pipe(transform)
    .pipe(texFilter)
    .pipe(tap(function(file) { console.log(file.contents.toString()) }))
    .pipe(latex())
    .pipe(pdftocairo({format: 'png'}))
    .pipe(gulp.dest('docs/images'))
    .pipe(tap(function(file) {
      transform.completeSync(file,function() {
        var img = '<img alt="'+this.alt+'" valign="middle" src="'+this.path+
                  '" width="'+this.width/2+'" height="'+this.height/2+'">'
        return this.inline ? img : '<p align="center">'+img+'</p>'
      })
    }))
    .pipe(texFilter.restore()).pipe(mdFilter)
    .pipe(gulp.dest('./'))
})
