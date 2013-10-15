context('AnnotationCommand tests')

test_that('DensCpGCommand gets its slots right', {
          foo <- DensCpGCommand('bar', 2032)
          expect_equal(foo@colName, 'bar')
          expect_equal(foo@windowSize, 2032)
})

test_that('DensCpGCommand refuses wrong data types', {
          expect_error(DensCpGCommand(2032, 2032))
          expect_error(DensCpGCommand('bar', 'foo'))
})

test_that('DensCpGCommand refuses invalid windowSize', {
          expect_error(DensCpGCommand('foo', -1))
          expect_error(DensCpGCommand('foo', 0))
          expect_error(DensCpGCommand('foo', 1))
})

          


