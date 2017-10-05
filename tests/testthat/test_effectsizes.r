context('Effectsizes tests')

library(epilab)

#
# cliffDelta tests
#

test_that('cliffDelta works only with numerical samples',
          {
            expect_error(cliffDelta(letters[1:10], 1:10))
            expect_error(cliffDelta(factor(c('A', 'B', 'C')), 1:10))
          })

test_that('cliffDelta works correctly on example',
          {
            expect_equal(cliffDelta(1:10, 2:11), -0.19)
            expect_equal(cliffDelta(2:11, 1:10), 0.19)
            expect_equal(cliffDelta(2:11, 2:11), 0)
          })

#
# camblorAC tests
#

test_that('camblorAC works only with numerical samples',
          {
            expect_error(camblorAC(letters[1:10], 1:10))
            expect_error(camblorAC(factor(c('A', 'B', 'C')), 1:10))
          })

test_that('camblorAC fails with wrong n values',
          {
            expect_error(camblorAC(1:10, 1:10, n=4)) 
            expect_error(camblorAC(1:10, 1:10, n=3)) 
            expect_error(camblorAC(1:10, 1:10, n=2)) 
            expect_error(camblorAC(1:10, 1:10, n=1)) 
            expect_error(camblorAC(1:10, 1:10, n=0)) 
            expect_error(camblorAC(1:10, 1:10, n=-1)) 
            expect_error(camblorAC(1:10, 1:10, n=-2)) 
          })

test_that('camblorAC works correctly on example',
          {
            expect_equal(camblorAC(1:10, 1:10), 1.0)
            expect_equal(camblorAC(-10:-1, 1:10), 0.05304149)
            expect_equal(camblorAC(1:10, 4:14), 0.5710375, tolerance=1e-7)
            expect_equal(camblorAC(1:10, 1e7:(1e7 + 100)), 0.0, tolerance=1e-7)
          })

#
# oddsRatio tests
#

mockMatrix <- matrix(c(1:4), ncol=2)
mockBigMatrix <- matrix(c(1:6), ncol=2)
mockWideMatrix <- matrix(c(1:6), ncol=3)

test_that('oddsRatio works only with 2x2 matrices',
          {
            expect_error(oddsRatio('foobar', TRUE))
            expect_error(oddsRatio(1:10, FALSE))
            expect_error(oddsRatio(factor(1:10), TRUE))
            expect_error(oddsRatio(list(1, 2, 3), FALSE))
            expect_error(oddsRatio(mockBigMatrix, TRUE))
          })

test_that('oddsRatio breaks on wrong changeRows values',
          {
            expect_error(oddsRatio(mockMatrix, 'foobar'))
            expect_error(oddsRatio(mockMatrix, 1:10))
            expect_error(oddsRatio(mockMatrix, c(TRUE, FALSE, TRUE)))
            expect_error(oddsRatio(mockMatrix, factor(1:3)))
          })

test_that('oddsRatio works correctly on example',
          {
            expect_equal(oddsRatio(mockMatrix, FALSE), 2 / 3)
            expect_equal(oddsRatio(mockMatrix, TRUE), 1.5)
          })

#
# oddsRatioLevel tests
#

rownames(mockMatrix) <- c('foo', 'bar')
rownames(mockBigMatrix) <- c('foo', 'bar', 'foobar')

test_that('oddsRatioLevel works only with nx2 matrices',
          {
            expect_error(oddsRatioLevel('foobar', 'foo'))
            expect_error(oddsRatioLevel(1:10, 'foo'))
            expect_error(oddsRatioLevel(factor(1:10), 'foo'))
            expect_error(oddsRatioLevel(list(1, 2, 3), 'foo'))
            expect_error(oddsRatioLevel(mockWideMatrix, 'foo'))
            expect_error(oddsRatioLevel(mockMatrix, 'foo'))
          })

test_that('oddsRatioLevel breaks on wrong aLevel values',
          {
            expect_error(oddsRatioLevel(mockMatrix, 'foobar'))
            expect_error(oddsRatioLevel(mockMatrix, 1:10))
            expect_error(oddsRatioLevel(mockMatrix, c(TRUE, FALSE, TRUE)))
            expect_error(oddsRatioLevel(mockMatrix, factor(1:3)))
          })

test_that('oddsRatioLevel works correctly on example',
          {
            expect_equal(oddsRatioLevel(mockBigMatrix, 'foo'), 11 / 20)
            expect_equal(oddsRatioLevel(mockBigMatrix, 'bar'), 1)
            expect_equal(oddsRatioLevel(mockBigMatrix, 'foobar'), 27 / 18)
          })

