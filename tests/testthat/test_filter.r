context('FilterCommand tests')

library(minfi)

#
# Mock objects for testing
#
mockDetectionP <- matrix(c(1, 0, 0.1,   0, 1,
                           0, 0,   0,   0, 1,
                           0, 1,   0,   0, 0,
                           0, 0,   0,   0, 1,
                           0, 1,   0, 0.1, 0
                           ), nrow=5, ncol=5, byrow=TRUE)

zeroDetectionP <- matrix(c(0, 0, 0, 0, 0,
                           0, 0, 0, 0, 0,
                           0, 0, 0, 0, 0,
                           0, 0, 0, 0, 0,
                           0, 0, 0, 0, 0
                           ), nrow=5, ncol=5, byrow=TRUE)

emptyMethylSet <- MethylSet()

mockMethylSet <- MethylSet(Meth=matrix(1:5, nrow=5, ncol=5), 
                           Unmeth=matrix(1:5, nrow=5, ncol=5, byrow=TRUE))
featureNames(mockMethylSet) <- paste0('F', 1:5)
sampleNames(mockMethylSet) <- paste0('S', 1:5)

#
# FilterCommand tests
#
test_that('FilterCommand derived classes refuse to run on empty objects',
          {
            foo <- kOverADetPFilterCommand(mockDetectionP, byRow=TRUE, k=2, a=0.07)
            expect_error(execute(foo, emptyMethylSet))
          })
          
#
# AtomicFilterCommand tests
#
test_that('AtomicFilterCommand accessors work well on derived classes',
          {
            foo <- kOverADetPFilterCommand(mockDetectionP, byRow=TRUE, k=2, a=0.07)
            expect_equal(getByRow(foo), TRUE)
            setByRow(foo) <- FALSE
            expect_equal(getByRow(foo), FALSE)
            expect_error(setByRow(foo) <- 42)
            expect_error(setByRow(foo) <- 'foobar')
            expect_error(setByRow(foo) <- c(TRUE, FALSE, TRUE))
          })

#
# DetPFilterCommand tests
#
test_that('DetPFilterCommand accessors work well on derived classes',
          {
            foo <- kOverADetPFilterCommand(mockDetectionP, byRow=TRUE, k=2, a=0.07)
            expect_equal(getDetectionP(foo), mockDetectionP)
            setDetectionP(foo) <- zeroDetectionP
            expect_equal(getDetectionP(foo), zeroDetectionP)
            expect_error(setDetectionP(foo) <- emptyMethylSet)
            expect_error(setDetectionP(foo) <- 42)
            expect_error(setDetectionP(foo) <- 'foobar')
            expect_error(setDetectionP(foo) <- TRUE)
          })

#
# KOverADetPFilterCommand tests
#
test_that('KOverADetPFilterCommand specific accessors work well',
          {
            foo <- kOverADetPFilterCommand(mockDetectionP, byRow=TRUE, k=2, a=0.07)
            expect_equal(getK(foo), 2)
            expect_equal(getA(foo), 0.07)
            setK(foo) <- 3
            expect_equal(getK(foo), 3)
            setA(foo) <- 0.00007
            expect_equal(getA(foo), 0.00007)
            expect_error(setK(foo) <- -1)
            expect_error(setK(foo) <- 9999)
            expect_error(setK(foo) <- c(4, 8, 15, 16, 23, 42))
            expect_error(setK(foo) <- NA)
            expect_error(setA(foo) <- -1)
            expect_error(setA(foo) <- 9999)
            expect_error(setA(foo) <- c(4, 8, 15, 16, 23, 42))
            expect_error(setA(foo) <- NA)
          })

test_that('KOverADetPfilterCommand percentage constructor works as expected',
          {
            foo <- kOverADetPFilterCommandFromFraction(mockDetectionP, byRow=TRUE, fraction=0.4, 
                                                       a=0.07)
            expect_equal(getK(foo), 2)
            expect_equal(getA(foo), 0.07)
            expect_error(foo <- kOverADetPFilterCommandFromFraction(mockDetectionP, byRow=TRUE, 
                                                                    fraction=-0.01, a=0.07))
            expect_error(foo <- kOverADetPFilterCommandFromFraction(mockDetectionP, byRow=TRUE, 
                                                                    fraction=1.01, a=0.07))
            expect_error(foo <- kOverADetPFilterCommandFromFraction(mockDetectionP, byRow=TRUE, 
                                                                    fraction=c(0.5, 0.7), a=0.07))

          })

cmd1 <- kOverADetPFilterCommand(mockDetectionP, byRow=TRUE, k=2, a=0.5)
cmd2 <- kOverADetPFilterCommand(mockDetectionP, byRow=TRUE, k=2, a=0.01)
cmd3 <- kOverADetPFilterCommand(mockDetectionP, byRow=FALSE, k=1, a=0.5)
cmd4 <- kOverADetPFilterCommand(mockDetectionP, byRow=FALSE, k=1, a=0.01)

test_that('KOverADetPFilterCommand execution breaks on wrong data types',
          {
            expect_error(execute(mockMethylSet, cmd1))
            expect_error(execute(mockMethylSet, cmd2))
            expect_error(execute(mockMethylSet, cmd3))
            expect_error(execute(mockMethylSet, cmd4))
          })

test_that('KOverADetPFilterCommand fails on empty example',
          {
            expect_error(execute(cmd1, emptyMethylSet))
            expect_error(execute(cmd2, emptyMethylSet))
            expect_error(execute(cmd3, emptyMethylSet))
            expect_error(execute(cmd4, emptyMethylSet))
          })

test_that('KOverADetPFilterCommand execution works correctly on examples',
          {
            bar1 <- execute(cmd1, mockMethylSet)
            expect_equal(featureNames(bar1), c('F2', 'F3', 'F4', 'F5'))
            expect_equal(sampleNames(bar1), c('S1', 'S2', 'S3', 'S4', 'S5'))
            expect_equal(nrow(bar1), c(Features=4))
            expect_equal(ncol(bar1), c(Samples=5))
            bar2 <- execute(cmd2, mockMethylSet)
            expect_equal(featureNames(bar2), c('F2', 'F3', 'F4'))
            expect_equal(sampleNames(bar2), c('S1', 'S2', 'S3', 'S4', 'S5'))
            expect_equal(nrow(bar2), c(Features=3))
            expect_equal(ncol(bar2), c(Samples=5))
            bar3 <- execute(cmd3, mockMethylSet)
            expect_equal(featureNames(bar3), c('F1', 'F2', 'F3', 'F4', 'F5'))
            expect_equal(sampleNames(bar3), c('S3', 'S4'))
            expect_equal(nrow(bar3), c(Features=5))
            expect_equal(ncol(bar3), c(Samples=2))
            bar4 <- execute(cmd4, mockMethylSet)
            expect_equal(featureNames(bar4), c('F1', 'F2', 'F3', 'F4', 'F5'))
            expect_equal(sampleNames(bar4), character(0))
            expect_equal(nrow(bar4), c(Features=5))
            expect_equal(ncol(bar4), c(Samples=0))
          })

