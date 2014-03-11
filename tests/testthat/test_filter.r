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
rownames(mockDetectionP) <- paste0('F', 1:5)
colnames(mockDetectionP) <- paste0('S', 1:5)

zeroDetectionP <- matrix(c(0, 0, 0, 0, 0,
                           0, 0, 0, 0, 0,
                           0, 0, 0, 0, 0,
                           0, 0, 0, 0, 0,
                           0, 0, 0, 0, 0
                           ), nrow=5, ncol=5, byrow=TRUE)
rownames(zeroDetectionP) <- paste0('F', 1:5)
colnames(zeroDetectionP) <- paste0('S', 1:5)

emptyMethylSet <- MethylSet()

mockMethylSet <- MethylSet(Meth=matrix(1:5, nrow=5, ncol=5), 
                           Unmeth=matrix(1:5, nrow=5, ncol=5, byrow=TRUE))
featureNames(mockMethylSet) <- paste0('F', 1:5)
sampleNames(mockMethylSet) <- paste0('S', 1:5)

mockMethylSetPlus <- MethylSet(Meth=matrix(1:6, nrow=5, ncol=6), 
                               Unmeth=matrix(1:6, nrow=5, ncol=6, byrow=TRUE))
featureNames(mockMethylSetPlus) <- paste0('F', 1:5)
sampleNames(mockMethylSetPlus) <- paste0('S', 1:6)

#
# FilterCommandIndices tests
#
mockRows <- c(TRUE, FALSE, TRUE, FALSE)
mockCols <- c(FALSE, TRUE, TRUE, FALSE)

test_that('FilterCommandIndices specific accessors work well',
          {
            foo <- filterCommandIndices(rows=mockRows, cols=mockCols)
            expect_equal(getRows(foo), mockRows)
            expect_equal(getCols(foo), mockCols)
            setRows(foo) <- mockCols
            setCols(foo) <- mockRows
            expect_equal(getRows(foo), mockCols)
            expect_equal(getCols(foo), mockRows)
          })

test_that('FilterCommandIndices refuses wrong data types',
          {
            expect_error(foo <- filterCommandIndices(rows=1:3, cols=mockCols))
            expect_error(foo <- filterCommandIndices(rows=mockRows, cols=1:7))
          })

#
# FilterCommand tests
#
test_that('FilterCommand derived classes refuse to run on empty objects',
          {
            foo <- kOverADetPFilterCommand(mockDetectionP, byRow=TRUE, k=2, a=0.07)
            expect_error(execute(foo, emptyMethylSet))
          })
          
test_that('FilterCommand fails on objects without dimensions',
          {
            foo <- kOverADetPFilterCommand(mockDetectionP, byRow=TRUE, k=2, a=0.07)
            expect_error(execute(foo, 'foobar'))
            expect_error(execute(foo, list()))
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

test_that('KOverADetPFilterCommand fails on different dimension names',
          {
            expect_error(execute(cmd1, mockMethylSetPlus),
                         regexp='must be included in dimension names')
          })

#
# FilterCommandList tests
#
test_that('FilterCommandList gets its slots right', 
          {
            foo <- filterCommandList(cmd1, cmd2, cmd3, cmd4)
            cmdList <- list(cmd1, cmd2, cmd3, cmd4)
            expect_equal(getCommandList(foo), cmdList)
          })

test_that('FilterCommandList refuses wrong data types', 
          {
            expect_error(filterCommandList(2032))
          })

cmdl1 <- filterCommandList(cmd1)
cmdl2 <- filterCommandList(cmd1, cmd2)
cmdl3 <- filterCommandList(cmd1, cmd2, cmd3)
cmdl4 <- filterCommandList(cmd1, cmd2, cmd3, cmd4)

test_that('FilterCommandList execution breaks on wrong data types', 
          {
            expect_error(execute(mockMethylSet, cmdl1))
          })

test_that('AnnotationCommandList execution works correctly on example', 
          {
            bar1 <- execute(cmdl1, mockMethylSet)
            expect_equal(featureNames(bar1), c('F2', 'F3', 'F4', 'F5'))
            expect_equal(sampleNames(bar1), c('S1', 'S2', 'S3', 'S4', 'S5'))
            expect_equal(nrow(bar1), c(Features=4))
            expect_equal(ncol(bar1), c(Samples=5))
            bar2 <- execute(cmdl2, mockMethylSet)
            expect_equal(featureNames(bar2), c('F2', 'F3', 'F4'))
            expect_equal(sampleNames(bar2), c('S1', 'S2', 'S3', 'S4', 'S5'))
            expect_equal(nrow(bar2), c(Features=3))
            expect_equal(ncol(bar2), c(Samples=5))
            bar3 <- execute(cmdl3, mockMethylSet)
            expect_equal(featureNames(bar3), c('F2', 'F3', 'F4'))
            expect_equal(sampleNames(bar3), c('S1', 'S3', 'S4'))
            expect_equal(nrow(bar3), c(Features=3))
            expect_equal(ncol(bar3), c(Samples=3))
            bar4 <- execute(cmdl4, mockMethylSet)
            expect_equal(featureNames(bar4), c('F2', 'F3', 'F4'))
            expect_equal(sampleNames(bar4), c('S1', 'S3', 'S4'))
            expect_equal(nrow(bar4), c(Features=3))
            expect_equal(ncol(bar4), c(Samples=3))
          })

