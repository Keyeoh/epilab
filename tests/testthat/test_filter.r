context('FilterCommand tests')

library(epilab)
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

mockPdata <- AnnotatedDataFrame(data.frame(foo=1:5))
mockFdata <- AnnotatedDataFrame(data.frame(bar=1:5))
mockMethylSet <- MethylSet(Meth=matrix(1:5, nrow=5, ncol=5), 
                           Unmeth=matrix(1:5, nrow=5, ncol=5, byrow=TRUE),
                           phenoData=mockPdata)
featureData(mockMethylSet) <- mockFdata
featureNames(mockMethylSet) <- paste0('F', 1:5)
sampleNames(mockMethylSet) <- paste0('S', 1:5)

mockPdataPlus <- AnnotatedDataFrame(data.frame(foo=1:6))
mockFdataPlus <- AnnotatedDataFrame(data.frame(bar=1:5))
mockMethylSetPlus <- MethylSet(Meth=matrix(1:6, nrow=5, ncol=6), 
                               Unmeth=matrix(1:6, nrow=5, ncol=6, byrow=TRUE),
                               phenoData=mockPdataPlus)
featureData(mockMethylSetPlus) <- mockFdataPlus
featureNames(mockMethylSetPlus) <- paste0('F', 1:5)
sampleNames(mockMethylSetPlus) <- paste0('S', 1:6)

mockBeta <- getBeta(mockMethylSet)

mockRanges <- GRanges(seqnames=c('chr1', 'chr7', 'chr3', 'chrX', 'chrY'),
                      ranges=IRanges(start=c(4, 8, 15, 16, 23),
                                     end=c(42, 50, 70, 32, 42)),
                      strand=c('+', '+', '-', '+', '-')
                      )
seqlengths(mockRanges) <- c(1000, 2000, 3000, 4000, 5000)

emptyGenomicSet <- GenomicMethylSet()

mockGenomicSet <- GenomicMethylSet(mockRanges, getMeth(mockMethylSet), getUnmeth(mockMethylSet), 
                                   pData(mockMethylSet), annotation(mockMethylSet), 
                                   preprocessMethod(mockMethylSet))
rownames(mockGenomicSet) <- paste0('F', 1:5)
colnames(mockGenomicSet) <- paste0('S', 1:5)

#
# FilterCommand tests
#
test_that('FilterCommand fails on objects without dimensions',
          {
            foo <- kOverAFilterCommand(mockDetectionP, byRow=TRUE, k=2, a=0.07)
            expect_error(execute(foo, 'foobar'))
            expect_error(execute(foo, list()))
          })
          
#
# AtomicFilterCommand tests
#
test_that('AtomicFilterCommand accessors work well on derived classes',
          {
            foo <- kOverAFilterCommand(mockDetectionP, byRow=TRUE, k=2, a=0.07)
            expect_equal(getByRow(foo), TRUE)
            setByRow(foo) <- FALSE
            expect_equal(getByRow(foo), FALSE)
            expect_error(setByRow(foo) <- 42)
            expect_error(setByRow(foo) <- 'foobar')
            expect_error(setByRow(foo) <- c(TRUE, FALSE, TRUE))
          })

#
# MatrixFilterCommand tests
#
test_that('MatrixFilterCommand accessors work well on derived classes',
          {
            foo <- kOverAFilterCommand(mockDetectionP, byRow=TRUE, k=2, a=0.07)
            expect_equal(getMatrix(foo), mockDetectionP)
            setMatrix(foo) <- zeroDetectionP
            expect_equal(getMatrix(foo), zeroDetectionP)
            expect_error(setMatrix(foo) <- 42)
            expect_error(setMatrix(foo) <- 'foobar')
            expect_error(setMatrix(foo) <- TRUE)
          })

#
# KOverAFilterCommand tests
#
test_that('KOverAFilterCommand specific accessors work well',
          {
            foo <- kOverAFilterCommand(mockDetectionP, byRow=TRUE, k=2, a=0.07)
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

test_that('KOverAfilterCommand percentage constructor works as expected',
          {
            foo <- kOverAFilterCommandFromFraction(mockDetectionP, byRow=TRUE, fraction=0.4, 
                                                       a=0.07)
            expect_equal(getK(foo), 2)
            expect_equal(getA(foo), 0.07)
            expect_error(foo <- kOverAFilterCommandFromFraction(mockDetectionP, byRow=TRUE, 
                                                                    fraction=-0.01, a=0.07))
            expect_error(foo <- kOverAFilterCommandFromFraction(mockDetectionP, byRow=TRUE, 
                                                                    fraction=1.01, a=0.07))
            expect_error(foo <- kOverAFilterCommandFromFraction(mockDetectionP, byRow=TRUE, 
                                                                    fraction=c(0.5, 0.7), a=0.07))

          })

koveracmd1 <- kOverAFilterCommand(mockDetectionP, byRow=TRUE, k=2, a=0.5)
koveracmd2 <- kOverAFilterCommand(mockDetectionP, byRow=TRUE, k=2, a=0.01)
koveracmd3 <- kOverAFilterCommand(mockDetectionP, byRow=FALSE, k=1, a=0.5)
koveracmd4 <- kOverAFilterCommand(mockDetectionP, byRow=FALSE, k=1, a=0.01)

test_that('KOverAFilterCommand execution breaks on wrong data types',
          {
            expect_error(execute(mockMethylSet, koveracmd1))
            expect_error(execute(mockMethylSet, koveracmd2))
            expect_error(execute(mockGenomicSet, koveracmd3))
            expect_error(execute(mockMethylSet, koveracmd4))
          })

test_that('KOverAFilterCommand fails on empty example',
          {
            expect_error(execute(koveracmd1, emptyGenomicSet))
            expect_error(execute(koveracmd3, emptyGenomicSet))
          })

test_that('KOverAFilterCommand execution works correctly on examples',
          {
            bar1 <- execute(koveracmd1, mockMethylSet)
            expect_equal(featureNames(bar1), c('F2', 'F3', 'F4', 'F5'))
            expect_equal(sampleNames(bar1), c('S1', 'S2', 'S3', 'S4', 'S5'))
            expect_equal(nrow(bar1), c(Features=4))
            expect_equal(ncol(bar1), c(Samples=5))
            bar2 <- execute(koveracmd2, mockMethylSet)
            expect_equal(featureNames(bar2), c('F2', 'F3', 'F4'))
            expect_equal(sampleNames(bar2), c('S1', 'S2', 'S3', 'S4', 'S5'))
            expect_equal(nrow(bar2), c(Features=3))
            expect_equal(ncol(bar2), c(Samples=5))
            bar3 <- execute(koveracmd3, mockMethylSet)
            expect_equal(featureNames(bar3), c('F1', 'F2', 'F3', 'F4', 'F5'))
            expect_equal(sampleNames(bar3), c('S3', 'S4'))
            expect_equal(nrow(bar3), c(Features=5))
            expect_equal(ncol(bar3), c(Samples=2))
            bar4 <- execute(koveracmd4, mockGenomicSet)
            expect_equal(rownames(bar4), c('F1', 'F2', 'F3', 'F4', 'F5'))
            expect_equal(colnames(bar4), character(0))
            expect_equal(nrow(bar4), 5)
            expect_equal(ncol(bar4), 0)
          })

test_that('KOverAFilterCommand fails on different dimension names',
          {
            expect_error(execute(koveracmd1, mockMethylSetPlus),
                         regexp='must be included in dimension names')
          })

#
# VarFilterCommand tests
#
test_that('VarFilterCommand specific accessors work well',
          {
            foo <- varFilterCommand(mockBeta, byRow=TRUE, type='quantile', threshold=0.25)
            expect_equal(getType(foo), 'quantile')
            expect_equal(getThreshold(foo), 0.25)
            setType(foo) <- 'absolute'
            expect_equal(getType(foo), 'absolute')
            setThreshold(foo) <- 0.00007
            expect_equal(getThreshold(foo), 0.00007)
            expect_error(setType(foo) <- -1)
            expect_error(setType(foo) <- 9999)
            expect_error(setType(foo) <- c(4, 8, 15, 16, 23, 42))
            expect_error(setType(foo) <- 'foobar')
            expect_error(setThreshold(foo) <- -1)
            setType(foo) <- 'quantile'
            expect_error(setThreshold(foo) <- 9999)
            expect_error(setThreshold(foo) <- c(4, 8, 15, 16, 23, 42))
            expect_error(setThreshold(foo) <- 'foobar')
            setType(foo) <- 'quantile'
            expect_error(setThreshold(foo) <- -0.000001)
            expect_error(setThreshold(foo) <- 1.0000001)
          })

varcmd1 <- varFilterCommand(mockBeta, byRow=TRUE, type='quantile', threshold=0.25)
varcmd2 <- varFilterCommand(mockBeta, byRow=FALSE, type='absolute', threshold=0.019)

test_that('VarFilterCommand execution breaks on wrong data types',
          {
            expect_error(execute(mockMethylSet, varcmd1))
            expect_error(execute(mockGenomicSet, varcmd2))
          })

test_that('VarFilterCommand fails on empty example',
          {
            expect_error(execute(varcmd2, emptyGenomicSet))
          })

test_that('VarFilterCommand execution works correctly on examples',
          {
            bar1 <- execute(varcmd1, mockGenomicSet)
            expect_equal(rownames(bar1), c('F1', 'F2', 'F3', 'F4'))
            expect_equal(colnames(bar1), c('S1', 'S2', 'S3', 'S4', 'S5'))
            expect_equal(nrow(bar1), 4)
            expect_equal(ncol(bar1), 5)
            bar2 <- execute(varcmd2, mockMethylSet)
            expect_equal(featureNames(bar2), c('F1', 'F2', 'F3', 'F4', 'F5'))
            expect_equal(sampleNames(bar2), c('S2', 'S3', 'S4'))
            expect_equal(nrow(bar2), c(Features=5))
            expect_equal(ncol(bar2), c(Samples=3))
          })

test_that('VarFilterCommand fails on different dimension names',
          {
            expect_error(execute(varcmd1, mockMethylSetPlus),
                         regexp='must be included in dimension names')
          })

#
# FilterCommandList tests
#
test_that('FilterCommandList gets its slots right', 
          {
            foo <- filterCommandList(koveracmd1, koveracmd2, koveracmd3, koveracmd4)
            koveracmdList <- list(koveracmd1, koveracmd2, koveracmd3, koveracmd4)
            expect_equal(getCommandList(foo), koveracmdList)
          })

test_that('FilterCommandList refuses wrong data types', 
          {
            expect_error(filterCommandList(2032))
          })

cmdl1 <- filterCommandList(koveracmd1)
cmdl2 <- filterCommandList(koveracmd1, koveracmd2)
cmdl3 <- filterCommandList(koveracmd1, koveracmd2, koveracmd3)
cmdl4 <- filterCommandList(koveracmd1, koveracmd2, koveracmd3, koveracmd4)

test_that('FilterCommandList execution breaks on wrong data types', 
          {
            expect_error(execute(mockMethylSet, cmdl1))
          })

test_that('FilterCommandList execution works correctly on example', 
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
            bar3 <- execute(cmdl3, mockGenomicSet)
            expect_equal(rownames(bar3), c('F2', 'F3', 'F4'))
            expect_equal(colnames(bar3), c('S1', 'S3', 'S4'))
            expect_equal(nrow(bar3), 3)
            expect_equal(ncol(bar3), 3)
            bar4 <- execute(cmdl4, mockMethylSet)
            expect_equal(featureNames(bar4), c('F2', 'F3', 'F4'))
            expect_equal(sampleNames(bar4), c('S1', 'S3', 'S4'))
            expect_equal(nrow(bar4), c(Features=3))
            expect_equal(ncol(bar4), c(Samples=3))
          })

