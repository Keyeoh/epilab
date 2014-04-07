context('Circos tests')

#
# Mock objects for testing
#
emptyRanges <- GRanges()
seqlengths <- c(249250621, 243199373, 198022430, 191154276, 180915260, 171115067, 159138663, 
                146364022, 141213431, 135534747, 135006516, 133851895, 115169878, 107349540, 
                102531392, 90354753, 81195210, 78077248, 59128983, 63025520, 48129895, 51304566, 
                155270560, 59373566)
names(seqlengths) <- paste0('chr', c(1:22, 'X', 'Y'))
mockRanges <- GRanges(seqnames=paste0('chr', c(6, 9, 11, 5, 7, 7, 13, 8, 19, 17)),
                      ranges=IRanges(start=c(1523751, 108209784, 1312565, 39074123, 91909243, 
                                             96655105, 26140724, 7404876, 14017309, 77833495),
                                     end=c(1523752, 108209785, 1312566, 39074124, 91909244, 
                                           96655106, 26140725, 7404877, 14017310, 77833496)),
                      strand='*', seqlengths=seqlengths)
names(mockRanges) <- c("cg13299743", "cg14341289", "cg07634195", "cg25020279", "cg19797536", 
                       "cg20250426", "cg12778938", "cg16762794", "cg26609550", "cg18502099")
mockRanges$dummy <- 1
mockRanges$notDummy <- 1:10

#
# averagePerBin tests
#

test_that('averagePerBin fails with wrong values of x',
          {
            expect_error(averagePerBin('foo', 999))
            expect_error(averagePerBin(1:17, 982))
            expect_error(averagePerBin(c(TRUE, FALSE), 666))
            expect_error(averagePerBin(matrix(rnorm(9), ncol=3), 654))
            expect_error(averagePerBin(emptyRanges, 777))
          })

test_that('averagePerBin fails with wrong values of binsize',
          {
            expect_error(averagePerBin(mockRanges, 0))
            expect_error(averagePerBin(mockRanges, -1))
            expect_error(averagePerBin(mockRanges, -Inf))
            expect_error(averagePerBin(mockRanges, Inf))
            expect_error(averagePerBin(mockRanges, NA))
          })

test_that('averagePerBin fails with wrong values of mcolnames',
          {
            expect_error(averagePerBin(mockRanges, 1e7, NA))
            expect_error(averagePerBin(mockRanges, 1e7, 1:10))
            expect_error(averagePerBin(mockRanges, 1e7, c(TRUE, FALSE)))
            expect_error(averagePerBin(mockRanges, 1e7, 'notSoDummy'))
            expect_error(averagePerBin(mockRanges, 1e7, c('dummy', 'aintHere', 'notDummy')))
          })

test_that('averagePerBin works correctly on a naive example',
          {
            foo <- averagePerBin(mockRanges, 1e7, c('dummy', 'notDummy'))
            expect_equal(length(foo), 322)
            fooDf <- as(foo, 'data.frame')
            bar <- aggregate(fooDf$dummy, list(fooDf$seqnames), length)
            expect_equal(bar$x, c(25, 25, 20, 20, 19, 18, 16, 15, 15, 14, 14, 14, 12, 11, 11, 10, 9,
                                  8, 6, 7, 5, 6, 16, 6))
            bar <- aggregate(fooDf$dummy, list(fooDf$seqnames), sum)
            expect_equal(bar$x, c(0, 0, 0, 0, 2e-07, 2e-07, 4e-07, 2e-07, 2e-07, 0, 2e-07, 0, 2e-07,
                                  0, 0, 0, 2e-07, 0, 2e-07, 0, 0, 0, 0, 0))
            bar <- aggregate(fooDf$notDummy, list(fooDf$seqnames), length)
            expect_equal(bar$x, c(25, 25, 20, 20, 19, 18, 16, 15, 15, 14, 14, 14, 12, 11, 11, 10, 9,
                                  8, 6, 7, 5, 6, 16, 6))
            bar <- aggregate(fooDf$notDummy, list(fooDf$seqnames), sum)
            expect_equal(bar$x, c(0, 0, 0, 0, 8e-07, 2e-07, 2.2e-06, 1.6e-06, 4e-07, 0, 6e-07, 0, 
                                  1.4e-06, 0, 0, 0, 2e-06, 0, 1.8e-06, 0, 0, 0, 0, 0))
          })

#
# generateCircosFromRanges tests
#

test_that('generateCircosFromRanges fails on wrong ranges input',
          {
            expect_error(generateCircosFromRanges('foo'))
            expect_error(generateCircosFromRanges(c(TRUE, FALSE)))
            expect_error(generateCircosFromRanges(1:17))
            expect_error(generateCircosFromRanges(matrix(rnorm(9), ncol=3)))
            expect_error(generateCircosFromRanges(emptyRanges))
          })

test_that('generateCircosFromRanges fails on wrong ids input',
          {
            expect_error(generateCircosFromRanges(mockRanges, 'foo'))
            expect_error(generateCircosFromRanges(mockRanges, c('cg13299743', 'foo')))
            expect_error(generateCircosFromRanges(mockRanges, NA))
            expect_error(generateCircosFromRanges(mockRanges, 7))
          })

test_that('generateCircosFromRanges fails on wrong values input',
          {
            expect_error(generateCircosFromRanges(mockRanges, values=letters[1:10]))
            expect_error(generateCircosFromRanges(mockRanges, values=matrix(rnorm(9), ncol=3)))
            expect_error(generateCircosFromRanges(mockRanges, values=1:11))
            expect_error(generateCircosFromRanges(mockRanges, values=1:7))
          })

test_that('generateCircosFromRanges works correctly on a naive example',
          {
            mockIds <- c('cg13299743', 'cg14341289', 'cg07634195', 'cg25020279')
            foo <- generateCircosFromRanges(mockRanges, ids=mockIds, values=1:4)
            expect_equal(rownames(foo), c('cg13299743', 'cg14341289', 'cg07634195', 'cg25020279'))
            expect_equal(foo$seqnames, c('hs6', 'hs9', 'hs11', 'hs5'))
            expect_equal(foo$start, c(1523751, 108209784, 1312565, 39074123))
            expect_equal(foo$end, c(1523752, 108209785, 1312566, 39074124))
            expect_equal(foo$values, 1:4)
          })    

