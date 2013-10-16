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
          expect_error(DensCpGCommand('foo', NA))
          expect_error(DensCpGCommand('foo', NaN))
})

cmd <- DensCpGCommand('foo', 2000)
emptyRanges <- GRanges()
seqlengths <- c(249250621, 243199373, 198022430, 191154276, 180915260, 
                171115067, 159138663, 146364022, 141213431, 135534747, 
                135006516, 133851895, 115169878, 107349540, 102531392, 90354753, 
                81195210, 78077248, 59128983, 63025520, 48129895, 51304566, 
                155270560, 59373566)
names(seqlengths) <- paste0('chr', c(1:22, 'X', 'Y'))
mockRanges <- GRanges(seqnames=paste0('chr', c(6, 9, 11, 5, 7, 7, 13, 8, 19, 
                                               17)),
                      ranges=IRanges(start=c(1523751, 108209784, 1312565, 
                                             39074123, 91909243, 96655105,
                                             26140724, 7404876, 14017309,
                                             77833495),
                                     end=c(1523752, 108209785, 1312566, 
                                           39074124, 91909244, 96655106, 
                                           26140725, 7404877, 14017310, 
                                           77833496)),
                      strand='*',
                      seqlengths=seqlengths)
names(mockRanges) <- c("cg13299743", "cg14341289", "cg07634195", "cg25020279", 
                       "cg19797536", "cg20250426", "cg12778938", "cg16762794", 
                       "cg26609550", "cg18502099")

test_that('DensCpGCommand execution breaks on wrong data types', {
          expect_error(execute(mockRanges, cmd))
})

test_that('DensCpGCommand execution works correctly on example', {
          bar <- execute(cmd, mockRanges)
          expect_equal(bar$foo, c(0.083, 0.102, 0.085, 0.158, 0.011, 0.063, 
                                  0.013, 0.065, 0.112, 0.056))
})

          


