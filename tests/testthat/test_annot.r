context('AnnotationCommand tests')

#
# Mock objects for testing
#
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

#
# DensCpGCommand tests
#
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
test_that('DensCpGCommand execution breaks on wrong data types', {
          expect_error(execute(mockRanges, cmd))
})

test_that('DensCpGCommand execution works correctly on example', {
          bar <- execute(cmd, mockRanges)
          expect_equal(bar$foo, c(0.083, 0.102, 0.085, 0.158, 0.011, 0.063, 
                                  0.013, 0.065, 0.112, 0.056))
})

test_that('DensCpGCommand execution fails on empty example', {
          expect_error(execute(cmd, emptyRanges))
})

#
# CPGICommand tests
#
test_that('CPGICommand gets its slots right', {
          foo <- CPGICommand('bar', TRUE)
          expect_equal(foo@colName, 'bar')
          expect_equal(foo@discardDirection, TRUE)
})

test_that('CPGICommand refuses wrong data types', {
          expect_error(CPGICommand(2032, TRUE))
          expect_error(CPGICommand('bar', 'foo'))
          expect_error(CPGICommand(-9, 666))
})

cmdFalse <- CPGICommand('foo', FALSE)
cmdTrue <- CPGICommand('foo', TRUE)

test_that('CPGICommand execution breaks on wrong data types', {
          expect_error(execute(mockRanges, cmdTrue))
          expect_error(execute(mockRanges, cmdFalse))
})

test_that('CPGICommand execution works correctly on example', {
          bar <- execute(cmdTrue, mockRanges)
          expect_equal(bar$foo, c('CGI', 'CGI-Shore', 'CGI', 'CGI', 'Non-CGI',
                                  'CGI-Shore', 'Non-CGI', 'CGI', 'CGI',
                                  'CGI-Shore'))
          bar <- execute(cmdFalse, mockRanges)
          expect_equal(bar$foo, c('CGI', 'CGI-N-Shore', 'CGI', 'CGI', 'Non-CGI',
                                  'CGI-S-Shore', 'Non-CGI', 'CGI', 'CGI',
                                  'CGI-N-Shore'))
})

test_that('CPGICommand execution fails on empty example', {
          expect_error(execute(cmd, emptyRanges))
})

#
# GapCommand tests
#
test_that('GapCommand gets its slots right', {
          foo <- GapCommand('bar')
          expect_equal(foo@colName, 'bar')
})

test_that('GapCommand refuses wrong data types', {
          expect_error(GapCommand(2032))
})

cmd <- GapCommand('foo')

test_that('GapCommand execution breaks on wrong data types', {
          expect_error(execute(mockRanges, cmd))
})

test_that('GapCommand execution works correctly on example', {
          bar <- execute(cmd, mockRanges)
          expect_equal(bar$fooCent, c(57306413, 57842104, 50331638, 7331516, 
                                      30854911, 35600773, 7140723, 36434009, 
                                      10664471, 52570488))
          expect_equal(bar$fooTelo, c(1513750, 32993645, 1302564, 39064122, 
                                      67219418, 62473556, 26130723, 7394875, 
                                      14007308, NA))
})

test_that('GapCommand execution fails on empty example', {
          expect_error(execute(cmd, emptyRanges))
})

#
# GenomicRegionCommand tests
#
test_that('GenomicRegionCommand gets its slots right', {
          foo <- GenomicRegionCommand('bar')
          expect_equal(foo@colName, 'bar')
})

test_that('GenomicRegionCommand refuses wrong data types', {
          expect_error(GenomicRegionCommand(2032))
})

cmd <- GenomicRegionCommand('foo')

test_that('GenomicRegionCommand execution breaks on wrong data types', {
          expect_error(execute(mockRanges, cmd))
})

test_that('GenomicRegionCommand execution works correctly on example', {
          bar <- execute(cmd, mockRanges)
          expect_equal(bar$fooProm, c(FALSE, TRUE, FALSE, FALSE, FALSE, TRUE, 
                                      TRUE, FALSE, TRUE, FALSE))
          expect_equal(bar$fooIntra, c(FALSE, FALSE, TRUE, TRUE, TRUE, FALSE, 
                                       TRUE, FALSE, FALSE, FALSE))
          expect_equal(bar$fooInter, c(TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, 
                                       FALSE, TRUE, FALSE, TRUE))
})

test_that('GenomicRegionCommand execution fails on empty example', {
          expect_error(execute(cmd, emptyRanges))
})

