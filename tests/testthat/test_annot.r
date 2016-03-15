context('AnnotationCommand tests')

#
# Mock objects for testing
#
emptyRanges <- GRanges()
seqlengths <- c(249250621, 243199373, 198022430, 191154276, 180915260, 171115067, 159138663,
                146364022, 141213431, 135534747, 135006516, 133851895, 115169878, 107349540,
                102531392, 90354753, 81195210, 78077248, 59128983, 63025520, 48129895, 51304566,
                155270560, 59373566)
names(seqlengths) <- paste0('chr', c(1:22, 'X', 'Y'))
mockRanges <- GRanges(seqnames=paste0('chr', c(6, 9, 11, 5, 7, 7, 13, 8, 19, 17,
                                               6, 5, 20, 6, 4)),
                      ranges=IRanges(start=c(1523751, 108209784, 1312565, 39074123, 91909243,
                                             96655105, 26140724, 7404876, 14017309, 77833495,
                                             42893443, 1122815, 61866960, 166418799, 866811),
                                     end=c(1523752, 108209785, 1312566, 39074124, 91909244,
                                           96655106, 26140725, 7404877, 14017310, 77833496,
                                           42893444, 1122816, 61866961, 166418800, 866812)),
                      strand='*', seqlengths=seqlengths)
names(mockRanges) <- c("cg13299743", "cg14341289", "cg07634195", "cg25020279", "cg19797536",
                       "cg20250426", "cg12778938", "cg16762794", "cg26609550", "cg18502099",
                       "cg16392865", "cg00395291", "cg09310185", "cg21749424", "cg03307425")

#
# DensCpGCommand tests
#
test_that('DensCpGCommand gets its slots right',
          {
            foo <- densCpGCommand('bar', 2032)
            expect_equal(getColName(foo), 'bar')
            expect_equal(getWindowSize(foo), 2032)
          })

test_that('DensCpGCommand refuses wrong data types',
          {
            expect_error(densCpGCommand(2032, 2032))
            expect_error(densCpGCommand('bar', 'foo'))
          })

test_that('DensCpGCommand refuses invalid windowSize',
          {
            expect_error(densCpGCommand('foo', -1))
            expect_error(densCpGCommand('foo', 0))
            expect_error(densCpGCommand('foo', 1))
            expect_error(densCpGCommand('foo', NA))
            expect_error(densCpGCommand('foo', NaN))
          })

cmd <- densCpGCommand('foo', 2000)

test_that('DensCpGCommand execution breaks on wrong data types',
          {
            expect_error(execute(mockRanges, cmd))
          })

test_that('DensCpGCommand execution works correctly on example',
          {
            bar <- execute(cmd, mockRanges)
            expect_equal(bar$foo, c(0.083, 0.102, 0.085, 0.158, 0.011, 0.063, 0.013, 0.065, 0.112,
                                    0.056, 0.048, 0.072, 0.06, 0.079, 0.077))
          })

test_that('DensCpGCommand execution fails on empty example',
          {
            expect_error(execute(cmd, emptyRanges))
          })

#
# CPGICommand tests
#
test_that('CPGICommand gets its slots right',
          {
            foo <- cpgiCommand('bar', TRUE)
            expect_equal(getColName(foo), 'bar')
            expect_equal(getDiscardDirection(foo), TRUE)
          })

test_that('CPGICommand refuses wrong data types',
          {
            expect_error(cpgiCommand(2032, TRUE))
            expect_error(cpgiCommand('bar', 'foo'))
            expect_error(cpgiCommand(-9, 666))
          })

cmdFalse <- cpgiCommand('foo', FALSE)
cmdTrue <- cpgiCommand('foo', TRUE)

test_that('CPGICommand execution breaks on wrong data types',
          {
            expect_error(execute(mockRanges, cmdTrue))
            expect_error(execute(mockRanges, cmdFalse))
          })

test_that('CPGICommand execution works correctly on example',
          {
            bar <- execute(cmdTrue, mockRanges)
            expect_equal(bar$foo, c('CGI', 'CGI-Shore', 'CGI', 'CGI', 'Non-CGI', 'CGI-Shore',
                                    'Non-CGI', 'CGI', 'CGI', 'CGI-Shore', 'CGI-Shore', 'CGI',
                                    'CGI-Shelf', 'CGI', 'CGI-Shore'))
            bar <- execute(cmdFalse, mockRanges)
            expect_equal(bar$foo, c('CGI', 'CGI-N-Shore', 'CGI', 'CGI', 'Non-CGI', 'CGI-S-Shore',
                                    'Non-CGI', 'CGI', 'CGI', 'CGI-N-Shore', 'CGI-S-Shore', 'CGI',
                                    'CGI-N-Shelf', 'CGI', 'CGI-S-Shore'))
          })

test_that('CPGICommand execution fails on empty example',
          {
            expect_error(execute(cmd, emptyRanges))
          })

#
# GapCommand tests
#
test_that('GapCommand gets its slots right',
          {
            foo <- gapCommand('bar')
            expect_equal(getColName(foo), 'bar')
          })

test_that('GapCommand refuses wrong data types',
          {
            expect_error(gapCommand(2032))
          })

cmd <- gapCommand('foo')

test_that('GapCommand execution breaks on wrong data types',
          {
            expect_error(execute(mockRanges, cmd))
          })

test_that('GapCommand execution works correctly on example',
          {
            bar <- execute(cmd, mockRanges)
            expect_equal(bar$fooCent, NULL)
            expect_equal(bar$fooTelo, c(1513750, 30174932, 1302564, 39064122, 67426729, 62680867,
                                        26130723, 7394875, 14007308, 5413945, 42883442, 1112814,
                                        2567206, 4377179, 856810))
          })

test_that('GapCommand execution fails on empty example',
          {
            expect_error(execute(cmd, emptyRanges))
          })

#
# GenomicRegionCommand tests
#
test_that('GenomicRegionCommand gets its slots right',
          {
            foo <- genomicRegionCommand('bar')
            expect_equal(getColName(foo), 'bar')
          })

test_that('GenomicRegionCommand refuses wrong data types',
          {
            expect_error(genomicRegionCommand(2032))
          })

cmd <- genomicRegionCommand('foo')

test_that('GenomicRegionCommand execution breaks on wrong data types',
          {
            expect_error(execute(mockRanges, cmd))
          })

test_that('GenomicRegionCommand execution works correctly on example',
          {
            bar <- execute(cmd, mockRanges)

            expect_equivalent(bar$fooProm, c(FALSE, TRUE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE,
                                             TRUE, FALSE, FALSE, FALSE, TRUE, FALSE, TRUE))
            expect_equivalent(bar$fooIntra, c(FALSE, FALSE, TRUE, TRUE, TRUE, FALSE, TRUE, FALSE,
                                              FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, TRUE))
            expect_equivalent(bar$fooInter, c(TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE,
                                              TRUE, FALSE, TRUE, FALSE, TRUE, FALSE, TRUE, FALSE))
          })

#
# DGenomicRegionCommand tests
#
test_that('DGenomicRegionCommand gets its slots right',
          {
            foo <- dGenomicRegionCommand('bar')
            expect_equal(getColName(foo), 'bar')
          })

test_that('DGenomicRegionCommand refuses wrong data types',
          {
            expect_error(dGenomicRegionCommand(2032))
          })

cmd <- dGenomicRegionCommand('foo')

test_that('DGenomicRegionCommand execution breaks on wrong data types',
          {
            expect_error(execute(mockRanges, cmd))
          })

test_that('DGenomicRegionCommand execution works correctly on example',
          {
            bar <- execute(cmd, mockRanges)

            expectedValues <- c('Intergenic', 'Promoter', 'Intron', 'Intron', 'Intron', 'Promoter',
                                'Promoter', 'Intergenic', 'Promoter', 'Intergenic', '3UTR',
                                'Intergenic', 'Promoter', 'Intergenic', 'Promoter')
            expect_equivalent(bar$foo, expectedValues)
          })

test_that('DGenomicRegionCommand execution fails on empty example',
          {
            expect_error(execute(cmd, emptyRanges))
          })

#
# NearestGeneCommand tests
#
# nearestGeneCommand
test_that('NearestGeneCommand gets its slots right',
          {
            foo <- nearestGeneCommand('bar')
            expect_equal(getColName(foo), 'bar')
          })

test_that('NearestGeneCommand refuses wrong data types',
          {
            expect_error(nearestGeneCommand(2032))
          })

cmd <- nearestGeneCommand('foo')

test_that('NearestGeneCommand execution breaks on wrong data types',
          {
            expect_error(execute(mockRanges, cmd))
          })

test_that('NearestGeneCommand execution works correctly on example',
          {
            bar <- execute(cmd, mockRanges)
            expect_equivalent(bar$fooGeneSymbol, c('FOXC1', 'FSD1L', 'TOLLIP', 'RICTOR', 'ANKIB1',
                                                   'DLX5', 'ATP8A2', 'FAM90A7P', 'CC2D1A', 'CBX4',
                                                   'PTCRA', 'SLC12A7', 'BIRC7', 'LINC00602', 'GAK'))
            expect_equivalent(bar$fooGeneId, c(2296, 83856, 54472, 253260, 54467, 1749, 51761,
                                               441317, 54862, 8535,
                                               171558, 10723, 79444, 441177, 2580))
#             expect_equivalent(bar$fooDTSS, c(-86928, -529, 18325, 376, 33694, -961, 194514, 27042,
#                                              352, -20281))
            expect_equivalent(bar$fooDistance, c(-86928, -529, 0, 0, 0, -961, 0, 8782, 0, -20281,
                                                 0, -10642, -314, 15695, 0))
          })

test_that('NearestGeneCommand execution fails on empty example',
          {
            expect_error(execute(cmd, emptyRanges))
          })

# nearestTSSGenCommand
test_that('NearestGeneCommand gets its slots right',
          {
            foo <- nearestTSSGeneCommand('bar')
            expect_equal(getColName(foo), 'bar')
          })

test_that('NearestGeneCommand refuses wrong data types',
          {
            expect_error(nearestTSSGeneCommand(2032))
          })

cmd <- nearestTSSGeneCommand('foo')

test_that('NearestGeneCommand execution breaks on wrong data types',
          {
            expect_error(execute(mockRanges, cmd))
          })

test_that('NearestGeneCommand execution works correctly on example',
          {
            bar <- execute(cmd, mockRanges)
            expect_equivalent(bar$fooGeneSymbol, c('FOXC1', 'FSD1L', 'TOLLIP', 'RICTOR', 'ANKIB1',
                                                   'DLX5', 'ATP8A2', 'FAM90A7P', 'CC2D1A', 'CBX4',
                                                   'CNPY3', 'SLC12A7', 'BIRC7', 'LINC00473', 'GAK'))
            expect_equivalent(bar$fooGeneId, c(2296, 83856, 54472, 253260, 54467, 1749, 51761,
                                               441317, 54862, 8535, 10695, 10723, 79444, 90632, 2580))
            expect_equivalent(bar$fooDistance, c(-86928, -529, 11406, 376, 33694, -961, -545, 27042,
                                                 352, -20281, -3415, -10642, -314, -17271, -249))
          })

test_that('NearestGeneCommand execution fails on empty example',
          {
            expect_error(execute(cmd, emptyRanges))
          })

# nearestTXGeneCommand
test_that('NearestGeneCommand gets its slots right',
          {
            foo <- nearestTXGeneCommand('bar')
            expect_equal(getColName(foo), 'bar')
          })

test_that('NearestGeneCommand refuses wrong data types',
          {
            expect_error(nearestTXGeneCommand(2032))
          })

cmd <- nearestTXGeneCommand('foo')

test_that('NearestGeneCommand execution breaks on wrong data types',
          {
            expect_error(execute(mockRanges, cmd))
          })

test_that('NearestGeneCommand execution works correctly on example',
          {
            bar <- execute(cmd, mockRanges)
            expect_equivalent(bar$fooGeneSymbol, c('FOXC1', 'FSD1L', 'TOLLIP', 'RICTOR', 'ANKIB1',
                                                   'DLX5', 'ATP8A2', 'FAM90A7P', 'CC2D1A', 'CBX4',
                                                   'PTCRA', 'SLC12A7', 'BIRC7', 'LINC00602', 'GAK'))
            expect_equivalent(bar$fooGeneId, c(2296, 83856, 54472, 253260, 54467, 1749, 51761,
                                               441317, 54862, 8535, 171558, 10723, 79444, 441177,
                                               2580))
            expect_equivalent(bar$fooDistance, c(-86928, -529, 0, 0, 0, -961, 0, 8782, 0, -20281, 0,
                                                 -10642, -314, 15695, 0))
          })

test_that('NearestGeneCommand execution fails on empty example',
          {
            expect_error(execute(cmd, emptyRanges))
          })

#
# AnnotationCommandList tests
#
cmdList <- list(nearestGeneCommand('foo1'), nearestGeneCommand('foo2'))

test_that('AnnotationCommandList gets its slots right',
          {
            foo <- annotationCommandList(cmdList[[1]], cmdList[[2]])
            expect_equal(getColName(foo), '')
            expect_equal(getCommandList(foo), cmdList)
          })

test_that('AnnotationCommandList refuses wrong data types',
          {
            expect_error(annotationCommandList(2032))
          })

cmd <- annotationCommandList(densCpGCommand('dcpg'),
                             cpgiCommand('cpgi'),
                             gapCommand('gap'),
                             genomicRegionCommand('genreg'),
                             nearestTSSGeneCommand('ng')
                             )

test_that('AnnotationCommandList prevents colName slot from being accidentally changed',
          {
            expect_error(setColName(cmd) <- 'foo')
          })

test_that('AnnotationCommandList execution breaks on wrong data types',
          {
            expect_error(execute(mockRanges, cmd))
          })

test_that('AnnotationCommandList execution works correctly on example',
          {
            bar <- execute(cmd, mockRanges)
            expect_equal(bar$dcpg, c(0.083, 0.102, 0.085, 0.158, 0.011, 0.063, 0.013, 0.065, 0.112,
                                     0.056, 0.048, 0.072, 0.06, 0.079, 0.077))
            expect_equal(bar$cpgi, c('CGI', 'CGI-N-Shore', 'CGI', 'CGI', 'Non-CGI', 'CGI-S-Shore',
                                     'Non-CGI', 'CGI', 'CGI', 'CGI-N-Shore', 'CGI-S-Shore', 'CGI',
                                     'CGI-N-Shelf', 'CGI', 'CGI-S-Shore'))
            expect_equal(bar$gapCent, NULL)
            expect_equal(bar$gapTelo, c(1513750, 30174932, 1302564, 39064122, 67426729, 62680867,
                                        26130723, 7394875, 14007308, 5413945, 42883442, 1112814,
                                        2567206, 4377179, 856810))
            expect_equivalent(bar$genregProm, c(FALSE, TRUE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE,
                                                TRUE, FALSE, FALSE, FALSE, TRUE, FALSE, TRUE))
            expect_equivalent(bar$genregIntra, c(FALSE, FALSE, TRUE, TRUE, TRUE, FALSE, TRUE, FALSE,
                                                 FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, TRUE))
            expect_equivalent(bar$genregInter, c(TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE,
                                                 TRUE, FALSE, TRUE, FALSE, TRUE, FALSE, TRUE, FALSE))
            expect_equivalent(bar$ngGeneSymbol, c('FOXC1', 'FSD1L', 'TOLLIP', 'RICTOR', 'ANKIB1',
                                                  'DLX5', 'ATP8A2', 'FAM90A7P', 'CC2D1A', 'CBX4',
                                                  'CNPY3', 'SLC12A7', 'BIRC7', 'LINC00473', 'GAK'))
            expect_equivalent(bar$ngGeneId, c(2296, 83856, 54472, 253260, 54467, 1749, 51761,
                                              441317, 54862, 8535, 10695, 10723, 79444, 90632, 2580))
            expect_equivalent(bar$ngDistance, c(-86928, -529, 11406, 376, 33694, -961, -545, 27042,
                                                352, -20281, -3415, -10642, -314, -17271, -249))
          })

test_that('AnnotationCommandList execution fails on empty example',
          {
            expect_error(execute(cmd, emptyRanges))
          })


