context('Probe-Gene relationship tests')

library(epilab)

# 
# A mock list of target ids
#
mockTids <- c('cg21268256', 'cg03472880', 'cg13583162', 'cg00079477', 'cg12302110', 'cg03046843',
              'cg07648498', 'cg23922289', 'cg26708970', 'cg27405400')
mockEntrez <- c('22864', '729041', '2175', '6092', '4000', '91107', '274')

#
# getGeneEntrezIds tests
#
test_that('getGeneEntrezIds fails on invalid tids input',
          {
            expect_error(getGeneEntrezIds(1:10, 'illumina', 2121))
            expect_error(getGeneEntrezIds(c(TRUE, FALSE, TRUE), 'ucsc19', 2121))
            expect_error(getGeneEntrezIds(factor(1:10), 'illumina', 98765))
          })

test_that('getGeneEntrezIds works correctly on example',
          {
            expect_error(getGeneEntrezIds(mockTids, 'illumina', 2000))
            foo <- getGeneEntrezIds(mockTids, 'ucsc19', 2000)
            expect_equal(foo, c('22864', '729041', '2175', '6092', '4000', '91107', '274'))
            foo <- getGeneEntrezIds(mockTids, 'ucsc19', 20000)
            expect_equal(foo, c('22864', '729041', '2175', '84501', '6092', '10444', '54662', 
                                '4000', '91107', '274'))
          })

test_that('getGeneEntrezIds works well on duplicated inputs',
          {
            doubleMockTids <- c(mockTids, mockTids)
            foo <- getGeneEntrezIds(mockTids, 'ucsc19', 3333)
            bar <- getGeneEntrezIds(doubleMockTids, 'ucsc19', 3333)
            expect_equal(foo, bar)
          })

#
# getSymbolsFromEntrezIds tests
#
test_that('getSymbolsFromEntrezIds fails on invalid tids input',
          {
            expect_error(getSymbolsFromEntrezIds(1:10))
            expect_error(getSymbolsFromEntrezIds(c(TRUE, FALSE, TRUE)))
            expect_error(getSymbolsFromEntrezIds(factor(1:10)))
          })

test_that('getSymbolsFromEntrezIds works on empty input',
          {
            foo <- getSymbolsFromEntrezIds(character(0))
            expect_equal(foo, character(0))
          })

test_that('getSymbolsFromEntrezIds works correctly on example',
          {
            foo <- getSymbolsFromEntrezIds(mockEntrez)
            expect_equivalent(foo, c('R3HDM2', 'FAAHP1', 'FANCA', 'ROBO2', 'LMNA', 'TRIM47', 'BIN1'))
          })

#
# getProbeGeneRelationship tests
#
test_that('getProbeGeneRelationship works on empty input',
          {
            expect_error(getProbeGeneRelationship(character(0), 'illumina', 2000))
            foo <- getProbeGeneRelationship(character(0), 'ucsc19', 3456)
            expect_equal(foo$probeId, character(0))
            expect_equal(foo$geneId, character(0))
            expect_equal(foo$geneSymbol, character(0))
          })

test_that('getProbeGeneRelationship fails on invalid promoterSize input',
          {
            expect_error(getProbeGeneRelationship(mockTids, 'illumina', -10))
            expect_error(getProbeGeneRelationship(mockTids, 'ucsc19', -10))
            expect_error(getProbeGeneRelationship(mockTids, 'ucsc19', 'foobar'))
            expect_error(getProbeGeneRelationship(mockTids, 'ucsc19', TRUE))
            expect_error(getProbeGeneRelationship(mockTids, 'ucsc19', c(100, 2000)))
          })

test_that('getProbeGeneRelationship fails on invalid tids input',
          {
            expect_error(getProbeGeneRelationship(1:10, 'illumina', 2121))
            expect_error(getProbeGeneRelationship(c(TRUE, FALSE, TRUE), 'ucsc19', 2121))
            expect_error(getProbeGeneRelationship(factor(1:10), 'illumina', 98765))
          })

test_that('getProbeGeneRelationship fails on unknown method',
          {
            expect_error(getProbeGeneRelationship(mockTids, 'illumina42', 2345))
          })

test_that('getProbeGeneRelationship works correctly on example',
          {
            expect_error(getProbeGeneRelationship(mockTids, 'illumina', 2000))
            foo <- getProbeGeneRelationship(mockTids, 'ucsc19', 2000)
            expect_equal(foo$probeId, c('cg00079477', 'cg03046843', 'cg07648498', 'cg12302110',
                                        'cg23922289', 'cg26708970', 'cg27405400'))
            expect_equal(foo$geneId, c('22864', '729041', '2175', '6092', '4000', '91107', '274'))
            expect_equal(foo$geneSymbol, c('R3HDM2', 'FAAHP1', 'FANCA', 'ROBO2', 'LMNA', 
                                           'TRIM47', 'BIN1'))
            foo <- getProbeGeneRelationship(mockTids, 'ucsc19', 20000)
            expect_equal(foo$probeId, c('cg00079477', 'cg03046843', 'cg07648498', 'cg07648498', 
                                        'cg12302110', 'cg13583162', 'cg13583162', 'cg23922289', 
                                        'cg26708970', 'cg27405400'))
            expect_equal(foo$geneId, c('22864', '729041', '2175', '84501', '6092', '10444', '54662', 
                                       '4000', '91107', '274'))
            expect_equal(foo$geneSymbol, c('R3HDM2', 'FAAHP1', 'FANCA', 'SPIRE2', 'ROBO2', 
                                           'ZER1', 'TBC1D13', 'LMNA', 'TRIM47', 'BIN1'))
          })

test_that('getProbeGeneRelationship works well on duplicated inputs',
          {
            doubleMockTids <- c(mockTids, mockTids)
            foo <- getProbeGeneRelationship(mockTids, 'ucsc19', 3333)
            bar <- getProbeGeneRelationship(doubleMockTids, 'ucsc19', 3333)
            expect_equal(foo$probeId, bar$probeId)
            expect_equal(foo$geneId, bar$geneId)
            expect_equal(foo$geneSymbol, bar$geneSymbol)
          })

