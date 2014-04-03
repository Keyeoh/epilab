context('Graphs tests')

#
# barGroupGraph tests
#

test_that('barGroupGraph fails with wrong dataList',
          {
            expect_error(barGroupGraph('foo'))
            expect_error(barGroupGraph(1:10))
            expect_error(barGroupGraph(TRUE))
            expect_error(barGroupGraph(matrix(1:10, ncol=2)))
          })

test_that('barGroupGraph works correctly on a naive example',
          {
            mockDataList <- list(g1=c('A', 'B', 'C'),
                                 g2=factor(c('A', 'A')),
                                 g3=c('B', 'C'))
            mockDataFrame <- data.frame(group=c('g1', 'g1', 'g1', 'g2', 'g2', 'g3', 'g3'),
                                        status=c('A', 'B', 'C', 'A', 'A', 'B', 'C'))
            foo <- barGroupGraph(mockDataList)
            expect_true(is(foo, 'gg'))
            expect_true(is(foo, 'ggplot'))
            expect_equal(length(names(foo)), 9)
            expect_equal(names(foo), c('data', 'layers', 'scales', 'mapping', 'theme', 
                                       'coordinates', 'facet', 'plot_env', 'labels'))
            expect_equal(foo$data, mockDataFrame)
            expect_equal(as.character(foo$mapping$x), 'group')
            expect_equal(as.character(foo$mapping$fill), 'status')
            expect_equal(foo$scales$scales[[1]]$scale_name, 'grey')
          })

#
# barGroupGraphFromMatrix tests
#

test_that('barGroupGraphFromMatrix fails with wrong matrix',
          {
            expect_error(barGroupGraphFromMatrix('foo'))
            expect_error(barGroupGraphFromMatrix(1:10))
            expect_error(barGroupGraphFromMatrix(TRUE))
            expect_error(barGroupGraphFromMatrix(list('a', 1, TRUE)))
          })

test_that('barGroupGraphFromMatrix works correctly on a naive example',
          {
            mockDataMatrix <- matrix(c(1, 1, 1, 2, 0, 0, 0, 1, 1), ncol=3)
            rownames(mockDataMatrix) <- c('A', 'B', 'C')
            colnames(mockDataMatrix) <- c('g1', 'g2', 'g3')
            mockDataFrame <- 
              data.frame(status=c('A', 'B', 'C', 'A', 'B', 'C', 'A', 'B', 'C'),
                         group=c('g1', 'g1', 'g1', 'g2', 'g2', 'g2', 'g3', 'g3', 'g3'),
                         value=c(1, 1, 1, 2, 0, 0, 0, 1, 1))
            foo <- barGroupGraphFromMatrix(mockDataMatrix)
            expect_true(is(foo, 'ggplot'))
            expect_equal(length(names(foo)), 9)
            expect_equal(names(foo), c('data', 'layers', 'scales', 'mapping', 'theme', 
                                       'coordinates', 'facet', 'plot_env', 'labels'))
            expect_equal(foo$data, mockDataFrame)
            expect_equal(as.character(foo$mapping$x), 'group')
            expect_equal(as.character(foo$mapping$fill), 'status')
            expect_equal(foo$scales$scales[[1]]$scale_name, 'grey')
          })

