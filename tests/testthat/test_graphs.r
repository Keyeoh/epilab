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
            wrongDataList <- list(g1=c('A', 'B', 'C'),
                                  g2=factor(c('A', 'A')),
                                  g3=c(4, 8, 15, 16, 23, 42))
            expect_error(barGroupGraph(wrongDataList))
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


#
# violinGroupGraph tests
#

test_that('violinGroupGraph fails with wrong dataList',
          {
            expect_error(violinGroupGraph('foo'))
            expect_error(violinGroupGraph(1:10))
            expect_error(violinGroupGraph(TRUE))
            expect_error(violinGroupGraph(matrix(1:10, ncol=2)))
            wrongDataList <- list(g1=c('A', 'B', 'C'),
                                  g2=factor(c('A', 'A')),
                                  g3=c(4, 8, 15, 16, 23, 42))
            expect_error(violinGroupGraph(wrongDataList))
          })

test_that('violinGroupGraph works correctly on a naive example',
          {
            mockDataList <- list(g1=1:5,
                                 g2=2:4,
                                 g3=c(4, 8, 15))
            mockDataFrame <- data.frame(group=c('g1', 'g1', 'g1', 'g1', 'g1', 'g2', 'g2', 'g2',
                                                'g3', 'g3', 'g3'),
                                        value=c(1, 2, 3, 4, 5, 2, 3, 4, 4, 8, 15))
            foo <- violinGroupGraph(mockDataList)
            expect_true(is(foo, 'gg'))
            expect_true(is(foo, 'ggplot'))
            expect_equal(length(names(foo)), 9)
            expect_equal(names(foo), c('data', 'layers', 'scales', 'mapping', 'theme', 
                                       'coordinates', 'facet', 'plot_env', 'labels'))
            expect_equal(foo$data, mockDataFrame)
            expect_equal(as.character(foo$mapping$x), 'group')
            expect_equal(as.character(foo$mapping$fill), 'group')
            expect_equal(foo$scales$scales[[1]]$scale_name, 'grey')
          })

#
# simpleBarGraph tests
#

test_that('simpleBarGraph fails with wrong dataVector',
          {
            expect_error(simpleBarGraph('foo'))
            expect_error(simpleBarGraph(TRUE))
            expect_error(simpleBarGraph(matrix(1:10, ncol=2)))
            wrongDataList <- list(g1=c('A', 'B', 'C'),
                                  g2=factor(c('A', 'A')),
                                  g3=c(4, 8, 15, 16, 23, 42))
            expect_error(simpleBarGraph(wrongDataList))
          })

test_that('simpleBarGraph works correctly on a naive example',
          {
            mockDataVector <- c(4, 8, 15, 16, 23, 42)
            names(mockDataVector) <- paste0('g', 1:6)
            foo <- simpleBarGraph(mockDataVector)
            expect_true(is(foo, 'gg'))
            expect_true(is(foo, 'ggplot'))
            expect_equal(length(names(foo)), 9)
            expect_equal(names(foo), c('data', 'layers', 'scales', 'mapping', 'theme', 
                                       'coordinates', 'facet', 'plot_env', 'labels'))
            expect_equal(as.character(foo$mapping$x), 'group')
            expect_equal(as.character(foo$mapping$fill), 'group')
            expect_equal(foo$scales$scales[[1]]$scale_name, 'grey')
          })

