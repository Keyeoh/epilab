context('BatchPlot tests')

#
# Mock objects for testing
#
mockPdata <- data.frame(
                        Q1=1:10,
                        Q2=rnorm(10),
                        N1=gl(2, 5),
                        N2=as.factor(rep(c(1, 2), 5))
                        )
mockSV <- matrix(c(1:10,
                   1:10,
                   10:1,
                   10:1
                   ),
                 nrow=10
                 )
mockNVars <- c('N1', 'N2')
mockQVars <- c('Q1', 'Q2')

mockSkeleton <- data.frame(
                           varname=gl(4, 4, labels=c('N1', 'N2', 'Q1', 'Q2')),
                           svname=gl(4,1,16,labels=paste0('SV', 1:4)),
                           type=gl(2, 8, labels=c('n', 'q'))
                           )

#
# Batch plot tests
#
test_that('batchPlot fails on wrong parameters', 
          {
            expect_error(batchPLot())
            expect_error(batchPlot(1, 2, 3, 4))
            expect_error(batchPlot(1, 2, 3, 4, 5))
            expect_error(batchPlot(1, mockSV, mockNVars, mockQVars))
            expect_error(batchPlot(mockPdata, mockSV, 2, mockQVars))
            expect_error(batchPlot(mockPdata, mockSV, mockNVars, 3))
            expect_error(batchPlot(mockPdata, 4, mockNVars, mockQVars))
          })

test_that('batchPlot generates a well-formed structure', 
          {
            foo <- batchPlot(mockPdata, mockSV, mockNVars, mockQVars)
            expect_true(all(foo[, 1:3] == mockSkeleton))
          })

test_that('batchPlot generates correct scores', 
          {
            foo <- batchPlot(mockPdata, mockSV, mockNVars, mockQVars)
            expect_equal(foo$score[1:4], rep(1, 4))
            expect_equal(foo$score[5:8], rep(0, 4))
            expect_equal(foo$score[9:10], c(1, 1))
            expect_equal(foo$score[11:12], c(-1, -1))
            expect_equal(foo$score[13:14], -foo$score[15:16])
          })

