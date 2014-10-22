context('General testing framework tests')

#
# rangeTestAgainstList tests
#
test_that('rangeTestAgainstList fails on incorrect or missing arguments',
          {
            expect_error(rangeTestAgainstList())
            expect_error(rangeTestAgainstList('one', 2))
            expect_error(rangeTestAgainstList(1, 'two', list(three=3)))
            expect_error(rangeTestAgainstList(NA, NA, NA))
            expect_error(rangeTestAgainstList(matrix(0, nrow=2, ncol=2), -2, -3))
          })

test_that('rangeTestAgainstList fails on empty ranges list',
          {
            mockRanges <- GRanges('foo', IRanges(1:3, 5:7))
            emptyRangesList <- GRangesList()
            expect_error(rangeTestAgainstList(mockRanges, mockRanges, emptyRangesList))
          })

test_that('rangesTestAgainstList fails if range list does not have names',
          {
            query <- GRanges('foo', IRanges(c(1,4,7), c(2, 5, 8)))
            background <- GRanges('foo', IRanges(c(6, 9), c(6, 10)))
            m1 <- GRanges('foo', IRanges(c(1, 4), c(2, 5)))
            m2 <- GRanges('foo', IRanges(6, 10))
            rglist <- GRangesList(m1, m2)
            expect_error(rangeTestAgainstList(query, background, rglist))
          })

test_that('rangesTestAgainstList works correctly on example',
          {
            query <- GRanges('foo', IRanges(c(1,4,7), c(2, 5, 8)))
            background <- GRanges('foo', IRanges(c(6, 9), c(6, 10)))
            m1 <- GRanges('foo', IRanges(c(1, 4), c(2, 5)))
            m2 <- GRanges('foo', IRanges(6, 10))
            rglist <- GRangesList(m1, m2)
            names(rglist) <- c('one', 'two')
            foo <- rangeTestAgainstList(query, background, rglist)
            expect_equal(unique(rowSums(foo[, 1:4])), 5)
            expect_equal(length(unique(rowSums(foo[, 1:4]))), 1)
#            expect_equal(
          })

#
# tids450kTestAgainstList tests
#

