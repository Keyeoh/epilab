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
            expect_equal(foo['one', 'Mark_Selected'], 2)
            expect_equal(foo['one', 'Mark_Background'], 0)
            expect_equal(foo['one', 'NoMark_Selected'], 1)
            expect_equal(foo['one', 'NoMark_Background'], 2)
            expect_equal(foo['one', 'P_Background'], 0)
            expect_equal(foo['one', 'OR'], Inf)
            expect_equal(foo['two', 'Mark_Selected'], 1)
            expect_equal(foo['two', 'Mark_Background'], 2)
            expect_equal(foo['two', 'NoMark_Selected'], 2)
            expect_equal(foo['two', 'NoMark_Background'], 0)
            expect_equal(foo['two', 'P_Background'], 1)
            expect_equal(foo['two', 'OR'], 0)
          })

#
# tids450kTestAgainstList tests
#
test_that('tids450kTestAgainstList fails on incorrect or missing arguments',
          {
            expect_error(tids450kTestAgainstList())
            expect_error(tids450kTestAgainstList('one'))
            expect_error(tids450kTestAgainstList(1, 'two'))
            expect_error(tids450kTestAgainstList(NA, NA))
            expect_error(tids450kTestAgainstList(matrix(0, nrow=2, ncol=2), -3))
          })

test_that('tids450kTestAgainstList fails on empty ranges list',
          {
            mockTids <- letters[1:10]
            emptyRangesList <- GRangesList()
            expect_error(tids450kTestAgainstList(mockTids, emptyRangesList))
          })

test_that('tids450kTestAgainstList fails if range list does not have names',
          {
            mockTids <- letters[1:10]
            m1 <- GRanges('foo', IRanges(c(1, 4), c(2, 5)))
            m2 <- GRanges('foo', IRanges(6, 10))
            rglist <- GRangesList(m1, m2)
            expect_error(tids450kTestAgainstList(mockTids, rglist))
          })

test_that('tids450kTestAgainstList works correctly on example',
          {
            mockTids <- c("cg04913815", "cg01686861", "cg05558259", "cg26978960", "cg03792876",
                          "cg09699726", "cg07549526", "cg02851049", "cg11876012", "cg14820573")
            m1 <- GRanges('chr16', IRanges(60000, 61000))
            m2 <- GRanges('chr16', IRanges(70000, 92000))
            rglist <- GRangesList(m1, m2)
            names(rglist) <- c('one', 'two')
            foo <- tids450kTestAgainstList(mockTids, rglist)
            expect_equal(unique(rowSums(foo[, 1:4])), 485577)
            expect_equal(length(unique(rowSums(foo[, 1:4]))), 1)
            expect_equal(foo['one', 'Mark_Selected'], 2)
            expect_equal(foo['one', 'Mark_Background'], 0)
            expect_equal(foo['one', 'NoMark_Selected'], 8)
            expect_equal(foo['one', 'NoMark_Background'], 485567)
            expect_equal(foo['one', 'P_Background'], 0)
            expect_equal(foo['one', 'OR'], Inf)
            expect_equal(foo['two', 'Mark_Selected'], 3)
            expect_equal(foo['two', 'Mark_Background'], 0)
            expect_equal(foo['two', 'NoMark_Selected'], 7)
            expect_equal(foo['two', 'NoMark_Background'], 485567)
            expect_equal(foo['two', 'P_Background'], 0)
            expect_equal(foo['two', 'OR'], Inf)
          })

