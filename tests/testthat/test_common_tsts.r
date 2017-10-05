context('General testing framework tests')

library(epilab)

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

#
# categoricalTest tests
#

test_that('categoricalTest fails on incorrect or missing arguments',
          {
            expect_error(categoricalTest())
            expect_error(categoricalTest('one'))
            expect_error(categoricalTest(NA, NA))
            expect_error(categoricalTest(NULL, NULL))
            expect_error(categoricalTest(matrix(0, nrow=2, ncol=2), 1:3))
            expect_error(categoricalTest(letters[1:5], matrix(1, nrow=3, ncol=4)))
            expect_error(categoricalTest(letters[1:5], NA))
          })

test_that('categoricalTest (factor, numeric) fails when indices are wrong',
          {
            mockFactor <- factor(letters[1:10])
            mockIndices <- c(1, 2, 3, 4, 11)
            expect_error(categoricalTest(mockFactor, mockIndices))
          })

test_that('categoricalTest (factor, logical) fails when logical indices are wrong',
          {
            mockFactor <- factor(letters[1:10])
            expect_error(categoricalTest(mockFactor, rep(TRUE, 3)))
            expect_error(categoricalTest(mockFactor, rep(TRUE, 12)))
            mockLogical <- rep(TRUE, 10)
            mockLogical[7] <- NA
            expect_error(categoricalTest(mockFactor, mockLogical))
          })

test_that('categoricalTest (factor, character) fails when character indices are wrong',
          {
            mockFactor <- factor(letters[1:10])
            mockIndices <- c('n1', 'n2', 'n3')
            expect_error(categoricalTest(mockFactor, mockIndices))
            names(mockFactor) <- paste0('n', 1:10)
            mockIndices <- c('n1', 'n2', 'n12')
            expect_error(categoricalTest(mockFactor, mockIndices))
          })

test_that('categoricalTest (character, numeric) fails when indices are wrong',
          {
            mockText <- letters[1:10]
            mockIndices <- c(1, 2, 3, 4, 11)
            expect_error(categoricalTest(mockText, mockIndices))
          })

test_that('categoricalTest (character, logical) fails when logical indices are wrong',
          {
            mockText <- letters[1:10]
            expect_error(categoricalTest(mockText, rep(TRUE, 3)))
            expect_error(categoricalTest(mockText, rep(TRUE, 12)))
            mockLogical <- rep(TRUE, 10)
            mockLogical[7] <- NA
            expect_error(categoricalTest(mockText, mockLogical))
          })

test_that('categoricalTest (character, character) fails when character indices are wrong',
          {
            mockText <- letters[1:10]
            mockIndices <- c('n1', 'n2', 'n3')
            expect_error(categoricalTest(mockText, mockIndices))
            names(mockText) <- paste0('n', 1:10)
            mockIndices <- c('n1', 'n2', 'n12')
            expect_error(categoricalTest(mockText, mockIndices))
          })

test_that('categoricalTest fails when id is not a single string',
          {
            mockFactor <- factor(c('A', 'A', 'B', 'B', 'C', 'C', 'C', 'C'))
            mockIndices <- 2:6
            expect_error(categoricalTest(mockFactor, mockIndices, 42))
            expect_error(categoricalTest(mockFactor, mockIndices, list()))
            expect_error(categoricalTest(mockFactor, mockIndices, matrix(0, nrow=2, ncol=2)))
            expect_error(categoricalTest(mockFactor, mockIndices, c('foo', 'bar')))
          })

test_that('categoricalTest works correctly on (factor, numeric) signature',
          {
            mockFactor <- factor(c('A', 'A', 'B', 'B', 'C', 'C', 'C', 'C'))
            mockIndices <- 2:6
            foo <- categoricalTest(mockFactor, mockIndices, testId='foo')
            expect_equal(foo$Id, 'foo')
            expect_equal(foo$A_In, 1)
            expect_equal(foo$B_In, 2)
            expect_equal(foo$C_In, 2)
            expect_equal(foo$A_Out, 1)
            expect_equal(foo$B_Out, 0)
            expect_equal(foo$C_Out, 2)
            expect_equal(foo$PValue, 0.449329, tolerance=0.0001)
            expect_equivalent(foo$OR_A, 0.5)
            expect_equivalent(foo$OR_B, Inf)
            expect_equivalent(foo$OR_C, 1 / 3)
            expect_equal(foo$P_A, 0.2)
            expect_equal(foo$P_B, 0.4)
            expect_equal(foo$P_C, 0.4)
          })

test_that('categoricalTest works correctly on (factor, logical) signature',
          {
            mockFactor <- factor(c('A', 'A', 'B', 'B', 'C', 'C', 'C', 'C'))
            mockIndices <- c(FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE)
            foo <- categoricalTest(mockFactor, mockIndices, testId='foo')
            expect_equal(foo$Id, 'foo')
            expect_equal(foo$A_In, 1)
            expect_equal(foo$B_In, 2)
            expect_equal(foo$C_In, 2)
            expect_equal(foo$A_Out, 1)
            expect_equal(foo$B_Out, 0)
            expect_equal(foo$C_Out, 2)
            expect_equal(foo$PValue, 0.449329, tolerance=0.0001)
            expect_equivalent(foo$OR_A, 0.5)
            expect_equivalent(foo$OR_B, Inf)
            expect_equivalent(foo$OR_C, 1 / 3)
            expect_equal(foo$P_A, 0.2)
            expect_equal(foo$P_B, 0.4)
            expect_equal(foo$P_C, 0.4)
          })

test_that('categoricalTest works correctly on (factor, character) signature',
          {
            mockFactor <- factor(c('A', 'A', 'B', 'B', 'C', 'C', 'C', 'C'))
            names(mockFactor) <- paste0('n', 1:8)
            mockIndices <- paste0('n', 2:6)
            foo <- categoricalTest(mockFactor, mockIndices, testId='foo')
            expect_equal(foo$Id, 'foo')
            expect_equal(foo$A_In, 1)
            expect_equal(foo$B_In, 2)
            expect_equal(foo$C_In, 2)
            expect_equal(foo$A_Out, 1)
            expect_equal(foo$B_Out, 0)
            expect_equal(foo$C_Out, 2)
            expect_equal(foo$PValue, 0.449329, tolerance=0.0001)
            expect_equivalent(foo$OR_A, 0.5)
            expect_equivalent(foo$OR_B, Inf)
            expect_equivalent(foo$OR_C, 1 / 3)
            expect_equal(foo$P_A, 0.2)
            expect_equal(foo$P_B, 0.4)
            expect_equal(foo$P_C, 0.4)
          })

test_that('categoricalTest works correctly on (character, numeric) signature',
          {
            mockFactor <- c('A', 'A', 'B', 'B', 'C', 'C', 'C', 'C')
            mockIndices <- 2:6
            foo <- categoricalTest(mockFactor, mockIndices, testId='foo')
            expect_equal(foo$Id, 'foo')
            expect_equal(foo$A_In, 1)
            expect_equal(foo$B_In, 2)
            expect_equal(foo$C_In, 2)
            expect_equal(foo$A_Out, 1)
            expect_equal(foo$B_Out, 0)
            expect_equal(foo$C_Out, 2)
            expect_equal(foo$PValue, 0.449329, tolerance=0.0001)
            expect_equivalent(foo$OR_A, 0.5)
            expect_equivalent(foo$OR_B, Inf)
            expect_equivalent(foo$OR_C, 1 / 3)
            expect_equal(foo$P_A, 0.2)
            expect_equal(foo$P_B, 0.4)
            expect_equal(foo$P_C, 0.4)
          })

test_that('categoricalTest works correctly on (character, logical) signature',
          {
            mockFactor <- c('A', 'A', 'B', 'B', 'C', 'C', 'C', 'C')
            mockIndices <- c(FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE)
            foo <- categoricalTest(mockFactor, mockIndices, testId='foo')
            expect_equal(foo$Id, 'foo')
            expect_equal(foo$A_In, 1)
            expect_equal(foo$B_In, 2)
            expect_equal(foo$C_In, 2)
            expect_equal(foo$A_Out, 1)
            expect_equal(foo$B_Out, 0)
            expect_equal(foo$C_Out, 2)
            expect_equal(foo$PValue, 0.449329, tolerance=0.0001)
            expect_equivalent(foo$OR_A, 0.5)
            expect_equivalent(foo$OR_B, Inf)
            expect_equivalent(foo$OR_C, 1 / 3)
            expect_equal(foo$P_A, 0.2)
            expect_equal(foo$P_B, 0.4)
            expect_equal(foo$P_C, 0.4)
          })

test_that('categoricalTest works correctly on (character, character) signature',
          {
            mockFactor <- c('A', 'A', 'B', 'B', 'C', 'C', 'C', 'C')
            names(mockFactor) <- paste0('n', 1:8)
            mockIndices <- paste0('n', 2:6)
            foo <- categoricalTest(mockFactor, mockIndices, testId='foo')
            expect_equal(foo$Id, 'foo')
            expect_equal(foo$A_In, 1)
            expect_equal(foo$B_In, 2)
            expect_equal(foo$C_In, 2)
            expect_equal(foo$A_Out, 1)
            expect_equal(foo$B_Out, 0)
            expect_equal(foo$C_Out, 2)
            expect_equal(foo$PValue, 0.449329, tolerance=0.0001)
            expect_equivalent(foo$OR_A, 0.5)
            expect_equivalent(foo$OR_B, Inf)
            expect_equivalent(foo$OR_C, 1 / 3)
            expect_equal(foo$P_A, 0.2)
            expect_equal(foo$P_B, 0.4)
            expect_equal(foo$P_C, 0.4)
          })

#
# continuousTest tests
#

test_that('continuousTest fails on incorrect or missing arguments',
          {
            expect_error(continuousTest())
            expect_error(continuousTest('one'))
            expect_error(continuousTest(NA, NA))
            expect_error(continuousTest(NULL, NULL))
            expect_error(continuousTest(matrix(0, nrow=2, ncol=2), 1:3))
            expect_error(continuousTest(letters[1:5], matrix(1, nrow=3, ncol=4)))
            expect_error(continuousTest(letters[1:5], NA))
            expect_error(continuousTest(factor(1, 10), 1:5))
            expect_error(continuousTest(letters[1:7], 1:5))
          })

test_that('continuousTest (indexed by numeric) fails when indices are wrong',
          {
            mockVar <- 1:10
            mockIndices <- c(1, 2, 3, 4, 11)
            expect_error(continuousTest(mockVar, mockIndices))
          })

test_that('continuousTest (indexed by logical) fails when logical indices are wrong',
          {
            mockVar <- 1:10
            expect_error(continuousTest(mockVar, rep(TRUE, 3)))
            expect_error(continuousTest(mockVar, rep(TRUE, 12)))
            mockLogical <- rep(TRUE, 10)
            mockLogical[7] <- NA
            expect_error(continuousTest(mockVar, mockLogical))
          })

test_that('continuousTest (indexed by character) fails when character indices are wrong',
          {
            mockVar <- 1:10
            mockIndices <- c('n1', 'n2', 'n3')
            expect_error(continuousTest(mockVar, mockIndices))
            names(mockVar) <- paste0('n', 1:10)
            mockIndices <- c('n1', 'n2', 'n12')
            expect_error(continuousTest(mockVar, mockIndices))
          })

test_that('continuousTest fails when id is not a single string',
          {
            mockVar <- 1:10
            mockIndices <- 2:6
            expect_error(continuousTest(mockVar, mockIndices, 42))
            expect_error(continuousTest(mockVar, mockIndices, list()))
            expect_error(continuousTest(mockVar, mockIndices, matrix(0, nrow=2, ncol=2)))
            expect_error(continuousTest(mockVar, mockIndices, c('foo', 'bar')))
          })

test_that('continuousTest works correctly on numeric indexing',
          {
            mockVar <- 1:10
            mockIndices <- 2:6
            foo <- continuousTest(mockVar, mockIndices, testId='foo')
            expect_equal(foo$Id, 'foo')
            expect_equivalent(foo$Median_In, 4)
            expect_equivalent(foo$Median_Out, 8)
            expect_equal(foo$PValue, 0.1507937, tolerance=0.0001)
            expect_equal(foo$AC, 0.1653926, tolerance=0.0001)
            expect_equal(foo$D, -0.6)
          })

test_that('continuousTest works correctly on logical indexing',
          {
            mockVar <- 1:10
            mockIndices <- c(FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE)
            foo <- continuousTest(mockVar, mockIndices, testId='foo')
            expect_equal(foo$Id, 'foo')
            expect_equivalent(foo$Median_In, 4)
            expect_equivalent(foo$Median_Out, 8)
            expect_equal(foo$PValue, 0.1507937, tolerance=0.0001)
            expect_equal(foo$AC, 0.1653926, tolerance=0.0001)
            expect_equal(foo$D, -0.6)
          })

test_that('continuousTest works correctly on character indexing',
          {
            mockVar <- 1:10
            names(mockVar) <- paste0('n', 1:10)
            mockIndices <- paste0('n', 2:6)
            foo <- continuousTest(mockVar, mockIndices, testId='foo')
            expect_equal(foo$Id, 'foo')
            expect_equivalent(foo$Median_In, 4)
            expect_equivalent(foo$Median_Out, 8)
            expect_equal(foo$PValue, 0.1507937, tolerance=0.0001)
            expect_equal(foo$AC, 0.1653926, tolerance=0.0001)
            expect_equal(foo$D, -0.6)
          })


