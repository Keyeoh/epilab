context('Fileparse tests')

#
# whichOne() tests
#
test_that('whichOne fails with wrong inputs',
          {
            foo <- c(4, 8, 15, 16, 23, 42)
            expect_error(whichOne(foo, 0))
            expect_error(whichOne(foo, -1))
          })

test_that('whichOne works correctly on naive example',
          {
            foo <- c(4, 8, 15, 16, 23, 42)
            expect_equal(whichOne(foo, 1), 4)
            expect_equal(whichOne(foo, 6), 42)
          })

#
# first() tests
#
test_that('first works correctly on naive example',
          {
            foo <- c(4, 8, 15, 16, 23, 42)
            expect_equal(first(foo), 4)
          })

#
# second() tests
#
test_that('second works correctly on naive example',
          {
            foo <- c(4, 8, 15, 16, 23, 42)
            expect_equal(second(foo), 8)
          })

#
# matchLength() tests
#
test_that('matchLength fails when attribute is not present',
          {
            foo <- c(4, 8, 15, 16, 23, 42)
            expect_error(matchLength(foo))
          })

test_that('matchLength works correctly on naive example',
          {
            foo <- c(4, 8, 15, 16, 23, 42)
            attr(foo, 'match.length') <- 666
            expect_equal(matchLength(foo), 666)
          })

#
# secondLength() tests
#
test_that('secondLength fails when attribute is not present',
          {
            foo <- c(4, 8, 15, 16, 23, 42)
            expect_error(secondLength(foo))
          })

test_that('secondLength works correctly on naive example',
          {
            foo <- c(4, 8, 15, 16, 23, 42)
            attr(foo, 'match.length') <- c(666, 999)
            expect_equal(secondLength(foo), 999)
          })

#
# getValuesForKey() tests
#
test_that('getValuesForKey fails with wrong inputs',
          {
            foo <- c('line1', 'line2', 'line3')
            expect_error(getValuesForKey(32, foo))
            expect_error(getValuesForKey(matrix(rnorm(9), ncol=3), foo))
            expect_error(getValuesForKey(c('key1', 'key2'), foo))
            expect_error(getValuesForKey('key', 1:7))
            expect_error(getValuesForKey('key', c(TRUE, FALSE, TRUE)))
            expect_error(getValuesForKey('key', matrix(rnorm(9), ncol=3)))
          })

test_that('getValuesForKey works correctly on naive example',
          {
            foo <- c('line1\tkey1=value1; key2=value2;', 'line2\tkey3=value3;')
            expect_equal(getValuesForKey('key1', foo), c('value1', NA))
            expect_equal(getValuesForKey('key2', foo), c('value2', NA))
            expect_equal(getValuesForKey('key3', foo), c(NA, 'value3'))
          })

#
# createFilesData() tests
#
test_that('createFilesData fails with wrong inputs',
          {
            expect_error(createFilesData(1:7))
            expect_error(createFilesData(c(TRUE, FALSE, TRUE)))
            expect_error(createFilesData(matrix(rnorm(9), ncol=3)))
          })

test_that('createFilesData works correctly on naive example',
          {
            foo <- c('line1.broadPeak.gz\tkey1=value1; key2=value2;', 
                     'line2.broadPeak.gz\tkey3=value3;')
            bar <- data.frame(stem=c('line1', 'line2'), stringsAsFactors=FALSE)
            expect_equal(createFilesData(foo), bar)
          })

#
# getAllKeys() tests
#
test_that('getAllKeys fails with wrong inputs',
          {
            expect_error(getAllKeys(1:7))
            expect_error(getAllKeys(c(TRUE, FALSE, TRUE)))
            expect_error(getAllKeys(matrix(rnorm(9), ncol=3)))
          })

test_that('getAllKeys works correctly on naive example',
          {
            foo <- c('line1.broadPeak.gz\tkey1=value1; key2=value2;', 
                     'line2.broadPeak.gz\tkey3=value3;')
            expect_equal(getAllKeys(foo), c('key1', 'key2', 'key3'))
          })

#
# addKey() tests
#
test_that('addKey fails with wrong inputs',
          {
            mockDataFrame <- data.frame(stem=c('line1', 'line2'))
            mockLines <- c('line1.broadPeak.gz\tkey1=value1; key2=value2;', 
                           'line2.broadPeak.gz\tkey3=value3;')
            expect_error(addKey(1:7, 'key1', mockLines))
            expect_error(addKey(matrix(1:9, ncol=3), 'key1', mockLines))
            expect_error(addKey(c(TRUE, FALSE), 'key1', mockLines))
            expect_error(addKey(factor(1:4), 'key1', mockLines))
            expect_error(addKey(mockDataFrame, 4, mockLines))
            expect_error(addKey(mockDataFrame, TRUE, mockLines))
            expect_error(addKey(mockDataFrame, matrix(1:9, ncol=3), mockLines))
            expect_error(addKey(mockDataFrame, factor(4:7), mockLines))
            expect_error(addKey(mockDataFrame, 'key1', 1:7))
            expect_error(addKey(mockDataFrame, 'key1', c(TRUE, FALSE, TRUE)))
            expect_error(addKey(mockDataFrame, 'key1', matrix(1:9, ncol=3)))          
          })

test_that('addKey works correctly on naive example',
          {
            mockDataFrame <- data.frame(stem=c('line1', 'line2'))
            mockLines <- c('line1.broadPeak.gz\tkey1=value1; key2=value2;', 
                           'line2.broadPeak.gz\tkey3=value3;')
            foo1 <- addKey(mockDataFrame, 'key1', mockLines)
            expect_equal(dim(foo1), c(2, 2))
            expect_equal(foo1$key1, c('value1', NA))
            foo2 <- addKey(mockDataFrame, 'key2', mockLines)
            expect_equal(dim(foo2), c(2, 2))
            expect_equal(foo2$key2, c('value2', NA))
            foo3 <- addKey(mockDataFrame, 'key3', mockLines)
            expect_equal(dim(foo3), c(2, 2))
            expect_equal(foo3$key3, c(NA, 'value3'))
          })

#
# addKeyFactory() tests
#
test_that('addKeyFactory fails with wrong inputs',
          {
            expect_error(addKeyFactory(1:7))
            expect_error(addKeyFactory(c(TRUE, FALSE)))
            expect_error(addKeyFactory(matrix(1:9, ncol=3)))
            expect_error(addKeyFactory(factor(1:4)))
          })

test_that('addKeyFactory works correctly on naive example',
          {
            mockDataFrame <- data.frame(stem=c('line1', 'line2'))
            mockLines <- c('line1.broadPeak.gz\tkey1=value1; key2=value2;', 
                           'line2.broadPeak.gz\tkey3=value3;')
            bar <- addKeyFactory(mockLines)
            foo1 <- bar(mockDataFrame, 'key1')
            expect_equal(dim(foo1), c(2, 2))
            expect_equal(foo1$key1, c('value1', NA))
            foo2 <- bar(mockDataFrame, 'key2')
            expect_equal(dim(foo2), c(2, 2))
            expect_equal(foo2$key2, c('value2', NA))
            foo3 <- bar(mockDataFrame, 'key3')
            expect_equal(dim(foo3), c(2, 2))
            expect_equal(foo3$key3, c(NA, 'value3'))
          })

#
# addAllKeys() tests
#
test_that('addAllKeys fails with wrong inputs',
          {
            mockDataFrame <- data.frame(stem=c('line1', 'line2'))
            mockLines <- c('line1.broadPeak.gz\tkey1=value1; key2=value2;', 
                           'line2.broadPeak.gz\tkey3=value3;')
            expect_error(addAllKeys(1:7, 'key1', mockLines))
            expect_error(addAllKeys(matrix(1:9, ncol=3), 'key1', mockLines))
            expect_error(addAllKeys(c(TRUE, FALSE), 'key1', mockLines))
            expect_error(addAllKeys(factor(1:4), 'key1', mockLines))
            expect_error(addAllKeys(mockDataFrame, 4, mockLines))
            expect_error(addAllKeys(mockDataFrame, TRUE, mockLines))
            expect_error(addAllKeys(mockDataFrame, matrix(1:9, ncol=3), mockLines))
            expect_error(addAllKeys(mockDataFrame, factor(4:7), mockLines))
            expect_error(addAllKeys(mockDataFrame, 'key1', 1:7))
            expect_error(addAllKeys(mockDataFrame, 'key1', c(TRUE, FALSE, TRUE)))
            expect_error(addAllKeys(mockDataFrame, 'key1', matrix(1:9, ncol=3)))          
          })

test_that('addAllKeys works correctly on naive example',
          {
            mockDataFrame <- data.frame(stem=c('line1', 'line2'), stringsAsFactors=FALSE)
            mockLines <- c('line1.broadPeak.gz\tkey1=value1; key2=value2;', 
                           'line2.broadPeak.gz\tkey3=value3;')
            foo <- addAllKeys(mockDataFrame, c('key1', 'key2', 'key3'), mockLines)
            expect_equal(dim(foo), c(2, 4))
            expect_equal(foo$stem, c('line1', 'line2'))
            expect_equal(foo$key1, c('value1', NA))
            expect_equal(foo$key2, c('value2', NA))
            expect_equal(foo$key3, c(NA, 'value3'))
          })

#
# readBroadHistoneFile() tests
#
test_that('readBroadHistoneFile fails with wrong inputs',
          {
            expect_error(readBroadHistoneFile(1:7))
            expect_error(readBroadHistoneFile(c(TRUE, FALSE)))
            expect_error(readBroadHistoneFile(matrix(1:9, ncol=3)))
            expect_error(readBroadHistoneFile(factor(1:4)))
          })

test_that('readBroadHistoneFile works correctly on naive example',
          {
            mockLines <- c('line1.broadPeak.gz\tkey1=value1; key2=value2;', 
                           'line2.broadPeak.gz\tkey3=value3;')
            foo <- readBroadHistoneFile(mockLines)
            expect_equal(dim(foo), c(2, 4))
            expect_equal(foo$stem, c('line1', 'line2'))
            expect_equal(foo$key1, c('value1', NA))
            expect_equal(foo$key2, c('value2', NA))
            expect_equal(foo$key3, c(NA, 'value3'))
          })

