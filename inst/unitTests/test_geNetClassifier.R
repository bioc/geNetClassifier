## geNetClassifier()

test_geNetClassifier <- function() 
{
	checkException(geNetClassifier(matrix(sample(50000,5*2),5, 2), c(rep("one", 2), rep("two", 2))), "The number of labels does not match the number of samples.")
	checkException(geNetClassifier(matrix(sample(50000,5*3),5, 3), c(rep("one", 2), rep("two", 1))), "There should be the same number of samples in each class.")
}
