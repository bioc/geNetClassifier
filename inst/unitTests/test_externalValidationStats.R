
##
## externalValidation.stats
## function(confussionMatrix, numDecimals=2)
##
## Calculates stats from the confussion matrix. 
# Input: 	Receives a confussion matrix (actual(rows) x prediction(cols)) --> Columns and rows in same order. "NotAssigned" last column
# Output: List (Sensitivity and specificity of the predictor for each class, global accuracy and call rate (rate of assigned samples))

test_externalValidationStats <- function()
{
	checkTrue(is.list(externalValidation.stats(table(c("A","B","C", "A", "B"), c("A","B","C", "B", "B")))))
}
