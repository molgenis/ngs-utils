SumOfSquares <- function(forward.column, samplefractions)
{
  squares <- (forward.column - samplefractions) ^ 2
  return(squares)
}

BestControlSet <-function(chromosomal.fractions, samplefractions, n.of.samples)
{

fractions.squared <- apply(chromosomal.fractions, 1, SumOfSquares, samplefractions = samplefractions)

sum.of.squares <- apply(fractions.squared[control.chromosomes,], 2, sum)

sorted.scores <- order(sum.of.squares)

return (sorted.scores[1:n.of.samples])
}





