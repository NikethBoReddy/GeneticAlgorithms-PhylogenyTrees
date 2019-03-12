package main

import(
  "fmt"
  "time"
  "strconv"
  "math/rand"
)

/*-------------------------------------------------------------------------------------------------------
 * The core implementation of the GA algorithm, performs score calculations, population selection along with
 * the necessary mutations and crossovers. Also prints the generated outputs and other info with a
 * predetermined saampling rate.
 *------------------------------------------------------------------------------------------------------*/
func RunGASimulations(startingPopulation []*node, speciesMap map[string]speciesGenome,
                      numGenerations, minStableGenerations int) *node {

  numSpecies := len(speciesMap)
  numSolutions := len(startingPopulation)
  fmt.Println("\nStarting the genetic algorithm with a population of " + strconv.Itoa(numSolutions) + " solutions over " + strconv.Itoa(numGenerations) + " generations")

  fittestSurvivalReproductionRate := LoadFloatConfig("ga.algo.params.selection.proliferation.fraction")

  // Limiting the max sequence lengths for likelihood estimation.
  var sequenceLength int
  for _, val := range speciesMap {
    sequenceLength = len(val.nucleotideSequence)
    break
  }
  sequenceLimit := LoadIntConfig("ga.algo.params.sequencedata.length.max")
  if sequenceLength > sequenceLimit {
    sequenceLength = sequenceLimit
  }

  maxLikelihoodScores := make([]float64, numGenerations)
  var bestSolution *node

  samplingRate := LoadIntConfig("ga.output.sampling.interval")
  rand.Seed(time.Now().UTC().UnixNano())

  for i:=0; i<numGenerations; i++ {

    printStatistics := false
    if i % samplingRate == 0 {
      printStatistics = true
      fmt.Println("Successfully completed " + strconv.Itoa(i) + " iterations.")
    }

    likelihoodScores := make([]float64, numSolutions)
    for j:=0; j<numSolutions; j++ {
      likelihoodScores[j] = CalculateMaxLikelihoodScores(startingPopulation[j], speciesMap, sequenceLength)
    }

    sortedLikelihoods, sortedPopulation := SortDescending(likelihoodScores, startingPopulation)
    maxLikelihoodScores[i] = sortedLikelihoods[0]
    bestSolution = sortedPopulation[0]
    if printStatistics {
      fmt.Print("The Likelihood score has been optimized to "); fmt.Println(sortedLikelihoods[0])
      PrintNewickFormatTree(NewickFormatTreeRepresentation(bestSolution))
    }

    futurePopulation := GenerateFuturePopulation(fittestSurvivalReproductionRate, sortedPopulation)
    MutateFuturePopulation(startingPopulation, futurePopulation, numSpecies)

    startingPopulation = futurePopulation

    if stabilityAcheived(maxLikelihoodScores, i, minStableGenerations) {
      fmt.Println("No significant variation in likelihood scores over the last " + strconv.Itoa(minStableGenerations) + " generations")
      fmt.Println("Terminating the GA Algorithm\n")
      return bestSolution
    }
  }

  return bestSolution
}

/*----------------------------------------------------------------------------------------------------
 * Since the dataset doesn't contain very large number of species, a very naive sorting algorithm is
 * implemented for sorting the likelihood scores and the population accordingly
 *---------------------------------------------------------------------------------------------------*/
func SortDescending(likelihoodScores []float64, startingPopulation []*node) ([]float64, []*node) {
  populationCount := len(likelihoodScores)
  sortedScores := make([]float64, populationCount)
  sortedPopulation := make([]*node, populationCount)
  for i:=0; i<populationCount; i++ {
    maxVal := likelihoodScores[0]
    maxIndex := 0
    for j:=0; j<populationCount-i; j++ {
      if likelihoodScores[j] > maxVal {
        maxVal = likelihoodScores[j]
        maxIndex = j
      }
    }
    sortedScores[i] = maxVal
    sortedPopulation[i] = startingPopulation[maxIndex]

    likelihoodScores = append(likelihoodScores[:maxIndex], likelihoodScores[maxIndex+1:]...)
    startingPopulation = append(startingPopulation[:maxIndex], startingPopulation[maxIndex+1:]...)
  }
  return sortedScores, sortedPopulation
}

/*-------------------------------------------------------------------------------------------------------
 * Given a sorted collection of trees according to their scores, we generate the next geenration of solutions
 * by selecting trees from the parent population with the appropriate probabilities
 *------------------------------------------------------------------------------------------------------*/
func GenerateFuturePopulation(fittestSurvivalReproductionRate float64, sortedPopulation []*node) []*node {
  futurePopulationCounter := 0
  numSolutions := len(sortedPopulation)
  futurePopulation := make([]*node, numSolutions)

  for j:=0; j<int(fittestSurvivalReproductionRate*float64(numSolutions)) && j<numSolutions; j++ {
    futurePopulation[j] = GenerateTreeCopy(sortedPopulation[0])
    futurePopulationCounter++
  }
  for j:=1; j<numSolutions && futurePopulationCounter<numSolutions; j++ {
    survivalRate := float64(2)/float64((j+1)*(j+2))
    if rand.Float64() < survivalRate {
      if futurePopulationCounter < numSolutions {
        futurePopulation[futurePopulationCounter] = GenerateTreeCopy(sortedPopulation[j])
        futurePopulationCounter++
      }
    }
    if j == numSolutions-1 {
      j = 0
    }
  }
  return futurePopulation
}

/*-------------------------------------------------------------------------------------------------------
 * It is necessary to generate a completely new tree for the offspring population as we do not want multiple
 * trees sharing the same pointers. Also helps in making the mutation and rearrangement function a lot easier
 * and managable. A simple recursive function for implementing the samplingRate
 *-------------------------------------------------------------------------------------------------------*/
func GenerateTreeCopy(root *node) *node {
  if root == nil {
    return nil
  }
  var newRoot node
  newRoot.name = root.name
  newRoot.leftChildDistance = root.leftChildDistance
  newRoot.rightChildDistance = root.rightChildDistance
  for i:=0; i<5; i++ {
    newRoot.nucleotideFrequencies[i] = root.nucleotideFrequencies[i]
    for j:=0; j<5; j++ {
      newRoot.conversionRatios[i][j] = root.conversionRatios[i][j]
    }
  }
  newRoot.leftChild = GenerateTreeCopy(root.leftChild)
  newRoot.rightChild = GenerateTreeCopy(root.rightChild)
  if newRoot.leftChild != nil {
    newRoot.leftChild.parent = &newRoot
  }
  if newRoot.rightChild != nil {
    newRoot.rightChild.parent = &newRoot
  }
  return &newRoot
}

/*----------------------------------------------------------------------------------------------
 * Function for checking the activity of the GA. Will terminate the program if no significant
 * variation is observed over a given set of iterations.
 *---------------------------------------------------------------------------------------------*/
func stabilityAcheived(maxLikelihoodScores []float64, currIndex, minStableGenerations int) bool {
  if currIndex+1 < minStableGenerations {
    return false
  } else {
    currValue := maxLikelihoodScores[currIndex]
    for i:=1; i<minStableGenerations; i++ {
        if maxLikelihoodScores[currIndex - i] != currValue {
          return false
        }
    }
  }
  return true
}
