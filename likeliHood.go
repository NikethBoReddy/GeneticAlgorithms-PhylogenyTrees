package main

import (
  "fmt"
  "math"
  "os"
)

/*--------------------------------------------------------------------------------------------------------------------------
 * Given a model phylogenic tree, the function calculates the probability that this model can explain the data that is observed
 * for a given sequence length. Although the problem consists of an exponential number of configuration for which probabilities
 * need to be calculated, there exists a dynamic programming solution that can solve it in polynomial time
 * Currently the slowest step in the algorithm; Needs to be improved.
 *------------------------------------------------------------------------------------------------------------------------*/
func CalculateMaxLikelihoodScores(solutionModel *node, speciesMap map[string]speciesGenome, sequenceLength int) float64 {
  var logScore float64 = 0
  for i:=0; i<sequenceLength; i++ {
    var score float64 = 0
    scoreMaps := make(map[*node]([5]float64))
    var currScore, currScore1, currScore2 float64
    var netScore float64
    for j:=0; j<5; j++ {
      currScore = solutionModel.nucleotideFrequencies[j]
      for k:=0; k<5; k++ {
        currScore1 = currScore*TransitionProbability(solutionModel, true, j, k)*
                                  CalculateRecursiveLikelihood(solutionModel.leftChild, k, scoreMaps, speciesMap, i)
        for l:=0; l<5; l++ {
          currScore2 = currScore1*TransitionProbability(solutionModel, false, j, l)*
                                  CalculateRecursiveLikelihood(solutionModel.rightChild, l, scoreMaps, speciesMap, i)
          netScore += currScore2
        }
      }
      score += netScore
    }
    logScore +=   math.Log(score)
  }
  return logScore
}

/*-------------------------------------------------------------------------------------------------------------------------
 * For a particular species, calculates its probability of making a change in nucleotide sequence from one base to another
 *-----------------------------------------------------------------------------------------------------------------------*/
func TransitionProbability(currNode *node, direction bool, nct1, nct2 int) float64 {
  var probability float64
  var branchLength float64

  if direction {
    branchLength = currNode.leftChildDistance
  } else {
    branchLength = currNode.rightChildDistance
  }

  if nct1 == nct2 {
    probability = (1.0-branchLength)
  } else {
    probability = branchLength
  }
  probability *= currNode.conversionRatios[nct1][nct2]
  return probability
}

/*---------------------------------------------------------------------------------------------------------------------------
 * The meoized recursive solution for estimating the maxlikelihood score, similar to it's parent function with the exception of
 * holding the solution for the base case.
 *--------------------------------------------------------------------------------------------------------------------------*/
func CalculateRecursiveLikelihood(currNode *node, nct int, scoreMaps map[*node]([5]float64), speciesMap map[string]speciesGenome, index int) float64 {
  if currNode == nil {
    fmt.Println("Invalid calculation of likelihood; null node detected!! ")
    os.Exit(1)
  }
  val, exists := scoreMaps[currNode]
  if exists {
    return val[nct]
  }
  if currNode.leftChild == nil && currNode.rightChild == nil {
    species, exists := speciesMap[currNode.name]
    if !exists {
      fmt.Println("Invalid Tree or Map present")
      os.Exit(1)
    }
    currNct := species.nucleotideSequence[index]
    var nctScores [5]float64
    if currNct == 'A' || currNct == 'a' {
      nctScores[0] = currNode.nucleotideFrequencies[0]
    } else if currNct == 'C' || currNct == 'c' {
      nctScores[1] = currNode.nucleotideFrequencies[1]
    } else if currNct == 'G' || currNct == 'g' {
      nctScores[2] = currNode.nucleotideFrequencies[2]
    } else if currNct == 'T' || currNct == 't' {
      nctScores[3] = currNode.nucleotideFrequencies[3]
    } else if currNct == '-' || currNct == '.' {
      nctScores[4] = currNode.nucleotideFrequencies[4]
    } else {
      fmt.Println("Invalid nucleotide base detected in input !!")
      os.Exit(1)
    }

    scoreMaps[currNode] = nctScores
    return nctScores[nct]
  }

  var nctScores [5]float64
  var currScore, currScore1, currScore2 float64
  var netScore float64
  for i:=0; i<5; i++ {
    currScore = currNode.nucleotideFrequencies[i]
    for j:=0; j<5; j++ {
      currScore1 = currScore*TransitionProbability(currNode, true, i, j)*
                                CalculateRecursiveLikelihood(currNode.leftChild, j, scoreMaps, speciesMap, index)
      for k:=0; k<5; k++ {
        currScore2 = currScore1*TransitionProbability(currNode, false, i, k)*
                                CalculateRecursiveLikelihood(currNode.rightChild, k, scoreMaps, speciesMap, index)
        netScore += currScore2
      }
    }
    nctScores[i] = netScore
  }

  scoreMaps[currNode] = nctScores
  return nctScores[nct]

}
