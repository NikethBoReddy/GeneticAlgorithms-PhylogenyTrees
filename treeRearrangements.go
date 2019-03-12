package main

import (
  "fmt"
  "os"
  "time"
  "math/rand"
  "gonum.org/v1/gonum/stat/distuv"
)

/*-----------------------------------------------------------------------------------------------------
 * The abstracted function represting the various types of mutations that are involved in the GA algo
 *---------------------------------------------------------------------------------------------------*/
func MutateFuturePopulation(parentPopulation, population []*node, numSpecies int) {
  numSolutions := len(population)

  branchMutationRate := LoadFloatConfig("ga.algo.params.mutation.branchlength")
  nucleotideMutationRate := LoadFloatConfig("ga.algo.params.mutation.nucleotide")
  topologyMutationRate := LoadFloatConfig("ga.algo.params.mutation.topology")
  recombinationProbability := LoadFloatConfig("ga.algo.params.crossover")

  rand.Seed(time.Now().UTC().UnixNano())
  for i:=1; i<numSolutions; i++ {
    MutateBranches(population[i], branchMutationRate)
    MutateNucleotideFrequencies(population[i], nucleotideMutationRate)
    population[i] = MutateTopology(population[i], topologyMutationRate, numSpecies)
    population[i] = PerformCrossOver(population[i], parentPopulation, recombinationProbability, numSpecies)
  }

}

/*-----------------------------------------------------------------------------------------------------
 * Branch Lengths are mutated with a given probability by multiplying them with a sample drawn from a
 * gamma distribution with alpha = 500 and mean = 1
 *----------------------------------------------------------------------------------------------------*/
func MutateBranches(solution *node, rate float64)  {
  if solution == nil || solution.rightChild == nil || solution.leftChild == nil {
    return
  }
  leftChance := rand.Float64()
  if leftChance < rate {
    solution.leftChildDistance *= GammaDistrubution(500, 500)
    if solution.leftChildDistance > 1 {
      solution.leftChildDistance = 1
    } else if solution.leftChildDistance <= 0.001 {
      solution.leftChildDistance = 0.001
    }
  }
  rightChance := rand.Float64()
  if rightChance < rate {
    solution.rightChildDistance *= GammaDistrubution(500, 500)
    // Care has to be taken such that branch lengths do not overflow or become too small
    if solution.rightChildDistance > 1 {
      solution.rightChildDistance = 1
    } else if solution.rightChildDistance < 0.001 {
      solution.rightChildDistance = 0.001
    }
  }

  MutateBranches(solution.leftChild, rate)
  MutateBranches(solution.rightChild, rate)

}

/*----------------------------------------------------------------------------------------------------
 * The nucleotide conversion ratios and base frequencies are altered in a similar fashion by multiplying
 * them with samples from the same gamma distribution curves and adjusted such that the net probability
 * for each transition type equals one
 *----------------------------------------------------------------------------------------------------*/
func MutateNucleotideFrequencies(currentNode *node, nucleotideMutationRate float64) {
  if currentNode == nil {
    return
  }
  for i:=0; i<5; i++ {
    chance := rand.Float64()
    if chance < nucleotideMutationRate {
      currentNode.nucleotideFrequencies[i] *= GammaDistrubution(500, 500)
    }
    for j:=0; j<5; j++ {
      chance := rand.Float64()
      if chance < nucleotideMutationRate {
          currentNode.conversionRatios[i][j] *= GammaDistrubution(500,500)
      }
    }
  }
  var frequencySum float64
  var rowSums [5]float64
  for i:=0; i<5; i++ {
    frequencySum += currentNode.nucleotideFrequencies[i]
    for j:=0; j<5; j++ {
      rowSums[i] += currentNode.conversionRatios[i][j]
    }
  }
  for i:=0; i<5; i++ {
    currentNode.nucleotideFrequencies[i] /= frequencySum
    for j:=0; j<5; j++ {
      currentNode.conversionRatios[i][j] /= rowSums[i]
    }
  }
  MutateNucleotideFrequencies(currentNode.leftChild, nucleotideMutationRate)
  MutateNucleotideFrequencies(currentNode.rightChild, nucleotideMutationRate)
}

/*----------------------------------------------------------------------------------------------------
 * An external library has been used for drawing samples from the gamma distribution
 *--------------------------------------------------------------------------------------------------*/
func GammaDistrubution(alpha, beta float64) float64 {
  var gammaDist = distuv.Gamma{Alpha: alpha, Beta: beta}
  return gammaDist.Rand()
}

/*----------------------------------------------------------------------------------------------------
 * Function used for changing the tree structures, when initiated, a subtree is picked at random, a new
 * location in the left over tree is chosen and the sub tree is attached at this site.
 *---------------------------------------------------------------------------------------------------*/
func MutateTopology(solution *node, rate float64, numSpecies int) *node {
  chance := rand.Float64()
  if chance < rate {
    totalNodes := 2*numSpecies - 1
    reqProb := float64(float64(1)/float64(totalNodes))
    var randSubtree *node
    randSubtree = nil

    for randSubtree == nil {
      randSubtree = PickRandomSubtree(solution, reqProb)
      if randSubtree == nil || randSubtree.parent == nil {
        randSubtree = nil
      }
    }
    var alteredBranchLength float64
    solution, alteredBranchLength = RemoveAndRestructureTree(solution, randSubtree)
    solution = MergeSubTrees(solution, randSubtree, numSpecies, alteredBranchLength)
  }
  return solution
}

/*----------------------------------------------------------------------------------------------------
 * Quite similar to the topology mutation, except that the random subtree is selected from a different
 * tree and the current tree is reorganized to accomadate the new subTree
 *---------------------------------------------------------------------------------------------------*/
func PerformCrossOver(solution *node, population []*node, recombinationProbability float64, numSpecies int) *node {
  chance := rand.Float64()
  if chance < recombinationProbability {
    numSolutions := len(population)
    secondParentIndex := rand.Intn(numSolutions)
    secondParent := population[secondParentIndex]

    totalNodes := 2*numSpecies - 1
    reqProb := float64(float64(1)/float64(totalNodes))
    var randSubtree *node
    randSubtree = nil
    for randSubtree == nil {
      randSubtree = PickRandomSubtree(secondParent, reqProb)
      if randSubtree == nil || randSubtree.parent == nil {
        randSubtree = nil
      }
    }
    // A fresh subtree is generated in order to avoid any complications with pointers and not adversly
    // affect the parent tree in case it is selected again for recombination.
    var alteredBranchLength float64
    randSubtree, alteredBranchLength = GenerateFreshSubTreeCopy(randSubtree)
    speciesSubList := CaptureSpecies(randSubtree, make([]string, 0))
    solution = RemoveSpecies(solution, solution, speciesSubList)
    solution = MergeSubTrees(solution, randSubtree, numSpecies, alteredBranchLength)
  }
  return solution
}

/*----------------------------------------------------------------------------------------------------
 * Given a tree, the function recursively searches for a node to be selected under the given probability
 * It has been generally ensured the probability reflects upon the total number of nodese  in the graph
 * Might return nil nodes at times, should be run multiple times until a proper node is returned
 *---------------------------------------------------------------------------------------------------*/
func PickRandomSubtree(root *node, prob float64) *node {
  if root == nil {
    return nil
  }
  var subTree *node
  chance := rand.Float64()
  if chance < prob {
    return root
  }

  searchDirection := 0.5
  chance = rand.Float64()
  if chance < searchDirection {
    subTree = PickRandomSubtree(root.leftChild, prob)
  } else {
    subTree = PickRandomSubtree(root.rightChild, prob)
  }
  if subTree != nil {
    return subTree
  }

  if chance < searchDirection {
    subTree = PickRandomSubtree(root.rightChild, prob)
  } else {
    subTree = PickRandomSubtree(root.leftChild, prob)
  }
  return subTree
}

/*----------------------------------------------------------------------------------------------------
 * Given the location of a node in the tree, makes the necessary adjustments to retain the binary structure
 * of the parent tree after its subgraph has been removed. We return the removed branch lenght as it has
 * to be used again while inserting the subtree again
 *---------------------------------------------------------------------------------------------------*/
func RemoveAndRestructureTree(root, randSubtree  *node) (*node, float64) {
  immParent := randSubtree.parent
  randSubtree.parent = nil
  var parentDist float64
  var direction int
  if immParent.leftChild == randSubtree {
    immParent.leftChild = nil
    parentDist = immParent.leftChildDistance
    immParent.leftChildDistance = 0
    direction = 0
  } else {
    immParent.rightChild = nil
    parentDist = immParent.rightChildDistance
    immParent.rightChildDistance = 0
    direction = 1
  }

  var edgeWeight float64
  // A special rearrangement is to be done if the parent of the random subtree is the root as the gneneral
  // assumption believes the existence of a grandparent node
  if immParent == root {
    if direction == 0 {
      root = immParent.rightChild
    } else {
      root = immParent.leftChild
    }
    root.parent = nil
    immParent = nil
  // The general case for rearrangement
  } else {
    if direction == 0 {
      immParent.rightChild.parent = immParent.parent
      edgeWeight = immParent.rightChildDistance
      if immParent == (immParent.parent.leftChild) {
        immParent.parent.leftChild = immParent.rightChild
        immParent.parent.leftChildDistance = (immParent.parent.leftChildDistance + edgeWeight)/2.0
      } else {
        immParent.parent.rightChild = immParent.rightChild
        immParent.parent.rightChildDistance = (immParent.parent.rightChildDistance + edgeWeight)/2.0
      }
    } else {
      immParent.leftChild.parent = immParent.parent
      edgeWeight = immParent.leftChildDistance
      if immParent == (immParent.parent.leftChild) {
        immParent.parent.leftChild = immParent.leftChild
        immParent.parent.leftChildDistance = (immParent.parent.leftChildDistance + edgeWeight)/2.0
      } else {
        immParent.parent.rightChild = immParent.leftChild
        immParent.parent.rightChildDistance = (immParent.parent.rightChildDistance + edgeWeight)/2.0
      }
    }
    immParent = nil
  }

  return root, parentDist
}


/*----------------------------------------------------------------------------------------------------
 * Given an imcomplete tree and a subTree, picks a random location on the incomplete tree and merges
 * the subtree at the given location
 *---------------------------------------------------------------------------------------------------*/
func MergeSubTrees(root, randSubtree *node, numSpecies int, alteredBranchLength float64) *node{
  var newRandLocation *node
  newRandLocation = nil
  reqProb := float64(1.0/float64(2*numSpecies - 1))
  for newRandLocation == nil {
    newRandLocation = PickRandomSubtree(root, reqProb)
  }
  root = ReorganizeTreeBranches(root, randSubtree, newRandLocation, alteredBranchLength)
  return root
}

/*-----------------------------------------------------------------------------------------------------
 * Given a subgraph and a  partial tree, attaches the subtree at  a random location in the partial tree
 * to complete it
 *---------------------------------------------------------------------------------------------------*/
func ReorganizeTreeBranches(root, randSubtree, newRandLocation *node, parentDist float64) *node  {
  if newRandLocation.leftChild == nil || newRandLocation.rightChild == nil {
    if newRandLocation.parent != nil {
      newRandLocation = newRandLocation.parent
    } else {
      var insertionNode node
      insertionNode.Initialize()
      insertionNode.name = "Ancestor"
      insertionNode.leftChild = randSubtree
      insertionNode.rightChild = newRandLocation
      randSubtree.parent, newRandLocation.parent = &insertionNode, &insertionNode
      insertionNode.leftChildDistance, insertionNode.rightChildDistance = parentDist, (rand.Float64()/10.0)
      return &insertionNode
    }
  }

  var insertionNode node
  insertionNode.name = "Ancestor"
  insertionNode.leftChild = randSubtree
  randSubtree.parent = &insertionNode
  insertionNode.parent = newRandLocation
  for i:=0; i<5; i++ {
    insertionNode.nucleotideFrequencies[i] = newRandLocation.nucleotideFrequencies[i]
    for j:=0; j<5; j++ {
      insertionNode.conversionRatios[i][j] = newRandLocation.conversionRatios[i][j]
    }
  }
  direction := rand.Float64()
  if direction < 0.5 {
    insertionNode.rightChild = newRandLocation.leftChild
    newRandLocation.leftChild.parent = &insertionNode
    newRandLocation.leftChild = &insertionNode
    insertionNode.leftChildDistance, insertionNode.rightChildDistance = parentDist, newRandLocation.leftChildDistance
  } else {
    insertionNode.rightChild = newRandLocation.rightChild
    newRandLocation.rightChild.parent = &insertionNode
    newRandLocation.rightChild = &insertionNode
      insertionNode.leftChildDistance, insertionNode.rightChildDistance = parentDist, newRandLocation.rightChildDistance
  }
  return root
}

/*----------------------------------------------------------------------------------------------------
 * Generates a copy of a subtree given its node pointer.
 *---------------------------------------------------------------------------------------------------*/
func GenerateFreshSubTreeCopy(subTree *node) (*node, float64) {

  var currNode node
  currNode.Initialize()
  if subTree == nil {
    return nil, 0
  } else {
    currNode.name = subTree.name
    currNode.leftChild, _ = GenerateFreshSubTreeCopy(subTree.leftChild)
    currNode.leftChildDistance = subTree.leftChildDistance
    currNode.rightChild, _ = GenerateFreshSubTreeCopy(subTree.rightChild)
    currNode.rightChildDistance = subTree.rightChildDistance
    if currNode.leftChild != nil {
      currNode.leftChild.parent = &currNode
    }
    if currNode.rightChild != nil {
      currNode.rightChild.parent = &currNode
    }
    for i:=0; i<5; i++ {
      currNode.nucleotideFrequencies[i] = subTree.nucleotideFrequencies[i]
      for j:=0; j<5; j++ {
        currNode.conversionRatios[i][j] = subTree.conversionRatios[i][j]
      }
    }
  }

  var parentDist float64
  parentNode := subTree.parent
  if parentNode.rightChild == subTree {
    parentDist = parentNode.rightChildDistance
    } else {
      parentDist = parentNode.leftChildDistance
    }
  return &currNode, parentDist
}

/*----------------------------------------------------------------------------------------------------
 * Given a node, captures all the dufferent species that this paritcular node is an ancestor of
 *---------------------------------------------------------------------------------------------------*/
func CaptureSpecies(randSubtree *node, speciesSubList []string) []string {
  if randSubtree == nil {
    fmt.Println("Invalid subTree given as input for SpeciesCapturing")
    os.Exit(1)
  }
  if randSubtree.leftChild == nil && randSubtree.rightChild == nil {
    speciesSubList = append(speciesSubList, randSubtree.name)
    return speciesSubList
  } else {
    speciesSubList = CaptureSpecies(randSubtree.leftChild, speciesSubList)
    speciesSubList = CaptureSpecies(randSubtree.rightChild, speciesSubList)
  }
  return speciesSubList
}

/*----------------------------------------------------------------------------------------------------
 * Given a iist of species to be deleted, the function removes and reogranized the graph after deleting
 * all the required leafs in a recursive manner
 *---------------------------------------------------------------------------------------------------*/
func RemoveSpecies(root, currNode *node, speciesSubList []string) *node {
  if currNode == nil || root == nil {
    fmt.Println("Invalid subTree given as input for sub-species Deletion")
    fmt.Println(NewickFormatTreeRepresentation(root))
    fmt.Println(NewickFormatTreeRepresentation(currNode))
    os.Exit(1)
  }
  speciesLen := len(speciesSubList)
  if currNode.leftChild == nil && currNode.rightChild == nil {
    matchFound := false
    for i:=0; i<speciesLen; i++ {
      if currNode.name == speciesSubList[i] {
        matchFound = true
        break
      }
    }
    if matchFound{
      root, _ = RemoveAndRestructureTree(root, currNode)
    }
  } else {
    root = RemoveSpecies(root, currNode.leftChild, speciesSubList)
    root = RemoveSpecies(root, currNode.rightChild, speciesSubList)
  }

  return root
}
