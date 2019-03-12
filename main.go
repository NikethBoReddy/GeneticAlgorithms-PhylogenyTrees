package main

import (
  "fmt"
  "os"
  //"time"
  "math/rand"
)

// The basic object that is used to represent a species in the phyolgenetic tree model
type node struct {
  name string
  conversionRatios [5][5]float64
  nucleotideFrequencies [5]float64
  leftChild, rightChild, parent *node
  leftChildDistance, rightChildDistance float64
}

// Initialization of the nucleotide conversion ratios and frequencies
func (n *node) Initialize() {
    n.conversionRatios = [5][5]float64{
      //A    C      G     T     Deletion
      {0.94, 0.02, 0.01, 0.02, 0.01}, // A
      {0.02, 0.94, 0.02, 0.01, 0.01}, // C
      {0.01, 0.02, 0.94, 0.02, 0.01}, // G
      {0.02, 0.01, 0.02, 0.94, 0.01}, // T
      {0.01, 0.01, 0.01, 0.01, 0.96}, // Insertion
      }
    n.nucleotideFrequencies = [5]float64 {0.24, 0.24, 0.24, 0.24, 0.04}
}

// A basic reperesentation for the aligned sequences of the given species
type speciesGenome struct {
  name string
  nucleotideSequence string
}

func main() {

  if len(os.Args) != 2 {
    fmt.Println("Please enter the dataset filepath for starting the program")
    os.Exit(1)
  }
  filename := os.Args[1]
  speciesList, speciesMap  := LoadDatasets(filename)
  fmt.Println("The required nucleotide sequences of species has been successfully obtained!!")

  numSolutions := LoadIntConfig("ga.algo.params.population.count")
  numGenerations := LoadIntConfig("ga.algo.params.generations.count")
  minStableGenerations := LoadIntConfig("ga.algo.params.generations.stable.limit")

  initialPopulation := GenerateRandomSolutions(speciesList, numSolutions)
  fmt.Println("Successfully generated a set of random Phylogenetic Trees for the Genetic Algorithm")

  bestPhylogenyModel := RunGASimulations(initialPopulation, speciesMap, numGenerations, minStableGenerations)
  fmt.Println(NewickFormatTreeRepresentation(bestPhylogenyModel))
}

/*----------------------------------------------------------------------------------------------------------------
 * Function for randomly generataing tree topologies and branchlengths. Works by recursively joining any two nodes
 * that do not have a parent together until there is only one such node present which becomes the root of the graph.
 *----------------------------------------------------------------------------------------------------------------*/
func GenerateRandomSolutions(speciesList []speciesGenome, numSolutions int) []*node {
  numSpecies := len(speciesList)
  population := make([]*node, numSolutions)

  for i:=0; i<numSolutions; i++ {

    treeConstructionBase := make([]*node, numSpecies)
    for j:=0; j<numSpecies; j++ {
      var leaf node
      leaf.Initialize()
      leaf.name = speciesList[j].name
      treeConstructionBase[j] = &leaf
    }

    for len(treeConstructionBase) != 1 {
      var ancestralNode node
      ancestralNode.Initialize()
      ancestralNode.leftChild, treeConstructionBase = GetAndRemoveRandomNodePointer(treeConstructionBase)
      ancestralNode.leftChildDistance = (rand.Float64()/10)
      ancestralNode.leftChild.parent = &ancestralNode
      ancestralNode.rightChild, treeConstructionBase = GetAndRemoveRandomNodePointer(treeConstructionBase)
      ancestralNode.rightChildDistance = (rand.Float64()/10)
      ancestralNode.rightChild.parent = &ancestralNode
      ancestralNode.name = "Ancestor"

      treeConstructionBase = append(treeConstructionBase, &ancestralNode)
    }

    population[i] = treeConstructionBase[0]
  }

  return population
}

/*-----------------------------------------------------------------------------------------------------------
 * Helper method for the above fucntion; performs the process of choosing a random node and maitaing the list
 * of parentless nodes at any point of time while the tree is being built
 *------------------------------------------------------------------------------------------------------------*/
func GetAndRemoveRandomNodePointer(treeConstructionBase []*node) (*node, []*node) {
  if len(treeConstructionBase) < 1 {
    fmt.Println("Invalid tree construction procedure, please double-check!!")
    os.Exit(1)
  }
  randNodeIndex := rand.Intn(len(treeConstructionBase))
  retNodePointer := treeConstructionBase[randNodeIndex]
  treeConstructionBase = append(treeConstructionBase[:randNodeIndex], treeConstructionBase[randNodeIndex+1:]...)
  return retNodePointer, treeConstructionBase
}
