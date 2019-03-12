package main

import (
  "strconv"
  "strings"
  "os"
  "github.com/fredericlemoine/gotree/io/newick"
	"github.com/fredericlemoine/gotree/tree"
  "github.com/fredericlemoine/gotree/draw"
)


func PrintNewickFormatTree(newickTree string) {
  var t *tree.Tree
  var err error
  t, err = newick.NewParser(strings.NewReader(newickTree)).Parse()
  if err != nil {
 	  panic(err)
  }

  width := LoadIntConfig("ga.output.draw.width")
  height := LoadIntConfig("ga.output.draw.height")

  d := draw.NewTextTreeDrawer(os.Stdout, width, height, 0)
  l := draw.NewNormalLayout(d, true, true, false, true)
  l.DrawTree(t)
}

/*--------------------------------------------------------------------------------------------------------------------------
 * A standard formta that is used for representing phyolgenetic trees in a recursive fashion using commas to separate children
 * and brackets to represnt an ancestral node. The output of this function can be plugges into online application that draw
 * the phylogenetic tree structure
 *------------------------------------------------------------------------------------------------------------------------*/
func NewickFormatTreeRepresentation(root *node) string {
  if root == nil {
    return ""
  }
  var newickFormat string
  if root.name == "Ancestor" {
    newickFormat =  "("
  }
  if root.leftChild == nil && root.rightChild == nil {
    newickFormat += root.name
  } else {
    newickFormat +=  NewickFormatTreeRepresentation(root.leftChild)  + ":" + strconv.FormatFloat(root.leftChildDistance, 'f', 5, 64) + ","
    newickFormat +=  NewickFormatTreeRepresentation(root.rightChild) + ":" + strconv.FormatFloat(root.rightChildDistance, 'f', 5, 64)
  }
  if root.name == "Ancestor" {
    newickFormat += ")"
  }
  if root.parent == nil {
    newickFormat += ";"
  }
  return newickFormat
}
