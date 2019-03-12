package main

import(
  "fmt"
  "os"
  "bufio"
  "strings"
  "strconv"
)

var stringConfigs map[string]string
var floatConfigs map[string]float64
var intConfigs map[string]int
var isLoaded bool

/*-------------------------------------------------------------------------------
 * A basic function for loading the information of the aligned gene sequences for
 * the phylogenetic tree reconstruction.
 *------------------------------------------------------------------------------*/
func LoadDatasets(filename string) ([]speciesGenome, map[string]speciesGenome) {
  file, err := os.Open(filename)
  if err != nil {
    fmt.Println("Something went wrong while trying to open the input file:" + filename)
    os.Exit(1)
  }
  scanner := bufio.NewScanner(file)
  scanner.Scan()
  speciesCount, err := strconv.Atoi(scanner.Text())
  if err != nil {
    fmt.Println("Invalid format of the input file; needs to begin with a number describing the number of aligned sequences")
    os.Exit(1)
  }

  speciesList := make([]speciesGenome, speciesCount)
  speciesMap := make(map[string]speciesGenome)
  for i:=0; i<speciesCount; i++ {
    scanner.Scan()
    lineInfo := strings.Split(scanner.Text(), " ")
    if len(lineInfo) != 2 {
      fmt.Println("Invalid format of the input file; Additional space characters found in line:" + strconv.Itoa(i+2))
      os.Exit(1)
    }
    speciesList[i] = speciesGenome{name:lineInfo[0], nucleotideSequence:lineInfo[1]}
    speciesMap[lineInfo[0]] = speciesGenome{name:lineInfo[0], nucleotideSequence:lineInfo[1]}
  }
  file.Close()
  if scanner.Err() != nil {
    fmt.Println("Something went wrong while trying to read the input file:" + filename)
    os.Exit(1)
  }
  return speciesList, speciesMap
}

/*---------------------------------------------------------------------------------
 * Since there have been a lot of variable parameter for the algorithm, I have decided
 * that moving forward with a config style approach would be more convenient rather than
 * using the traditional command line approach
 *------------------------------------------------------------------------------------*/
func LoadConfigs() {
  file, err := os.Open("config.txt")
  if err != nil {
    fmt.Println("Something went wrong while trying to read the input file:config.txt")
    os.Exit(1)
  }
  stringConfigs = make(map[string]string)
  floatConfigs  = make(map[string]float64)
  intConfigs    = make(map[string]int)

  scanner := bufio.NewScanner(file)
  for scanner.Scan(){
    property := strings.Split(scanner.Text(),"=")
    if len(property) != 2 {
      fmt.Println("Invalid format for configFile of line:" + scanner.Text())
      os.Exit(1)
    }
    val := strings.Split(property[1], ",")
    if len(val) != 2 {
      fmt.Println("Invalid format for configFile of line:" + scanner.Text())
      os.Exit(1)
    }

    switch val[1] {
    case "int":
      intVal, err := strconv.Atoi(val[0])
      if err != nil {
        fmt.Println("Error while writing int value for: " + scanner.Text())
      }
      intConfigs[property[0]] = intVal

    case "float64":
      floatVal, err := strconv.ParseFloat(val[0], 64)
      if err != nil {
        fmt.Println("Error while writing int value for: " + scanner.Text())
      }
      floatConfigs[property[0]] = floatVal

    case "string":
      stringConfigs[property[0]] = val[0]
    }
  }
  file.Close()
  if scanner.Err() != nil {
    fmt.Println("Something went wrong while trying to read the input file:config.txt")
    os.Exit(1)
  }
  fmt.Println("Succesfully loaded data from config file into program")
}

/*--------------------------------------------------------------------------------------
 * Smaller set of functions to fetch different value types
 *-------------------------------------------------------------------------------------*/
func LoadIntConfig(property string) int {
  if !isLoaded {
    LoadConfigs()
    isLoaded = true
  }
  retVal, exists := intConfigs[property]
  if !exists {
    fmt.Println("Invalid int property requested: " + property)
    os.Exit(1)
  }
  return retVal
}

func LoadFloatConfig(property string) float64 {
  if !isLoaded {
    LoadConfigs()
    isLoaded = true
  }
  retVal, exists := floatConfigs[property]
  if !exists {
    fmt.Println("Invalid int property requested: " + property)
    os.Exit(1)
  }
  return retVal
}

func LoadStringConfig(property string) string {
  if !isLoaded {
    LoadConfigs()
    isLoaded = true
  }
  retVal, exists := stringConfigs[property]
  if !exists {
    fmt.Println("Invalid int property requested: " + property)
    os.Exit(1)
  }
  return retVal
}
