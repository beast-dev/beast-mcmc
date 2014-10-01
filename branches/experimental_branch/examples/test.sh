#!/bin/bash
for file in release{/,/Benchmarks/,/calibrations/,/clockModels/,/continuousTraits/,/discreteTraits/,/EmpiricalCodonModels/,/microsatellite/,/Phylogeography/,/recombination/,/starBEAST/,/testXML/,/VirusPractical/}*\.xml 
do 
java -jar ../build/dist/beast.jar -overwrite $file 
done