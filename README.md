# BioinformaticsProject2020
ShortestDistanceTadFinder V1.0

//Input Matrix is an n*n matrix based on Interaction Frequency in a chromosome. alpha is a value between .1 and 2
//impacting the distance matrix created. compare is the XYZ matrix of the actual structure to use for correlation comparisons
MyTAD('input matrix destination', alpha, compare)
IE MyTAD('regular90.txt', .1, 'regularstructre.txt')
//saves data to 'currentdirectory/output/simulated_output/MyTAD'

//Full2Sparse takes a full matrix and creates a sparse matrix 
//Full2Sparse(inputMatrix, fileDestination)

//MyTADFindBestAlpha finds the best alpha depending on the n*n input matrix and a comparison matrix location
//MyTADFindBestAlpha(inputMatrixLocation, compareMatrixLocation)

//CompareOthers takes an input XYZ matrix, a comparison XYZ matrix, an alpha value, a start location for the input XYZ matrix, and an output location
//CompareOthers(inputXYZmat, compareXYZmat, startColInput, outputLocationFile, optionalAlpha)
//saves data to 'currentdirectory/output/simulated_output/MyTAD'


