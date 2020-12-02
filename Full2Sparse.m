%Program I wrote to take in a sparse matrix vector from matlab and converts
%it to a sparse matrix file, seperated by spaces. Input Matrix
%Sparse2Sparse(inputMatrixDirectory, filePathToSave)
function Full2Sparse(inputMat, filePath)

inputMatrix = readmatrix(inputMat);

fileID = fopen(filePath, 'w');
for i = 1:size(inputMatrix, 1)
    for j = i:size(inputMatrix, 2)
        if(inputMatrix(i, j)~=0)        
            fprintf(fileID, "%d %d %d\n", i, j, inputMatrix(i, j));
        end
    end
end
