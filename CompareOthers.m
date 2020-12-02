%input is Compare('directory of sparsematrix1', 'directory of true
%structure', where the actual XYZ coordinates start in the matrix)

function CompareOthers(mat1, mat2, start, folderName, alpha)

    x = str2num(start);
    MatrixTry = readmatrix(mat1);
    MatrixTrue = readmatrix(mat2);

    [filepath, name, ext] = fileparts(mat1);
    
    myxyzMatrixDist.X = MatrixTry(:, x);
    myxyzMatrixDist.Y = MatrixTry(:, x+1);
    myxyzMatrixDist.Z = MatrixTry(:, x+2);
    
    myxyzMatrixRaw.X = MatrixTrue(:, 1);
    myxyzMatrixRaw.Y = MatrixTrue(:, 2);
    myxyzMatrixRaw.Z = MatrixTrue(:, 3);
    
    inputMat = [myxyzMatrixDist.X myxyzMatrixDist.Y myxyzMatrixDist.Z];
    inputMat2 = [myxyzMatrixRaw.X myxyzMatrixRaw.Y myxyzMatrixRaw.Z];
    
    
    rhoTot = corr(inputMat, inputMat2);
    %Pearson correlations
    rho1X = corr(myxyzMatrixDist.X, myxyzMatrixRaw.X, 'Type', 'Pearson');
    rho1Y = corr(myxyzMatrixDist.Y, myxyzMatrixRaw.Y, 'Type', 'Pearson');
    rho1Z = corr(myxyzMatrixDist.Z, myxyzMatrixRaw.Z, 'Type', 'Pearson');

    %Spearman correlations, does not work if negative
    try
        rho2X = corr(myxyzMatrixDist.X, myxyzMatrixRaw.X, 'Type', 'Spearman');
        rho2Y = corr(myxyzMatrixDist.Y, myxyzMatrixRaw.Y, 'Type', 'Spearman');
        rho2Z = corr(myxyzMatrixDist.Z, myxyzMatrixRaw.Z, 'Type', 'Spearman');

        %catching if the inputs are complex, if so then skip the Spearman
        %correlation
    catch exception
        disp('exception thrown for having complex inputs, no spearman correlation');
        rho2X = 0;
        rho2Y = 0;
        rho2Z = 0;
        exceptionThrown = 1; 
    end

    %RMSE calculation based on XYZ coordinates.
    RMSE = sqrt(mean(((myxyzMatrixDist.X-myxyzMatrixRaw.X)+(myxyzMatrixDist.Y-myxyzMatrixRaw.Y)+(myxyzMatrixDist.Z-myxyzMatrixRaw.Z)).^2));

    string1 = sprintf('Simulated_output/%s/%s.log',folderName, name);
    string2 = sprintf('Simulated_output/%s/%s.pdb',folderName, name);
    %Converting the matrix to a pdb, from 3DMax
    mat2pdb(myxyzMatrixDist);
    
    %writing to log file
    fid = fopen(string1, 'w');
    
    %writing to log file
    fprintf(fid, "Input File: %s \n", mat1);
    fprintf("Input File: %s \n", mat1);
    try
        fprintf(fid, "Convert Factor: %d \n", alpha);
        fprintf("Convert Factor: %d \n", alpha);
    catch exception
        fprintf(fid, "Convert Factor: NA \n");
        fprintf("Convert Factor: NA \n");
    end
    fprintf(fid, "RMSE: %d \n", RMSE);
    fprintf("RMSE: %d \n", RMSE);
    fprintf(fid, "AVG Spearman correlation Dist vs. Reconstructed Dist: %d \n", (rho2X+rho2Y+rho2Z)/3);
    fprintf("AVG Spearman correlation Dist vs. Reconstructed Dist: %d \n", (rho2X+rho2Y+rho2Z)/3);
    fprintf(fid, "AVG Pearson correlation Dist vs. Reconstructed Dist: %d \n", (rho1X+rho1Y+rho1Z)/3);
    fprintf("AVG Pearson correlation Dist vs. Reconstructed Dist: %d \n", (rho1X+rho1Y+rho1Z)/3);

    movefile('mat2PDB.pdb', string2);
end