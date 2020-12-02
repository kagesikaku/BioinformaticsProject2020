%
%Format of MyTAD
%TADPather(input matrix directory, alpha value, comparison file
%Primarily used to have MyTAD run against an n*n matrix as input
%outputs PDB and .log to output/simulated_output/MyTAD/
%
function MyTAD(input, alpha, compare)

    global matrixX;
    global exceptionThrown;
    
    %breaking down files to save later
    [filepath, name, ext] = fileparts(input);
    
    %reading in matrix
    matrixX = readmatrix(input);

    %trying to set alpha
    try
        alphaLevel = alpha;
    catch exception
        fprintf("Alpha level not set, setting to 1");
        alphaLevel = 1;
    end

    %obtaining distance from IF and alpha value
    matrixRaw = 1./(matrixX.^alphaLevel);
    
    %converting to sparse matrix
    dSparse = sparse(matrixRaw);
    
    %running Johnson's shortest path algorithm
    d = graphallshortestpaths(dSparse);
 
    %Converting dist matrix and shortest path matrix to XYZ coordinates
    [xyzMatrixDist] = cmdscale(d, 3);

    %Reading actual structure matrix
    xyzMatrixRaw = readmatrix(compare);
    
    %getting XYZ from Method and valid struct
    myxyzMatrixDist.X = xyzMatrixDist(:, 1);
    myxyzMatrixDist.Y = xyzMatrixDist(:, 2);
    myxyzMatrixDist.Z = xyzMatrixDist(:, 3);
    
    myxyzMatrixRaw.X = xyzMatrixRaw(:, 1);
    myxyzMatrixRaw.Y = xyzMatrixRaw(:, 2);
    myxyzMatrixRaw.Z = xyzMatrixRaw(:, 3);
    
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

    %Converting the matrix to a pdb, from 3DMax
    mat2pdb(myxyzMatrixDist);

    %RMSE calculation based on XYZ coordinates.
    RMSE = sqrt(mean(((myxyzMatrixDist.X-myxyzMatrixRaw.X)+(myxyzMatrixDist.Y-myxyzMatrixRaw.Y)+(myxyzMatrixDist.Z-myxyzMatrixRaw.Z)).^2));

    string1 = sprintf('utput/simulated_output/MyTad/%s.log', name);
    
    %writing to log file
    fid = fopen(string1, 'w');
    
    fprintf(fid, "Input File: %s \n", input);
    fprintf("Input File: %s \n", input);
    fprintf(fid, "Convert Factor: %d \n", alpha);
    fprintf("Convert Factor: %d \n", alpha);
    fprintf(fid, "RMSE: %d \n", RMSE);
    fprintf("RMSE: %d \n", RMSE);
    fprintf(fid, "AVG Spearman correlation Dist vs. Reconstructed Dist: %d \n", (rho2X+rho2Y+rho2Z)/3);
    fprintf("AVG Spearman correlation Dist vs. Reconstructed Dist: %d \n", (rho2X+rho2Y+rho2Z)/3);
    fprintf(fid, "AVG Pearson correlation Dist vs. Reconstructed Dist: %d \n", (rho1X+rho1Y+rho1Z)/3);
    fprintf("AVG Pearson correlation Dist vs. Reconstructed Dist: %d \n", (rho1X+rho1Y+rho1Z)/3);
    
    %printing out PDB files
    string2 = sprintf('simulated_output/MyTad/%s.pdb', name);
    movefile('mat2PDB.pdb', string2);
    
end

