%Takes an input matrix and compares alpha values from .1-2 with the MyTAD
%method. Uses the compare value as the matrix to compare it to

function MyTADFindBestAlpha(input, compare)
for i = 1:20    
    global matrixX;
    global euclidMatrix;
    global exceptionThrown;
    
    %reading in matrix
    matrixX = readmatrix(input);

    try
        alphaLevel = alpha;
    catch exception
        fprintf("Alpha level not set, setting to 1");
        alphaLevel = i/10;
    end

    %obtaining distance from IF and alpha value
    matrixRaw = 1./(matrixX.^alphaLevel);
    dSparse = sparse(matrixRaw);
    
    %running Johnson's shortest path algorithm
    d = graphallshortestpaths(dSparse);
 
    %Converting dist matrix and shortest path matrix to XYZ coordinates
    [xyzMatrixDist] = cmdscale(d, 3);
    
    %obtaining matrix to compare it to
    xyzMatrixRaw = readmatrix(compare);
    
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

    %creates arrays to derive best score from after program is run
    Pearson(i) = (rho1X+rho1Y+rho1Z)/3;
    Spearman(i) = (rho2X+rho2Y+rho2Z)/3;
    RMSEA(i) = RMSE;
    
end


