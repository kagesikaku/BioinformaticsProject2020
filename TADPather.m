function TADPather(input)
    global matrixX;
    global euclidMatrix;
    global exceptionThrown;
    
    %reading in matrix
    matrixX = readmatrix(input);

    try
        alphaLevel = alpha;
    catch exception
        fprintf("Alpha level not set, setting to 1");
        alphaLevel = .8;
    end

    %obtaining distance from IF and alpha value
    matrixRaw = 1./(matrixX.^alphaLevel);

    %creating euclidean distance matrix for comparison later
    euclidMatrix = dist(matrixX);
    
    [rows cols] = size(matrixRaw);
    for i = 1:rows
        for j = 1:cols
            if (matrixRaw(i, j) == Inf)
                matrixRaw(i, j) = 1;
            end
        end
    end

    %normalizing based on 0-1 scaling, the equivalent of SCN
    %obtaining mix value in matrix
    miner = min(matrixRaw, [], 'all');
    %obtaining max
    maxer = max(matrixRaw, [], 'all');
    %normalizing the matrix to a scale of 0-1
    matrixNorm = rescale(matrixRaw, 'InputMin', miner, 'InputMax', maxer);

    %turning raw and normalized matrix into a graph
    graphRaw = graph(matrixRaw, 'upper');
    graphNorm = graph(matrixNorm, 'upper');

    %running Djikistra's shortest algorithm on both graphs
    d = distances(graphRaw, 'method', 'positive');

    
    %Converting dist matrix and shortest path matrix to XYZ coordinates
    [xyzMatrixDist] = XYZ(d);
    [xyzMatrixRaw] = XYZ(matrixRaw);
    
    inputMat = [xyzMatrixDist.X xyzMatrixDist.Y xyzMatrixDist.Z];
    inputMat2 = [xyzMatrixRaw.X xyzMatrixRaw.Y xyzMatrixRaw.Z];
    
    rhoTot = corr(inputMat, inputMat2);
    %Pearson correlations
    rho1X = corr(xyzMatrixDist.X, xyzMatrixRaw.X, 'Type', 'Pearson');
    rho1Y = corr(xyzMatrixDist.Y, xyzMatrixRaw.Y, 'Type', 'Pearson');
    rho1Z = corr(xyzMatrixDist.Z, xyzMatrixRaw.Z, 'Type', 'Pearson');

    %Spearman correlations, does not work if negative
    try
        rho2X = corr(xyzMatrixDist.X, xyzMatrixRaw.X, 'Type', 'Spearman');
        rho2Y = corr(xyzMatrixDist.Y, xyzMatrixRaw.Y, 'Type', 'Spearman');
        rho2Z = corr(xyzMatrixDist.Z, xyzMatrixRaw.Z, 'Type', 'Spearman');

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
    mat2pdb(xyzMatrixDist);

    %RMSE calculation based on XYZ coordinates.
    RMSE = sqrt(mean(((xyzMatrixDist.X-xyzMatrixRaw.X)+(xyzMatrixDist.Y-xyzMatrixRaw.Y)+(xyzMatrixDist.Z-xyzMatrixRaw.Z)).^2));

    %writing to log file
    fid = fopen('tadPather.log', 'w');
    fprintf(fid, "Input File: %s \n", input);
    fprintf("Input File: %s \n", input);
    fprintf(fid, "Convert Factor: %d \n", alphaLevel);
    fprintf("Convert Factor: %d \n", alphaLevel);
    fprintf(fid, "RMSE: %d \n", RMSE);
    fprintf("RMSE: %d \n", RMSE);
    fprintf(fid, "AVG Spearman correlation Dist vs. Reconstructed Dist: %d \n", (rho2X+rho2Y+rho2Z)/3);
    fprintf("AVG Spearman correlation Dist vs. Reconstructed Dist: %d \n", (rho2X+rho2Y+rho2Z)/3);
    fprintf(fid, "AVG Pearson correlation Dist vs. Reconstructed Dist: %d \n", (rho1X+rho1Y+rho1Z)/3);
    fprintf("AVG Pearson correlation Dist vs. Reconstructed Dist: %d \n", (rho1X+rho1Y+rho1Z)/3);

    
%Shrec3D XYZ coordinates calculator
function [xyzObj] = XYZ(d)
    %converting to xyz coordinates, taken from Shrec3D
    n=length(d);
    center=zeros(n,1);
    for i=1:n
        for j=1:n
        center(i)=center(i)+d(i,j)^2-1/n*(d(j,j:n)*d(j:n,j));
        end
    end
    center(center<0)=0;
    center=sqrt(center/n);
    distmat=1/2*(repmat(center.*center,1,n)+repmat((center.*center)',n,1)-d.*d);
    [V,D]=eigs(distmat,4);
    lambda=diag(D);
    XYZ=[V(:,1)*sqrt(lambda(1)) V(:,2)*sqrt(lambda(2)) V(:,3)*sqrt(lambda(3))];
    scale=100/max(max(XYZ));
    XYZ=XYZ*scale;

    xyzObj.X = XYZ(:, 1);
    xyzObj.Y = XYZ(:, 2);
    xyzObj.Z = XYZ(:, 3);

    end

end


