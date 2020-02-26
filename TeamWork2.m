%{
load Sparse1Pair_psb420.dist;
data = Sparse1Pair_psb420;
%}

format longEng % Highest precision 

files = ["Sparse1Pair_psb420.dist", "SparseNoisy1Pair_psb420.dist";
         "Sparse2Pair_psb420.dist", "SparseNoisy2Pair_psb420.dist";
         "Sparse3Pair_psb420.dist", "SparseNoisy3Pair_psb420.dist";
         "Sparse4Pair_psb420.dist", "SparseNoisy4Pair_psb420.dist";
         "Sparse5Pair_psb420.dist", "SparseNoisy5Pair_psb420.dist";
         "Sparse6Pair_psb420.dist", "SparseNoisy6Pair_psb420.dist";
         "Sparse7Pair_psb420.dist", "SparseNoisy7Pair_psb420.dist";
         "Sparse8Pair_psb420.dist", "SparseNoisy8Pair_psb420.dist";
         "Sparse9Pair_psb420.dist", "SparseNoisy9Pair_psb420.dist"];
     
Error_X = zeros(1, 10);
ErrorNoisy_X = zeros(1, 10);
Error_X_Small = zeros(1, 10);
ErrorNoisy_X_Small = zeros(1, 10);

% If you just want to test on 1 file, change from i=1:10 to i=x:x
for eVal=1:2
for i=1:9   
  for fileIndex=1:2
    clean_file = files(i,1);
    noisy_file = files(i,2);
        
    %% load data
    
    if fileIndex == 1
        disp(["Calculating Gram Matrices for ", clean_file])
        clean = fopen(clean_file,'r');
        clean_data = fscanf(clean, "%g");
        fclose(clean);
    else
        disp(["Calculating Gram Matrices for ", noisy_file])
        clean = fopen(noisy_file,'r');
        clean_data = fscanf(clean, "%g");
        fclose(clean);
    end 
    
    
    n = clean_data(1);

    m = clean_data(2);
    distances_clean = clean_data(3:end);
    distances_clean = reshape(distances_clean, 3, []);
    R_clean = zeros(n,n);

    for k=1:size(distances_clean, 2)
        il = distances_clean(1, k);
        jl = distances_clean(2, k);
        dl = distances_clean(3, k);
        R_clean(il, jl) = dl;
    end
 
    
    R_clean = R_clean.^2; 
    %% Part 1
    
    eps = [1 0.1];

    [~, idx] = sort(R_clean);
    smallest = idx(2,:)';
    theta = [[1:100]' smallest];
    e = zeros(n,1);
    size(distances_clean,2);
    
    cvx_begin sdp
    variable G(n,n) semidefinite;
    minimize trace(G)
    subject to
        G == G';
        G*ones(n,1) == zeros(n,1);
        for e=1:size(distances_clean,2)
            E=zeros(n,1);
            E(distances_clean(1,e),1)=1;
            E(distances_clean(2,e),1)=-1;
            abs(dot(G*E, E)-distances_clean(3,e)) <= eps(1, eVal)
        end 
    cvx_end

    %% Part 2
    % #### Estimating G using Alg. 1 ####

    % for clean data
    clean = fopen('Pair_psb420.dist','r');
    dataClean = fscanf(clean, "%g");
    fclose(clean);

    nClean = dataClean(1);
    distancesClean = dataClean(2:end);
    R = reshape(distancesClean, n, n);
    S = R.^2;

    % Following algorithim 1 for Gram matrix
    ones_vec = ones(n,1);

    rho_clean = 1/(2*n) * ones_vec' * S * ones_vec;
    v_clean = 1/n*(S*ones_vec - rho_clean*ones_vec);
    G_true = 1/2*v_clean*ones_vec' + 1/2*ones_vec*v_clean' - 1/2.*S;

    % for noisy data
    noisy = fopen('MeasuredPair_psb420.dist','r');
    dataNoisy = fscanf(clean, "%g");
    fclose(noisy);

    n = dataNoisy(1);
    distancesNoisy = dataNoisy(2:end);
    RN = reshape(distancesNoisy, n, n);
    SN = RN.^2;

    ones_vec = ones(n,1);

    rhoNoisy = 1/(2*n) * ones_vec' * SN * ones_vec;
    v_noisy = 1/n*(SN*ones_vec - rhoNoisy*ones_vec);
    G_noisy = 1/2*v_noisy*ones_vec' + 1/2*ones_vec*v_noisy' - 1/2.*SN;

        %% Part 3

        if eVal == 1
            if fileIndex == 1
                Error = norm(G-G_true, 'fro')
                Error_X(i) = Error;
            else
                ErrorNoisy = norm(G-G_noisy, 'fro')
                ErrorNoisy_X(i) = ErrorNoisy;
            end
        else 
            if fileIndex == 1
                Error = norm(G-G_true, 'fro')
                Error_X_Small(i) = Error;
            else
                ErrorNoisy = norm(G-G_noisy, 'fro')
                ErrorNoisy_X_Small(i) = ErrorNoisy;
            end 
        end
  end
end
end
%% Part 4
% We had some issues with the cvx modeling. It worked initially, then
% returned the same matrix for all G's calculated, hence the straight lines.
hold on 
plot(0:9, Error_X);
plot(0:9, ErrorNoisy_X);
hold off
figure
hold on 
plot(0:9, Error_X_Small);
plot(0:9, ErrorNoisy_X_Small);
hold off
%% Written Responses

% Question 1:
% One could certainly determine the decay rate, but our optimization
% is not following the correct rate of decay to do so. We could fit the
% errors of the fit to calculate the least squares. Matlabs poly-fit would
% most likely provide an adequate solution.

% Question 2
% A smaller epsilon should lead to an overall smaller error, as the
% tolerance is narrower. This comes at the cost of calculation time.

%Question 3
% It might take a long time but there should be a solutoin as the equation
% is convex and continous throughout (semi-definite matrix). However,
% computers may not be able to have high precision/handle floating point
% values well or complex numbers. 