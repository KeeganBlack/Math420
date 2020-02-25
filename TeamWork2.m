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
    disp("--------------------------------------------------------------------------------")
    clean_file = files(i,1);
    noisy_file = files(i,2);
        
    disp(["Calculating Gram Matrices for ", clean_file, "and ", noisy_file, ": "])
    %% load data
    if fileIndex == 1
        clean = fopen(clean_file,'r');
        clean_data = fscanf(clean, "%g");
        fclose(clean);
    else
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
    
    E = zeros(n,n);
    [~, idx] = sort(R_clean);
    smallest = idx(2,:)';
    theta = [[1:100]' smallest];
    
    %{
    Note:
    - |<Ge_{i,j}, e_{i,j}> - d_{i,j}^2| is what is used in the slides. The
      values did not seem equal when i played around with dummy values.
    - The e_{i,j} used is the point closest to point i. I dont think that
      this is the correct way to do this, so we will need to fix this.
    %}
    
    disp("#### Minimizing trace G ####")
        disp(["** Using tolerance: ", eps(eVal)])
        cvx_begin sdp
        variable G(n,n) semidefinite;
        minimize trace(G);
        subject to
            G == G';
            G*ones(n,1) == zeros(n,1);
            for k=1:n
                idx = theta(k,:);
                e = zeros(n,1);
                e(idx(1)) = 1;
                e(idx(2)) = -1;
                d = R_clean(idx);
                abs(dot(G*e, e)- d) <= eps(1,eVal);
            end   
        cvx_end
    
    %% Part 2
    disp("#### Estimating G using Alg. 1 ####")
    
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

   
%     
%     rho_clean = 1/(2*n) * ones_vec' * R_clean * ones_vec;
%     v_clean = 1/n*(R_clean*ones_vec - rho_clean*ones_vec);
%     G_true = 1/2*v_clean*ones_vec' + 1/2*ones_vec*v_clean' - 1/2.*R_clean;
    
    
    
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
    
    disp("#### Calculating Error ####")
    Error = norm(G-G_true, 'fro')
    ErrorNoisy = norm(G-G_noisy, 'fro')
   
    if eVal == 1
        if fileIndex == 1
         Error_X(i) = Error;
        else
        ErrorNoisy_X(i) = ErrorNoisy;
        end
    else 
        if fileIndex == 1
            Error_X_Small(i) = Error;
        else
            ErrorNoisy_X_Small(i) = ErrorNoisy;
        end 
    end
    
    
    disp("--------------------------------------------------------------------------------")
  end
end
end
%% Part 4

hold on 
plot(0:9, Error_X);
plot(0:9, ErrorNoisy_X);
hold off
figure
hold on 
plot(0:9, Error_X_Small);
plot(0:9, ErrorNoisy_X_Small);
hold off
