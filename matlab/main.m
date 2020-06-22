%% Analyzing data from barcoding-experiment 
% This is the main pipeline for processing barcoded data
% It starts from the expression matrix with cells as rows and
% genes (opt.),barcodes,timepoints,treatments and states (opt.) as columns.

%% Step 1: Input data
load('expressionMatrices.mat')

%% Step 2a: Assign states to cells (optional)
% This step can be skipped if state assignment has been made separately
% Consult stateAssignment for input requirements and options

states=stateassignment(expMatOriginal(:,1:end-2),'gene');
expMatOriginal=[expMatOriginal states];

%% Step 2b: Construct the original network

% Create the X and design matrix
[X,design]=datatobarcodedata(expMatOriginal);

% Manually create the samplefractions matrix
samplefractions=[1 7 0.8; 1 14 0.8; 1 21 1];

% Other input variables
sparsity_joint=1; 
sparsity_individual=1; %only used for treatment exp but needs to be set
tolerated_deviation_from_prior=100000; %prior weight
Apri=zeros(6,6); %prior matrix

[AhatOrig,~]=netfrombarcode(X,design,samplefractions,sparsity_joint,sparsity_individual,tolerated_deviation_from_prior,Apri);

figure(1)
plotmynetwork(AhatOrig,0.01);

%% Step 3: Construct the treatment network

[X,design]=datatobarcodedata(expMat);

samplefractions=[1 0 0.7;1 8 1*(1/3);2 0 0.7;2 8 1*(1/3);3 0 0.7;3 8 1*(1/3)]; 
sparsity_joint=1.1; 
sparsity_individual=0.9; % closer to zero gives few therapy-specific links
Apri=AhatOrig;
tolerated_deviation_from_prior=0.01; 

[Ahat,dAhat]=netfrombarcode(X,design,samplefractions,sparsity_joint,sparsity_individual,tolerated_deviation_from_prior,Apri);

figure(2)
plotmynetwork(Ahat,0.01);

%% Step 4: Calculate and plot growth rates
growth=sum(Ahat);
growthBMP4=sum(Ahat+dAhat{2});
growthTMZ=sum(Ahat+dAhat{3});

figure(3)
growths=[growth' growthBMP4' growthTMZ'];
bar(growths)

%% Step 5: Optimize eigenvalues (optional)

% Find the minimal changes to the A-matrix (dA) that leads to only
% negative real eigenvalues. 

%Condition 2 is from eq 11 in Zavlanos et al

A=Ahat+dAhat{2};
sparsity=0.3;
cvx_begin quiet
    variables dA(6,6)
    minimize norm((A+dA)-A,'fro')
        subject to
            sum(sum(abs((dA)-diag(diag((dA))))))<= sparsity
            diag((A+dA))<=-sum(abs(((A+dA)-diag(diag((A+dA))))'))'
cvx_end

