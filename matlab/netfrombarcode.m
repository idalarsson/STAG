
function [net_joint,nets_individual,I]=netfrombarcode(X,design,samplefractions,sparsity_joint,sparsity_individual,tolerated_deviation_from_prior,Apri)

% INPUT
% -----
% X={xij} is a matrix of states (rows 1,2,...k) times barcodes (columns)
% each element xij is the number of cells that are in state i and carry
% barcode j. IMPORTANT: each sample is represented by the same barcodes
% in the same order.
%
% design is a matrix of
% row 1 = barcoding experiment index
% row 2 = barcoding experiment timepoint
%
% samplefractions is a numerical matrix with the structure
%
% EXP   TIME    
% 1     1       0.8
% 1     2       0.8
% 1     3       1.0
% 2     1       0.75 etc
% where the first column is the experiment index
% and the second column is the timepoint index
% and the third column is the fraction f of cells that went to 10X analysis
% We assume that 1-f is the fraction of cells that were kept in culture
%
% OUTPUT
% ------
% net_joint is a k*k matrix of parameters, estimated as the average
% model across all conditions
% nets_individual is a 3d matrix of k*k*t parameters where each
% layer contains parameters for treatment t
%
%
% 


% Step 1: input sanity check

if(length(unique(hist(design(1,:),1:max(design(1,:)))))~=1)
    errorVar('Your data X contains samples with different size. Each barcode should be present in each sample in the same order');
end

% Step 2: get relevant constants that define the size of the data
num_conditions=max(design(1,:));
num_samplefractions=size(samplefractions,1);
num_barcodes=size(X,2)/num_samplefractions;
num_states=size(X,1);
num_treatments=max(samplefractions(:,1));

% Step 3: find the interval of column indices in X that defines each sample
for i=1:num_samplefractions
    I(i,:)=[(i-1)*num_barcodes+1 i*num_barcodes];
end

% Step 4: define a string that expresses the error function specific to
% this problem, corresponding to global error function in manuscript
errorVar='';
for i=1:num_treatments
    f=find(samplefractions(:,1)==i); % f points to all data for one treatment
    for g=2:length(f)
        i1=f(g-1);                   % i1 points to 'before' time point
        i2=f(g);                     % i2 points to 'after' time point
        I_before=[I(i1,1) I(i1,2)];  %
        I_after=[I(i2,1) I(i2,2)];
        t_before=samplefractions(i1,2);
        t_after=samplefractions(i2,2);
        dt=t_after-t_before;
        constant=num2str(dt*(1-samplefractions(i1,3))/(samplefractions(i1,3))*samplefractions(i2,3),5);
        constant2=num2str((1-samplefractions(i1,3))/(samplefractions(i1,3))*samplefractions(i2,3),5);
        err=['norm(X(:,' num2str(I_after(1)) ':' num2str(I_after(2)) ') - (I*' constant2 ' + ' constant '*(A+dA' num2str(i) '))*X(:,' num2str(I_before(1)) ':' num2str(I_before(2)) '),''fro'')'];
        errorVar=[errorVar '+' err];
    end
end
errorVar=errorVar(2:end);

%

% Step 5: express optimisation of the global error function as a CVX script
% which is stored in the code file cvxrun. 


vari=['A(' num2str(num_states) ',' num2str(num_states) ')'];
I=eye(num_states);
for i=1:num_treatments
    vari=[vari ' dA' num2str(i) '(' num2str(num_states) ',' num2str(num_states) ')'];
end

unix('rm cvxrun.m');
fid=fopen('cvxrun.m','w');
fprintf(fid,'cvx_begin quiet\n','cvx_begin');
fprintf(fid,'variables %s\n',vari);
fprintf(fid,'minimize %s\n',errorVar);
fprintf(fid,'subject to\n');

fprintf(fid,'A-diag(diag(A))>=0\n');

%fprintf(fid,'sum(sum(abs(A-diag(diag(A)))))<sparsity_joint\n',errorVar);

% this is where we define the style of the solution

%fprintf(fid,'sum(sum(alpha*(abs(A-diag(diag(A))))+(1-alpha)*((A-diag(diag(A))-Apri-diag(diag(Apri))).^2)))<sparsity_joint\n');

fprintf(fid,'sum(sum(abs(A-diag(diag(A)))))<= sparsity_joint\n');

fprintf(fid,'norm(A-Apri,''fro'') <= tolerated_deviation_from_prior \n');

fprintf(fid,'dA1==0\n');
for i=2:num_treatments
    %tag=['dA' num2str(i)];
    tag2=['(A + dA' num2str(i) ')'];
    fprintf(fid,'%s-diag(diag(%s)) >= 0 \n',tag2,tag2);
    fprintf(fid,'sum(sum(abs(%s - diag(diag(%s))))) <= sparsity_individual\n',tag2,tag2);
end
fprintf(fid,'cvx_end\n','cvx_begin');
fclose(fid)



% Step 6: run the CVX script. It will find the joint matrix A, and 
% each treatment specific deltaA1, deltaA2, ... assuming that deltaA1=0
% because that is the negative control (baseline). 

cvxrun

% Step 7: return the output from the function


net_joint=A;
nets_individual={};
for i=1:num_treatments
        tag=['dA' num2str(i)];
    eval(['nets_individual{i}=' tag ';'])
end

