function[states]=stateAssignment(expMat,scoreMethod)

%INPUT
%-----
%expMat: Should be table with cells as rows and genes as columns.
%
%scoreMethod: One of 'kmeans', 'geneset' or 'gene'. Will ask for user
%input, make sure to write this in quotation marks (e.g. 'GDF15'), except
%for number of states for kmeans. 
%
%OUTPUT
%-----
%states: vector with state assignment

%Step1:Define states
if strcmp(scoreMethod,'kmeans')
    prompt='How many states? ';
    x=input(prompt);
    [a,b]=size(expMat);
    y=zeros(a,b-3);
    for i=1:b-3
        i;
        y(:,i)=expMat.(expMat.Properties.VariableNames{i+1});
    end
    %cluster using kmeans and add to matrix
    f=find(sum(y)>100);
    %[u,s,v]=svd(y(:,f),'econ');
    %rng(1);
    states=kmeans(y(:,f),x,'replicates',20);
    %states_v2=kmeans(y(:,f)-u(:,1)*s(1,1)*v(:,1)',5,'replicates',10);
    states=array2table(states);
    %expMat=[expMat, states];
elseif strcmp(scoreMethod,'geneset')
    prompt='What is the name and location of your geneset-file? ';
    x=input(prompt);
    geneSets=readtable(x);
    prompt2='What direction is your genesets? ';
    y=input(prompt2);
    %Detta fungerar inte än
    [~,b]=size(geneSets);
    [c,d]=size(expMat);
    d=(d-3);
    scoreFrame=array2table(zeros(c,b));
    for i=1:b
        score=genesetScore(expMat(:,1:d),geneSets(:,i),y);
        scoreFrame(:,i)=score(:,1);
        scoreFrame.Properties.VariableNames{i}=char(score.Properties.VariableNames);
    end
    scoreFrame.Properties.RowNames=score.Properties.RowNames;
    states=num2cell(zeros(c,2));
    for j=1:c
        [~,I]=max(table2array(scoreFrame(i,:)));
        states(i,1)=scoreFrame.Properties.RowNames{i};
        states(i,1)=scoreFrame.Properties.VariableNames{I};
    end
    %expMat=[expMat states];
elseif strcmp(scoreMethod,'gene')
    prompt = 'What gene should define the states? ';
    x = input(prompt);
    prompt2 = 'Specify mean or quantile ';
    y= input(prompt2);
    states=singleGeneScore(expMat,x,y);
    %expMat=[expMat states];
end

end 
