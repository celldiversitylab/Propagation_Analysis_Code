% Propagation Analysis Code version 1.0 © GPL-3.0, 2019, ETH Zurich, Institute of Biochemistry, Ulrich Berge
function[observedDivisionMatrix, M2,pobserved, P]=PropagationAnalysis(data, permutations) 

% input:
% data - matrix without column or row names
%     1. column: cell ID in the form of 1 is the mother cell of the two daughters 11 and 12 etc. (int)
%     2. column: real time of cell division/death
%     3. column: cell cycle duration
%     4. column: generation attribute of the cell within one lineage (int)
%     9. column: lineage ID (int)
% permutations - number of intended permutations to be run (int)
% 
% output:
% each 3x3 division matrix:
% S-->S/0 S-->0/0 S-->L/0
% 0-->S/0 0-->0/0 0-->L/0
% L-->S/0 L-->0/0 L-->L/0
% where S, L, and 0 refer to S-, L-, and 0-cells
% 
% observedDivisionMatrix - 3x3
% M2 - cell with each element being the result of a permuted 3x3 division matrix;
%     is the basis for the estimated probability distribution (CDF)
% pobserved - probability of the observed division type count based on the 
%             estimated probability distribution 
% P - cell with each element being the result of a permuted 3x3 division matrix
%     is the estimation of the probability distribution function of obtaining 
%     any possible occurrence value per division type
    
    
rng % sets random number generator back to its default settings

d2 = ODM(data,2);
%generates the data as described in Figure 3b
[observedDivisionMatrix]=makeAbsolTransMatrix(d2);

%generates the data as described in Figure 3c
for i = 1:permutations
        D=[];
            for a = 1:max(data(:,4))    %randomise data from all lineages within each generation...
                d3=data(find(data(:,4)==a),:);
                n = find(isnan(d3(:,3))~=1);
                idx = randperm(length(n));
                resamp=d3(n(idx),3);
                d3(n,3)=resamp;
                D=[D; d3];
            end
        temp=ODM(D,1);                  %...and do the clustering for each permutation
        [M]=makeAbsolTransMatrix(temp); %...and create a division matrix
        M2{i}=M;
end

%create for each division type the probability distributions of the
%permutations
 T = sum(sum(observedDivisionMatrix));
    X=[];
    for a = 1:length(M2{1}(:,1))
        for b = 1:length(M2{1}(1,:))
            s = size(M2);
            temp = [];temp2=[];
            for c = 1:s(2)
                temp = [temp;  M2{c}(a,b)];
            end
            [x y]=hist(temp, 0:1:T);
            X=[X; x, a b];
        end
    end

% find in this probability distribution the observed value
    for a = 1:length(M2{1}(:,1))
        for b = 1:length(M2{1}(1,:))
            temp = X(find(X(:,end-1)==a & X(:,end)==b),1:end-2);
            tempreal = observedDivisionMatrix(a,b);
            total = sum(temp);
            pobserved(a,b) = (sum(temp(1:tempreal))+temp(tempreal+1)/2)/sum(temp);
        end
    end

%generates the data as described in Figure 3c
for j = 1:permutations
        D=[];
            for a = 1:max(data(:,4))
                d3=data(find(data(:,4)==a),:);
                n = find(isnan(d3(:,3))~=1);
                idx = randperm(length(n));
                resamp=d3(n(idx),3);
                d3(n,3)=resamp;
                D=[D; d3];
            end
        temp=ODM(D,1);
        [M]=makeAbsolTransMatrix(temp);
        REAL=M;

    for a = 1:length(M2{1}(:,1))
        for b = 1:length(M2{1}(1,:))
            temp = X(find(X(:,end-1)==a & X(:,end)==b),1:end-2);
            tempreal = REAL(a,b);
            total = sum(temp);
            p(a,b) = (sum(temp(1:tempreal))+temp(tempreal+1)/2)/sum(temp);
        end
    end
    P{j}=p;
end
    
    
    