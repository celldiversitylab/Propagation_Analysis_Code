% Propagation Analysis Code version 1.0 © GPL-3.0, 2019, ETH Zurich, Institute of Biochemistry, Ulrich Berge
function  [observedDivisionMatrix]=makeAbsolTransMatrix(clustered)
% 
% input:
% clustered: - matrix without column or row names
%     1. column: lineage ID (int)
%     2. column: cell ID in the form of 1 is the mother cell of the two daughters 11 and 12 etc. (int)
%     3. column: cell cycle duration
%     4. column: identifier of which kind a cell is an outlier
%         -3: S-cell
%         3: L-cell
%         1: the three remaining granddaughter set cells of an S- or L-cell
%         -4/4: 2:2 granddaughter set where the outlier cells are not siblings
%         -2/2: 2:2 granddaughter set where the outlier cells are siblings
% output:
% empirical 3x3 division matrix: counts how often each division type occurred
% % S-->S/0 S-->0/0 S-->L/0
% % 0-->S/0 0-->0/0 0-->L/0
% % L-->S/0 L-->0/0 L-->L/0
% % where S, L, and 0 refer to S-, L-, and 0-cells


%this analysis only focuses on 3:L and 3:S granddaugher sets; 2:2
%granddaughter sets are defined as 0-cells
clustered(find(clustered(:,4)~=3 & clustered(:,4)~=-3),4)=0;    

M=zeros(3,3);
for a = 1: max(clustered(:,1))                      %for each lineage
    temp = clustered(find(clustered(:,1)==a),:);
    for b = 1:length(temp(:,1))                     %find for each cell its two daughter cells...
         mom = num2str(temp(b,2));
         d1 = str2num([mom '1']);
         d2 = str2num([mom '2']);
         mom = str2num(mom);
         if length(find(temp(:,2)==mom))>0 & length(find(temp(:,2)==d1))>0 & length(find(temp(:,2)==d1))>0
             temp2=[temp(find(temp(:,2)==mom),4) temp(find(temp(:,2)==d1),4) temp(find(temp(:,2)==d2),4)];
             if temp2(1)==-3 & sum(temp2(2:3))== -3  %...and add up for each division type (3x3 matrix; here M)...
                 M(1,1)=M(1,1)+1;                    %...the number of occurrences
             elseif temp2(1)==-3 & sum(temp2(2:3))== 0    
                 M(1,2)=M(1,2)+1;
             elseif temp2(1)==-3 & sum(temp2(2:3))== 3    
                 M(1,3)=M(1,3)+1;
             elseif temp2(1)==0 & sum(temp2(2:3))== -3    
                 M(2,1)=M(2,1)+1;
             elseif temp2(1)==0 & sum(temp2(2:3))== 0    
                 M(2,2)=M(2,2)+1;
             elseif temp2(1)==0 & sum(temp2(2:3))== 3    
                 M(2,3)=M(2,3)+1;    
             elseif temp2(1)==3 & sum(temp2(2:3))== -3    
                 M(3,1)=M(3,1)+1;
             elseif temp2(1)==3 & sum(temp2(2:3))== 0    
                 M(3,2)=M(3,2)+1;
             elseif temp2(1)==3 & sum(temp2(2:3))== 3    
                 M(3,3)=M(3,3)+1;    
             end
         end
    end
end

observedDivisionMatrix=M;
end     


                 
                 
                 