% Propagation Analysis Code version 1.0 © GPL-3.0, 2019, ETH Zurich, Institute of Biochemistry, Ulrich Berge
function [clustered] = ODM(data,threshold)

% input:
% data - matrix without column or row names
%     1. column: cell ID in the form of 1 is the mother cell of the two daughters 11 and 12 etc. (int)
%     2. column: real time of cell division/death
%     3. column: cell cycle duration
%     4. column: generation attribute of the cell within one lineage (int)
%     9. column: lineage ID (int)
% threshold - integer number
%     the threshold to classify cells into L-, S-, or 0-cells; S- and L-cells have the same
%     threshold

% output:
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

threshold=threshold-1;
clustered=[];
for l = 1:max(data(:,9)) %per tree
    tree = data(find(data(:,9)==l),:);
    if length(tree(:,1))>0
       temp=[];
       for a = 1:length(tree(:,1)) %finds granddaughter sets
                mom = tree(a,1);
                gd1 = str2num([num2str(mom), '11']);
                gd2 = str2num([num2str(mom), '12']);
                gd3 = str2num([num2str(mom), '21']);
                gd4 = str2num([num2str(mom), '22']);
                temp = [tree(find(tree(:,1)==gd1),3);...
                    tree(find(tree(:,1)==gd2),3);...
                    tree(find(tree(:,1)==gd3),3);...
                    tree(find(tree(:,1)==gd4),3)];
                nans = find(isnan(temp(:,1))==1);
                temp(nans)=[];
                if isempty(temp)==0
                    if length(temp(:,1))==4 %only if a complete granddaughter set is found...
                        temp = [[gd1;gd2;gd3;gd4] temp];
                        temp = sortrows(temp,2);
                        if length(unique(temp(:,2)))>1 %...perform the clustering
                            da = pdist(temp(:,2),'mahalanobis');
                            Z = linkage(da);
                            [H,T,outperm]=dendrogram(Z);
                            set(gcf,'visible','off')
                            Z(find(Z(:,3)>1),1:2);
                            [r c]=find(Z(:,1:3)==5);
                            Z(r,c)=str2num([num2str(Z(1,1)) num2str(Z(1,2))]);
                            [r c]=find(Z(:,1:3)==6);
                            Z(r,c)=str2num([num2str(Z(2,1)) num2str(Z(2,2))]);
                            if sum(Z(1:2,3))<=Z(3,3) & (Z(3,3)/Z(2,3))>threshold
                                sc = Z(find(Z(:,3)>1),1:2);
                                sc = sc(end,:);
                            else 
                                sc=[];
                            end
                            if length(sc)>0
                                sc1=[];sc2=[];
                                for d = 1:length(num2str(sc(1)))
                                    t = num2str(sc(1));
                                    sc1 = [sc1 str2num(t(d))];
                                end
                                for d = 1:length(num2str(sc(2)))
                                    t = num2str(sc(2));
                                    sc2 = [sc2 str2num(t(d))];
                                end
                                if length(sc1)==length(sc2)
                                    if  abs(diff(temp(sc1,1)))~=1   %2:2 outlier are not siblings
                                        if sum(temp(sc1,2)>temp(sc2,2))==2
                                            temp(sc1,3) = 4;
                                            temp(sc2,3) = -4;
                                        else
                                            temp(sc1,3) = 4;
                                            temp(sc2,3) = -4;
                                        end
                                    else                            %2:2 outliers are siblins
                                        temp(sc1,3) = 2;
                                        temp(sc2,3) = -2;
                                    end
                                elseif sc1==4                       %3:L outlier
                                    temp(sc1,3) = 3;
                                    temp(sc2,3) = 1;
                                elseif sc1==1                       %3:S outlier
                                    temp(sc1,3) = -3;
                                    temp(sc2,3) = 1;
                                end
                            else
                                temp(:,3)=0;
                            end
                            clustered=[clustered; repmat(l,4,1) temp];
                        else
                        clustered=[clustered; repmat(l,4,1) temp zeros(4,1)];
                        end
                    end
            end
        end
end
end