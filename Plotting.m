% Propagation Analysis Code version 1.0 © GPL-3.0, 2019, ETH Zurich, Institute of Biochemistry, Ulrich Berge
function Plotting(observedDivisionMatrix, M2,pobserved, P)


% input:
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



    T = sum(sum(observedDivisionMatrix));
    X=[];d=1;L=[];

figure(1)
        for a = 1:3
            for b = 1:3
                s = size(M2);
                temp = [];
                for c = 1:s(2)
                    temp = [temp;  M2{c}(a,b)];
                end
                [x y]=hist(temp, 0:1:T);
                xrel = x./sum(x);xcum=[];
                for c = 1:length(y)
                    xcum(c,:) = sum(xrel(1:c-1))+0.5*xrel(c);
                end
                subplot(3,3,d)
                plot(y,xcum, 'o-'), hold on
                plot(observedDivisionMatrix(a,b),pobserved(a,b), 'rx')
                axis([0 T 0 1])
                text(1,0.5, num2str(pobserved(a,b)))
                d=d+1;
            end
        end
        

 figure(2)
      X=[];d=1;L=[];L2=[];
         for a = 1:3
            for b = 1:3
                s = size(P);
                temp = [];
                for c = 1:s(2)
                    temp = [temp;  P{c}(a,b)];
                end
                [x y]=hist(temp, 0:0.01:1);
                X=[X; x, a b];
                [xreal yreal]=hist(pobserved(a,b), 0:0.01:1); xreal2 = find(xreal>0);
                subplot(3,3,d)
                bar(y, x/sum(x))
                hold on
                bar(yreal, xreal./5, 'r')
                axis([-0.1 1.1 0 1])
                d=d+1;
                if pobserved(a,b) >=0.5
                    L(a,b)=length(find(temp>=pobserved(a,b)))/length(temp)
                else
                    L(a,b)=length(find(temp<=pobserved(a,b)))/length(temp)
                end
                text(0.1,0.2, num2str(L(a,b)))
            end
         end  
    

    l2=0;
    for c = 1:length(P)
        temp = [P{c}];
            l(c)=length(find(temp<0.025 | temp>0.975))
        if l(c) > 0
            l2 = l2+1;
        end
    end

    figure(3)
    [x y]=hist(l, 0:5);
    bar(y, x./sum(x))
    for a = 1:length(y)
        xs = sum(x(1:a))/sum(x)
        text(a-1,0.9,num2str(xs))
    end
    ylim([0 1])
end
        
        
        
        
        
        
        
        
        
        
        
        