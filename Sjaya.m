
% Sjaya

function [BestValue,XTarget, BestCost]= Sjaya(nPop,VarMin,VarMax,nVar, fobj, MaxIt)

%%Input parameters
%%fobj----------------objective function
%%nPop---------------population size 
%%nVar---------------the number of variables
%%VarMin-------------the lower boundaries of variables
%%VarMin-------------the upper boundaries of variables
%%MaxIt--------------the maximum number of iterations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Output parameters
%%BestCost-----------convergence curve
%%BestValue----------the optimal fitness value
%%XTarget------------the optimal solution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  FES = 0; 
dim = nVar;
npop = nPop;


for i=1:nPop
    X(i,:)=VarMin+(VarMax-VarMin).*rand(1,nVar); 
end

for i=1:nPop
    f(i) = fobj( X(i,:) );
    FES = FES + 1;
end








%%  Main Loop
gen=1;
[BestCost(1),ind]=min(f);
XTarget=X(ind,:);
now_npop = npop;


while(gen +1<= MaxIt)

%     if FES >=  50000
%         disp it；
%         disp(gen);
%                     
%         break;
%     end

k = FES / 50000;



[f, sortedIndex] = sort(f);
X = X(sortedIndex, :);









    col = dim;
    Best=X(1,:);
    worst=X(now_npop,:);
    xnew=zeros(now_npop,col);

    
    for i=1:now_npop
        k1 = FES/50000;

        a1 = 1;
        a2 = 1;
        
        if rand < (1-k1)/2
            a1 = -1;
        end

        if rand < (1-k1)/2
            a2 = -1;
        end       



        rn1 = randn + 0.5;
        rn2 = randn + 0.5;


        
        fi=rand;

        if fi<=1/3
            xnew(i,:)=X(i,:)+ a1*rn1*(Best- X(i,:) )-  a2*rn2*(worst-X(i,:));
        elseif fi <= 2/3
            temp_X_index = select_from_history(npop, k, 3);
            temp_a = temp_X_index(1);
            temp_b = temp_X_index(2);
            temp_c = temp_X_index(3);
            xnew(i,:)=X(i,:)+ a1*rn1*(  X(temp_a,:)- X(i,:) ) - a2*rn2*( X(temp_b,:)- X(temp_c,:) ); 
        else
            k =  gen / MaxIt;
            X_A = X(1:floor( (1-k) * npop * 0.5));
            X_B = X(ceil(0.5*(1+k)*npop ) : npop );
            xnew(i,:)=X(i,:)+ a1*rn1*(mean(X_A)-X(i,:))-a1*rn1*(mean(X_B)-X(i,:)); 
        end
    
    end



    for i=1:npop
        xnew(i,:) = max(min(xnew(i,:),VarMax),VarMin);
        fnew(i) = fobj(xnew(i,:));
        FES = FES + 1;
    end
    
    for i=1:nPop
        if(fnew(i)<f(i))
            X(i,:) = xnew(i,:);
            f(i) = fnew(i);
        end
    end
    gen = gen+1;
    [BestCost(gen),ind]=min(f);
    XTarget=X(ind,:);

    

end
BestValue=min(f);
%%

end

function selected_indices = select_from_history(npop, k, select_number)
    % 计算根据高斯分布选择历史种群个体的概率
    mu = npop * (1 - k); % 均值
    sigma = npop / 6; % 标准差，这里选择 6 是基于经验值
    % 生成正态分布概率密度函数
    pdf = normpdf(1:npop, mu, sigma);
    % 归一化概率密度函数
    pdf = pdf / sum(pdf);
    
    % 从概率密度函数中根据不重复概率随机选择一半个体的索引
    selected_indices = zeros(1, select_number);
    for i = 1:(select_number)
        selected_indices(i) = randsample(1:npop, 1, true, pdf);
        pdf(selected_indices(i)) = 0; % 将已选中的概率设为0，确保不重复选择
        pdf = pdf / sum(pdf); % 重新归一化概率密度函数
    end
end