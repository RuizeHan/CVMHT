function [nodes, cost, timeSpent] = GMMCP_Solver_ADN(Net_Cost_Mat, NN, NC, isFindMax, method, maxIter, maxTime,showDebugInfo, Kcliques, dummyWeight)

% Variable List
% Net_Cost_Mat is the weights of the GMCP graph! Please include the dummy
% nodes in the Net_Cost_Mat. LAST NODE OF EACH CLUSTER IS THE DUMMY NODE
%
% NN: Number of Nodes in each Cluster (including the dummy node). It can be a single number
% (For Graphs with the same number of NNs in each Cluster) or can be a
% vector of numbers which shows number of members for each Cluster
%
% NC: Number of Clusters
%
% isFindMax: 1 if finding Maximum, 0 if finding Minimum
%
% method: 0 if Linear Programmin
%         1 if Binary Linear Programming
%         2 if mixed Linear Programming

method = 2;
% 
% if(isunix)
%     addpath('/MOT_toolbox/cplex/')
% end

Net_Cost_Mat = (Net_Cost_Mat + Net_Cost_Mat')/2;
Net_Cost_Mat = Net_Cost_Mat - dummyWeight;

infin = 1000000;

if ~exist('maxIter','var')
    maxIter = inf;
end

if ~exist('maxTime','var')
    maxTime = inf;
end

if ~exist('showDebugInfo','var')
    showDebugInfo = 0;
end

if length(NN) == 1
    NN = NN*ones(1,NC);
end

if length(NN) ~= NC
    warning('Length of NN needs to be consistant with the provided NC');
    keyboard;
end

if showDebugInfo
    disp('GMCP LP Solver has started making the matrixes');
end

startT = tic;

VN = sum(NN);% Number of all tracklets
NNcumsum = cumsum(NN); % Cumsum Number of  tracklets

Enum = zeros(VN,VN); % Number of edges
Vnum = zeros(NC, max(NN)); % Number of clusters x Number of max nodes number  0 : not a nodes

EN = 0;
ClusterLabel = [];
counterV = 0;  % Number of nodes
for i = 1:NC
    Vnum(i,1:NN(i)) = counterV+1:counterV+NN(i);
    counterV = counterV + NN(i);
    ClusterLabel = [ClusterLabel;i*ones(NN(i),1)];
    EN = EN + (VN-NC-(NN(i)-1))*(NN(i)-1);  %(NN(i)-1):number of no-dummy nodes in cluster i; (VN-NC-(NN(i)-1)): number of no-dummy nodes in other clusters except i
end
EN = EN/2;

if VN ~= size(Net_Cost_Mat,1)
    warning('Size of Net_Cost_Mat has to be the same as number of nodes.');
    keyboard;
    Net_Cost_Mat = Net_Cost_Mat(1:VN,1:VN);
end


E = 0;
fE = zeros(EN + NC, 1);
x0 = zeros(EN+ NC, 1);
ctype = char([ones(1,EN)*('B'*1), ones(1,NC)*('C'*1)]);
last = [4,1,2,3];
next = [2,3,4,1];

% The constraint # 3 
AeqR1 = 0;
Aeq1 = zeros(NC, EN + NC);
beq1 = ones(NC, 1)* Kcliques * 2;

for i = 1:VN
    Ci = ClusterLabel(i);
    if (i==NNcumsum(Ci))  % The last node in each cluster is dummy node
        Aeq1(Ci, EN+Ci) = 1; % last NC clos represent The dummy nodes
        continue;
    end
    for j = 1:VN
        
        Cj = ClusterLabel(j);
        if (j==NNcumsum(Cj))
            continue;
        end
        
        if Ci == Cj
            continue;
        end
        
        %  The energy function C
        if i<j &&  (Cj - 1 == mod(Ci,4) || Ci - 1 == mod(Cj,4))
            E = E + 1;
            Enum(i,j) = E; % Record the new position in the vector of the original position (i,j) in the matrix 
            Enum(j,i) = E;  % All the possible connect between two real nodes
            Elist(E,1) = i;
            Elist(E,2) = j;
            fE(E) = Net_Cost_Mat(i,j);
            if isnan(fE(E))
                fE(E) = - infin;
            end
        else
        end
        
        if Cj - 1 == mod(Ci,4)
            Aeq1(Ci,Enum(i,j)) = 1;
        end
    end
end

Aeq1_s = Aeq1;

for Ci = 1 : size(Aeq1,1) 
     Aeq1_s(Ci,1:EN) = Aeq1(Ci,1:EN) + Aeq1(last(Ci),1:EN) ;   
end


AR1 = 0;
A1 = sparse(NC*(NC-1)*(NC-2)/2*max(NN)^3, EN + NC);

AR2 = 0;AR3 = 0;
A2 = zeros(VN*2, EN+NC);
b2 = ones(VN*2,1);

% The constraint # 1 
for i = 1:4
   j = next(i);             
   m = next(j);      
   n = next(m);       
   for ii = 1:NN(i)-1
        for jj = 1:NN(j)-1
            for mm = 1:NN(m)-1
                for nn = 1:NN(n)-1
                AR1 = AR1 + 1;
                A1(AR1, Enum(Vnum(i,ii),Vnum(j,jj))) = 1;
                A1(AR1, Enum(Vnum(j,jj),Vnum(m,mm))) = 1;
                A1(AR1, Enum(Vnum(m,mm),Vnum(n,nn))) = 1;
                A1(AR1, Enum(Vnum(n,nn),Vnum(i,ii))) = -1;
                end
            end
        end
   end

end

% The constraint # 2
for m = 1: NC
    n = next(m);
    for mm = 1:NN(m)-1  % Cluster i Node ii
        AR2 = AR2 + 1;
        A2(AR2, Enum(Vnum(m,mm), Vnum(n,1:NN(n)-1))) = 1;
    end
    for nn = 1:NN(n)-1  
        AR2 = AR2 + 1;
        A2(AR2, Enum(Vnum(m,1:NN(m)-1), Vnum(n,nn))) = 1;
    end
end

A1 = A1(1:AR1,:);
b1 = ones(AR1, 1) * 2;

Aeq = [Aeq1_s];
beq = [beq1];

% A = [sparse(A2)];
% b = [b2];
A = [A1;sparse(A2)];
b = [b1;b2];

if showDebugInfo
    disp('Optimization has started');
end
    
    %x = solveMBIP(fE,A,b,Aeq,beq,ctype,x0,isFindMax);
    
        options = cplexoptimset();
    
        if isFindMax
            fE = - fE;
        end
        if showDebugInfo
            disp('everything is ready now, Optimization has started');
        end
        if method==1
            [x, fval, exitflag, output] = cplexbilp(fE, A, b, Aeq, beq, ...
                x0, options);
        elseif method ==0
            devideBy = 5;
            [x, fval, exitflag, output] = cplexlp(fE*50, A, b, Aeq, beq/devideBy, ...
                zeros(size(fE)),ones(size(fE))/devideBy,[ ], options);
        elseif method ==2
            [x, fval, exitflag, output] = cplexmilp(fE,A,b,Aeq,beq,[],[],[],zeros(size(fE)),[],ctype,x0,options);
        elseif method == 3
            %It has already calculated before...
        end

if any(abs(x-round(x))>10e-8)
    warning('There are some variables which differ at least 10e-12 from the integer numbers. Be careful about it, I will round the solution');
    keyboard;
end
x = round(x);

x2 = x(1:EN);
x2 = sort(find(x2==1));
x2 = Elist(x2,:);
x2 = unique(x2(:));

counter = 0;
nodes = zeros(Kcliques, NC);
while length(x2)
    counter = counter + 1;
    nodes(counter,ClusterLabel(x2(1))) = x2(1);
    tmp = Enum(x2(1),:);
    tmp = nonzeros(tmp);
    tmp = Elist(tmp(x(tmp)>0.9)',2);
    tmp2 = Enum(tmp(1),:);
    tmp2 = nonzeros(tmp2);
    tmp2 = Elist(tmp2(x(tmp2)>0.9)',2);   
    nodes(counter,ClusterLabel(tmp)) = tmp;
    nodes(counter,ClusterLabel(tmp2)) = tmp2;
    x2 = x2(~ismember(x2,nodes(counter,:)));
end

Net_Cost_Mat = Net_Cost_Mat + dummyWeight;
Net_Cost_Mat(isnan(Net_Cost_Mat)) = 0;

x2 = find(~ismember(1:VN, NNcumsum));
x2 = x2(~ismember(x2,nonzeros(nodes(:))));
[~,sInd] = sort(sum(Net_Cost_Mat(x2,:),2));
x2 = x2(sInd);

for i = counter + 1 : min(Kcliques, counter + length(x2))
    nodes(i,ClusterLabel(x2(i-counter))) = x2(i-counter);
end

for i = 1:NC
    zeroInd = (nodes(:,i) == 0);
    nodes(zeroInd, i) = NNcumsum(i);
end

cost = 0;
for i = 1:size(nodes,1)
    tmp = Net_Cost_Mat(nodes(i,:),nodes(i,:));
    tmp = tmp .* (1-eye(size(tmp)));
    cost = cost + sum(sum(tmp));
end
     cost = cost / 2;

timeSpent = toc(startT);

