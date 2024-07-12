%%%%%%%%%%%%%%%%%%%%%%%%
%           WOA-PID source codes              %
%                                                             %
%     Developed in MATLAB R2019a          %
%                                                             %
%   Author and programmer: Qin Jibiao    %
%                                                             %
%         e-Mail:321829896@qq.com          %
%%%%%%%%%%%%%%%%%%%%%%%%

global yd y timef

SearchAgents_no=30;%SearchAgents_no：种群规模
dim=3;%dim维数：优化参数个数
ub = [20,10,10];%ub优化变量上限
lb = [0,0,0];%lb优化变量下限
Max_iter = 100;%Max_iter最大迭代次数

% The Whale Optimization Algorithm - PID
Bsj=0;

% initialize position vector and score for the leader
Leader_pos=zeros(1,dim);
Leader_score=Fitness_Func(Leader_pos,Bsj); %change this to -inf for maximization problems


%Initialize the positions of search agents
Positions=initialization(SearchAgents_no,dim,ub,lb);

Convergence_curve=zeros(1,Max_iter);

t=1;% Loop counter

% Main loop
while t<Max_iter
    time(t)=t;
    for i=1:size(Positions,1)
        
        % Return back the search agents that go beyond the boundaries of the search space
        Flag4ub=Positions(i,:)>ub;
        Flag4lb=Positions(i,:)<lb;
        Positions(i,:)=(Positions(i,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;
        Bsj=Fitness_Func(Positions(i,:),Bsj);
        
        % Calculate objective function for each search agent        
        fitness=Fitness_Func(Positions(i,:),Bsj);

        
        % Update the leader
        if fitness<Leader_score % Change this to > for maximization problem
            Leader_score=fitness; % Update alpha
            Leader_pos=Positions(i,:);
        end
        
    end
    
    a=2-t*((2)/Max_iter); % a decreases linearly fron 2 to 0
    
    % a2 linearly dicreases from -1 to -2 to calculate
    a2=-1+t*((-1)/Max_iter);
    
    % Update the Position of search agents 
    for i=1:size(Positions,1)
        r1=rand(); % r1 is a random number in [0,1]
        r2=rand(); % r2 is a random number in [0,1]
        
        A=2*a*r1-a;  
        C=2*r2;      
        
        
        b=1;              
        l=(a2-1)*rand+1;   
        p = rand();        
        
        for j=1:size(Positions,2)
            
            if p<0.5   
                if abs(A)>=1
                    rand_leader_index = floor(SearchAgents_no*rand()+1);
                    X_rand = Positions(rand_leader_index, :);
                    D_X_rand=abs(C*X_rand(j)-Positions(i,j)); 
                    Positions(i,j)=X_rand(j)-A*D_X_rand;      
                    
                elseif abs(A)<1
                    D_Leader=abs(C*Leader_pos(j)-Positions(i,j)); 
                    Positions(i,j)=Leader_pos(j)-A*D_Leader; 
                end
                
            elseif p>=0.5
              
                distance2Leader=abs(Leader_pos(j)-Positions(i,j));
                
                Positions(i,j)=distance2Leader*exp(b.*l).*cos(l.*2*pi)+Leader_pos(j);
                
            end
            
        end
    end
    Bsj=Leader_score;
    BestS=Leader_pos;
    
    BestS
    kp(t)=BestS(1);
    ki(t)=BestS(2);
    kd(t)=BestS(3);
    
    Bsj
    BsJ_l(t)=Bsj;
     
    t=t+1;

end

disp('kp,ki,kd');
BestS

%绘图
figure(1);%创建绘图窗口
plot(timef,yd,'r',timef,y,'c:','linewidth',2);%绘图参数
xlabel('Time(s)');ylabel('yd,y');%坐标轴
legend('STEP','WOA-PID');%图例

figure(2);
plot(time,BsJ_l,'r','linewidth',2);
xlabel('Times','FontSize',8);ylabel('Best J');

figure(3);
plot(time,kp,'r',time,ki,'g',time,kd,'b','linewidth',2);
xlabel('Time(s)');ylabel('kp,ki,kd');
legend('kp','ki','kd');


function Positions=initialization(SearchAgents_no,dim,ub,lb)

Boundary_no= size(ub,2); % numnber of boundaries

% If the boundaries of all variables are equal and user enter a signle
% number for both ub and lb
if Boundary_no==1
    Positions=rand(SearchAgents_no,dim).*(ub-lb)+lb;
end

% If each variable has a different lb and ub
if Boundary_no>1
    for i=1:dim
        ub_i=ub(i);
        lb_i=lb(i);
        Positions(:,i)=rand(SearchAgents_no,1).*(ub_i-lb_i)+lb_i;       
    end
end
end
