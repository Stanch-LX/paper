%_____________________________________________________________________________________________ %
%  Butterfly Optimization Algorithm (BOA) source codes demo V1.0                               %
%                                                                                              %
%  Author and programmer: Sankalap Arora                                                       %
%                                                                                              %
%         e-Mail: sankalap.arora@gmail.com                                                     %
%                                                                                              %
%  Main paper: Sankalap Arora, Satvir Singh                                                    %
%              Butterfly optimization algorithm: a novel approach for global optimization	   %
%              Soft Computing, in press,                                                       %
%              DOI: https://doi.org/10.1007/s00500-018-3102-4                                  %
%___________________________________________________________________________________________   %
%
function [fmin,best_pos,Convergence_curve]=BOA(n,N_iter,Lb,Ub,dim,fobj)
display('BOA is now tackling your problem');
tic
% n is the population size
% N_iter represnets total number of iterations
p=0.2;                       % probabibility switch
power_exponent=0.1;
sensory_modality=0.01;
c=(Lb+Ub)/2;
wMax=0.001;
wMin=0.001;
F=0.2;
CR=0.3;
%Initialize the positions of search agents
Sol=initialization(n,dim,Ub,Lb);
%Sol=Lb+(Ub-Lb).*rand(n,dim);
for i=1:n,
    Fitness(i)=fobj(Sol(i,:));
end

% Find the current best_pos
[fmin,I]=min(Fitness);
best_pos=Sol(I,:);
S=Sol; 

% Start the iterations -- Butterfly Optimization Algorithm 
for t=1:N_iter,
  
        for i=1:n, % Loop over all butterflies/solutions
         
          %Calculate fragrance of each butterfly which is correlated with objective function
          Fnew=abs(fobj(S(i,:)));
          FP=(sensory_modality*(Fnew^power_exponent));   
          w=wMin+t*((wMax-wMin)/N_iter);
         % w=-(exp(t/N_iter))^2+0.1;
        if rand<F
          %Global or local search
         if rand<p,    
               A=Levy(1.5,dim);
            dis = rand * rand * best_pos - Sol(i,:);        %Eq. (2) in pape
            S(i,:)=Sol(i,:)+dis*FP;
            So(i,:)=Lb+Ub-S(i,:);
            if S(i,:)<c
                S1(i,:)=c+(So(i,:)-c)*rand;
            else
                S1(i,:)=So(i,:)+(c-So(i,:))*rand;
            end
              fitness(i)=fobj(S(i,:));
              fitness1(i)=fobj(S1(i,:));
              if fitness(i)<fitness1(i)
                  S(i,:)=S(i,:);
              else
                  S(i,:)=S1(i,:);
              end
            
           else
              % Find random butterflies in the neighbourhood
              epsilon=rand;
              JK=randperm(n);
              A=Levy(1.5,dim);
                  dis=epsilon*epsilon*Sol(JK(1),:)-Sol(JK(2),:);
                  S(i,:)=w*Sol(i,:)+dis*FP;          
         end
          
        else
            if abs(fobj(Sol(i,:))-fobj(best_pos))>t||rand<CR  %可以让100动态地变化
                for j=1:size(Sol,2)
                   ran1=ceil(size(Sol,2)*rand);
                   ran2=ceil(size(Sol,2)*rand);
                   ran3=ceil(size(Sol,2)*rand);
                   while(ran1==ran2||ran2==ran3||ran1==ran3)
                        ran1=ceil(size(Sol,2)*rand);
                        ran2=ceil(size(Sol,2)*rand);
                        ran3=ceil(size(Sol,2)*rand);
                   end
                   if  j==ran1||j==ran2||j==ran3
                       S(i,j)=best_pos(j);
                   else 
                       JK=randperm(n);
                       S(i,j)=Sol(JK(1),j);
                   end
                end
                
            else
               JK=randperm(n);
                for j=1:size(Sol,2)
                 S(i,j)=Sol(i,j)+0.5*(Sol(JK(2),j)-Sol(JK(3),j));
                end
                 
            end
           for j=1:size(Sol,2)
                  ran1=ceil(size(Sol,2)*rand);
                  if  j==ran1
                     A=Levy(1.5,1);
                     S(i,j)=A*Sol(i,j);
                  end
                      
                   
           end
                   
        end
           
            % Check if the simple limits/bounds are OK
            %S(i,:)=simplebounds(S(i,:),Lb,Ub);
           Flag4ub=S(i,:)>Ub;
        Flag4lb=S(i,:)<Lb;
        S(i,:)=(S(i,:).*(~(Flag4ub+Flag4lb)))+Ub.*Flag4ub+Lb.*Flag4lb;
            % Evaluate new solutions
            Fnew=fobj(S(i,:));  %Fnew represents new fitness values
            
            % If fitness improves (better solutions found), update then
            if (Fnew<=Fitness(i)),
                Sol(i,:)=S(i,:);
                Fitness(i)=Fnew;
            end
           
           % Update the current global best_pos
           if Fnew<=fmin,
                best_pos=S(i,:);
                fmin=Fnew;
           end
         end
            
         Convergence_curve(t,1)=fmin;
         
         %Update sensory_modality
        %  sensory_modality=sensory_modality_NEW(sensory_modality, N_iter);
end
toc

% Boundary constraints
function s=simplebounds(s,Lb,Ub)
  % Apply the lower bound
  ns_tmp=s;
  I=ns_tmp<Lb;
  ns_tmp(I)=Lb;
  
  % Apply the upper bounds 
  J=ns_tmp>Ub;
  ns_tmp(J)=Ub;
  % Update this new move 
  s=ns_tmp;

  
function y=sensory_modality_NEW(x,Ngen)
y=x+(0.025/(x*Ngen));



