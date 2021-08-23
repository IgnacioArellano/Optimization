function Genetic_Algoritm_2D_test
%% mínimo U(x*) = 0   x*=[0,0]
U = @(x) 20 + x(1).^2 -10*cos (2*pi.*x(1))+ x(2).^2-10*cos(2*pi.*x(2));  % multimodal function
U2 = @(x1,x2) 20 + x1.^2 -10*cos (2*pi.*x1)+ x2.^2-10*cos(2*pi.*x2);


%%initial parameters
n_of_v=2;
generation_max = 1000000; % maximum number of generations
size_popu = 30; % population size
range(1,:) = [-5.12 5.12]; % variable bound
range(2,:) = [-5.12 5.12];
lengthS(1) = 20; % number of bits
lengthS(2) = 20;
method = 1;  % 1 for roulette and  ~1 for Tournament
cross_prob = 0.9; % crossover probability
multi_crossover = 0;% use multi-crossover
mut_prob = 0.1; % mutation probability
tourni_flag = 1; % use roulette wheel
epsilon = 1e-7; % function tolerance
flag = 0; % stall generations flag
flag1 = 0; % scalin flag
stall_gen = 500; % stall generations for termination
n_of_c = 0; % for constraint handling


% INITIALIZATION OF STRINGS
genotype=cell(1,n_of_v);
for i = 1:n_of_v
    genotype(i)= {randi([0,1],[size_popu,lengthS(i),])};
end

fenotype = cell(n_of_v,1);

for generation=1:generation_max
    % Decoded and normalized value
    for i = 1:2
        fenotype{i} = bi2de(genotype{i},'left-msb');
        fenotype{i} = range(i,1)+((range(i,2)-range(i,1))/(2^lengthS(i)-1))*fenotype{i};
    end
    
    Fitness = U2(fenotype{1},fenotype{2});
    if sum(Fitness < 0)>0
        Fitness = Fitness-min(Fitness);
    end
    
    
    
    
    if method == 1
        % Roulette Wheel "Natural Selection"
        Fitness2 = 1./(1+Fitness); %% only for minimization
        Fitness2 = Fitness2/sum(Fitness2);
        [~,ind] = max(Fitness2);
        
        %%% The best value of generation
        cum_p=cumsum(Fitness2);
        dup = rand([size_popu,1])';
        Rep = cum_p'<dup(:);
        Rep = sum(Rep,2)+1;
        parent_set1 = genotype{1}(Rep,:);
        parent_set2 = genotype{2}(Rep,:);
        
        
    else
        % Tournament "Natural Selection"
        [~,ind] = min(Fitness);
        dup = randi([1 size_popu],[size_popu,1])';
        parent_set1=zeros(size_popu,lengthS(1));
        parent_set2=zeros(size_popu,lengthS(2));
        for i=1:size_popu
             %%% avoid incest
            j=0;
            if dup(i)==i
                while j==0
                    dup(i) = randi([1 size_popu],[1,1])';
                    if dup(i)~=i
                        j=1;
                    end
                end
            end
            [~,Sel]=min([Fitness(i),Fitness(dup(i))],[],2);
            if Sel==1
                parent_set1(i,:) = genotype{1}(i,:);
                parent_set2(i,:) = genotype{2}(i,:); 
            else
                parent_set1(i,:) = genotype{1}(dup(i),:);
                parent_set2(i,:) = genotype{2}(dup(i),:);
            end
        end

    end
    
     X_best_curr=[fenotype{1}(ind),fenotype{2}(ind)];
     Fx_best_curr= U2(fenotype{1}(ind),fenotype{2}(ind));
    if generation==1
         fprintf('f(x*)= 0 at x*=[0,0], best fitness found: f(x*)= %8.6f x=[ %8.6f, %8.6f]  \n',Fx_best_curr,X_best_curr(1),X_best_curr(2))
         X_best_prev=X_best_curr;
         Fx_best_prev=Fx_best_curr;
    else
        if Fx_best_curr<Fx_best_prev
             fprintf('f(x*)= 0 at x*=[0,0], best fitness found: f(x*)= %8.6f x=[ %8.6f, %8.6f]  \n',Fx_best_curr,X_best_curr(1),X_best_curr(2))
              X_best_prev=X_best_curr;
              Fx_best_prev=Fx_best_curr;
        else
             fprintf('f(x*)= 0 at x*=[0,0], best fitness found: f(x*)= %8.6f x=[ %8.6f, %8.6f]  \n',Fx_best_prev,X_best_prev(1),X_best_prev(2))
        end
        
    end
        
        
     
    
    % CROSSOVER
    dup = randi([1 size_popu],[size_popu,1])';
    dup2 = randi([1 size_popu],[size_popu,1])';
    pointC=randi([2,min(lengthS)-1],[1,size_popu]);
    %     for k = 1:2:size_popu
    %
    parentvar1M = parent_set1(dup,:);
    parentvar2M = parent_set2(dup,:);
    parentvar1F = parent_set1(dup2,:);
    parentvar2F = parent_set2(dup2,:);
    for i=1:2:size_popu
        genotype{1}(i,1:pointC(i))= parentvar1M(i,1:pointC(i));
        genotype{1}(i,pointC(i)+1:lengthS(1))= parentvar1F(i,pointC(i)+1:lengthS(1));
        genotype{1}(i+1,1:pointC(i))= parentvar1F(i,1:pointC(i));
        genotype{1}(i+1,pointC(i)+1:lengthS(1))= parentvar1M(i,pointC(i)+1:lengthS(1));
        
        genotype{2}(i,1:pointC(i))= parentvar2M(i,1:pointC(i));
        genotype{2}(i,pointC(i)+1:lengthS(2))= parentvar2F(i,pointC(i)+1:lengthS(2));
        genotype{2}(i,1:pointC(i))= parentvar2F(i,1:pointC(i));
        genotype{2}(i,pointC(i)+1:lengthS(2))= parentvar2M(i,pointC(i)+1:lengthS(2));
        
    end
    mutacion= (random('unif',0,1,[size_popu,lengthS(1),2])<mut_prob);
    genotype{1}=abs(mutacion(:,:,1)-genotype{1});
    genotype{2}=abs(mutacion(:,:,2)-genotype{2});

end


[X1,X2]=meshgrid(-5.12:0.1:5.12,-5.12:0.1:5.12);
f=U2(X1,X2);
figure('color',[1 1 1]);
hold off
surfc(X1,X2,f),
shading interp
hold on
plot3(X_best_curr(1),X_best_curr(2),Fx_best_curr,'ro')
