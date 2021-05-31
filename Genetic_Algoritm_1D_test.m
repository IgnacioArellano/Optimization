function Genetic_Algoritm_1D_test
U = @(T) 204165.5/(330-2*T) + 10400/(T-20);
%% mínimo U = 1225.2   T= 55.0835
lon=17;
pobla=15;
Timin =40;
Timax =90;
Pmut=0.02;
flag_mut=0;
Iterations=3000;

x=zeros(pobla,lon);
for i=1:pobla
x(i,:) = randi([0 1],[lon,1])';
end

fprintf('----------------------------------\n')
%%%%%%%%% Población inicial
% x= [1,1,0,1,1,0,0,1,1,0,0,0,1,0,1; 1,0,0,0,0,1,0,1,0,1,1,1,0,1,0; 0,0,0,1,1,0,1,0,1,1,1,0,1,0,1;...
%     1,0,0,0,0,0,1,1,0,0,1,1,1,0,1;0,0,0,0,1,1,1,0,0,1,0,0,1,1,1; 1,0,0,1,0,0,1,0,1,0,1,1,0,0,0;...
%     0,1,0,1,1,0,1,0,0,1,1,0,0,0,1; 1,0,0,1,1,0,1,0,1,1,1,0,0,1,1; 1,1,1,1,0,0,1,0,0,0,1,1,0,1,0;...
%     0,0,1,1,1,1,1,0,0,1,1,1,0,0,1];
for itera=1:Iterations
    
    %%%%%%%%% fitness
    f=zeros(pobla,1);
    for i=1:pobla
        DV=bi2de(x(i,:),'left-msb');
        Ti= Timin+ ((Timax-Timin)*DV)/(2^lon-1);
        f(i)=U(Ti);
    end
    
    %%%%%%%%%%%%%%%%%%
    fprintf('Minimum = %8.3f  the best fitness: %8.6f muta= %d \n',1225.2 ,min(f),flag_mut)
    flag_mut=0;
    
    %%%%%%%%%%%%%%%%%%%%%%%%% Ruleta selección natural
    if min(f)<0
        f=f-min(f);
    end
    
    F=1./(1+f);
    ProF= F/sum(F);
    
    
%     %%%% Ruleta selección
     CumP=cumsum(ProF)';
     %duplicamos la población
     dup = rand([pobla,1])';
     Rep = CumP<dup(:);
     Rep = sum(Rep,2)+1;
     Rpob = f(Rep);
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%% Selección por Torneo natural
    dup = randi([1 pobla],[pobla,1])';
   % dup=[8,4,2,9,10,7,8,1,4,2];
    
    [~,Tp] = min([f,f(dup)],[],2);
    progs(:,:,1)=x;
    progs(:,:,2)=x(dup,:);
    nprogs=zeros(size(x));
    for i=1:pobla
        nprogs(i,:)=progs(i,:,Tp(i));
    end
    %%%%%%%%%%%%%%%%%%%%%%%%% Cruce
    pointC=randi([2,lon-1],[1,pobla]);
    ite=1;
    while ite<=pobla
        t=randi([1,pobla],[1,2]);
        if t(1)~=t(2)
            x(ite,1:pointC(ite))=nprogs(t(1),1:pointC(ite));
            x(ite,pointC(ite)+1:lon)=nprogs(t(2),pointC(ite)+1:lon);
            ite=ite+1;
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%% Mutación
    mut =find( random('unif',0,1,[1,pobla])<Pmut);
    
    if ~isempty(mut)
        pointm=randi([1,lon],[1,length(mut)]);
        x(mut,pointm)=~x(mut,pointm);
        flag_mut=1;
    end
    
    
end















