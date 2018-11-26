%% Cargo transport, get diff rates (attach and detach)


clear;  clc; format shortG;
close all;

% rng(1); % random number generator

kin.spr = 0.3; % kinesin spring constant, Vale/Coppin 1996 and Nedelec 2002.

dyn.spr = 0.065; % Oiwa/Sakakibara 2005, dyn.spr = 0.01 from Gennerich NatCom 2015
dyn.type = 1; %1=kunwar maxVel=-212nm/s, 2=yildiz with cargo, maxVel=-513, 3 == yildiz DDB, no cargo
dyn.catch = 0; % Dynein's catch bond

F=0; % Extenal force
runs = 5000; % number of iterations
% ----------------------------------------------------

% intialize vectors
A = zeros(1,6); % A = [kin.num,dyn.num, Fd, Fs, mean(cargovel), SD]
B = A; % B = [kin.num,dyn.num, Fd, Fs, mean(RL), SD]
D = zeros(runs,1); % vector of all run length
E = D; % vector of all run time
allvel = E; % vector that contains <vel>
G = struct; % structure that contains {mot attach, cargo position, RT, mot type, mot force, cargo Pos}
meanvel = 0; % vector that contain mean velocity of allvel
emp = struct; % structure that contains empirical kinetics
fd = 1;fs = 0;dd = 0; % counters

% - - - - - - - - - - -- - - -
sav = 0; % save file?
sendtext=0; % send text once Simulation Completed?
category = 1;
traceplot = 0;
% if traceplot ~= 0
%     category = 3; % traces only
% end
saveTraceYN = 0;
%{
1=emp, 2=diff, 3=trace_for_emp(regular tug of war, 4=never detach motors,
5 = cv vs (Nkin or Ndyn)


kfor = [64/3, 318/9, 50+318/9*0.25, 75+318/9*0.25, 100+318/9*0.25 125+318/9*0.25];
% 100a, 212b, 400c, 600d, 800e, 1000f nm/s

dyn.kback = 318/9*0.25; %1/s, 212 nm/s @ 0pN || 8*(kfor-kback) = 212nm
dyn.kdet = [0.05b, 0.1,0.2,0.5,1];



%}
A = struct;
kin.num = 1;

dyn.num = 1;
dyn.catch = 0;

% %kfor = [64/3, 318/9, 50+318/9*0.25, 75+318/9*0.25, 100+318/9*0.25 125+318/9*0.25] ;%64/3; % fixed, 100 nm/s
% kfor = [400 600 800 1000]/8+318/9*0.25;
% %kfor=1000/8+318/9*0.25;
% kdet =[0.05, 0.1, 0.2, 0.5, 1];% fixed
% %kdet =0.05;
% dyn.kback = 318/9*0.25; % fixed
% velList=100:100:1000;
% kdetList=[0.05 0.1 0.2 0.5 1 2 5];

% Fds=1:8;
% Fss=1:8;
velList=400; kdetList=0.5;Fds=1;Fss=1;
kintypes=[1 2];
dynspList=[0.3 0.1 0.065 0.03 0.01 0.003 0.001]; %default0.065
kattList=[1 2 5 10 20 50]; %default 5
kbackLists=[0, 1, 2, 4, 8.8, 16]; %default 8.8
meanvelcargo=zeros(length(velList),length(kdetList),length(Fds),length(Fss),length(kintypes),length(dynspList),length(kattList),length(kbackLists)); %Fd,Fs,kintype
varvelcargo=meanvelcargo;
%tic;
%ratiobacktofor=0.25;
%kforList=velList/8/(1-ratiobacktofor);
for ii = 1:length(kdetList)
    dyn.kdet = kdetList(ii);
    for jj = 1:length(velList)
        %dyn.kfor = kforList(jj);
        %dyn.kback=ratiobacktofor*dyn.kfor;
        for type = kintypes
            kin.type = type;
            for Fd = Fds
                dyn.Fd = Fd;
                for Fs = Fss
                    dyn.Fs = Fs;
                    %name = sprintf('kfor%g_kdet%gFs%gFd%gtype%g\n',dyn.kfor, dyn.kdet, dyn.Fs, dyn.Fd, kin.type);
                    %fprintf(name);
                    %h = waitbar(0,name);
                    emp = struct;
                    for isp=1:length(dynspList)
                        dyn.spr=dynspList(isp);
                            for ikatt=1:length(kattList)
                                dyn.kon=kattList(ikatt);
                                kin.kon=kattList(ikatt);
                                    for ikback=1:length(kbackLists)
                                        dyn.kback=kbackLists(ikback);
                                        dyn.kfor=velList(jj)/8+dyn.kback;
                    
 %%%%%%%%%%%%%%%%%%%%MC run                   
                    parfor irun = 1:runs
                        %tic
                        Simu = tOhashi34_varykforvel_mod(F,dyn,kin,category);
                        %toc
%                         if category == 1 || 4
%                             G = getsimpleG2(G,Simu,ii); % G(z2,ii) = G(dyn.num,runs)
%                             emp = getEmpRate3(emp,ii,G);
%                         end
%                         if traceplot == 1
%                             trace(Simu,name,kin,z)
%                         elseif traceplot == 2
%                             SimuCounter = 1;
%                             trace2(Simu,name,kin,saveTraceYN,ii,runs,SimuCounter,G);
%                         end
                        
                        D(irun) = Simu(end-1).cX - Simu(1).cX;  % Get the final position of the cargo
                        E(irun) = Simu(end).t;
                        %h = waitbar(irun/runs);
                        %                             close(h);
                        
                    end
                    
                    allvel(:) = D(:)./E(:);  % mean velocity of all obtained Simu data.
                    meanvel = mean(allvel(:));
%                     A(type).vel(Fs,Fd) = meanvel;
%                     A(type).kn(Fs,Fd) = kin.num;
%                     A(type).dn(Fs,Fd) = dyn.num;
%                     A(type).std(Fs,Fd) = std(allvel(:),E(:));
%                     A(type).sem(Fs,Fd) = std(allvel(:),E(:)) / sqrt(runs);
                    meanvelcargo(jj,ii,Fd,Fs,kin.type,isp,ikatt,ikback)=meanvel;
                    varvelcargo(jj,ii,Fd,Fs,kin.type,isp,ikatt,ikback)=std(allvel(:),E(:));
                    %close(h)
                    %%%%%%%%%%%%%%%%%%%%%%%%end of MC run
                                    end
                            end
                    end
                    
                end
            end
        end
    
    
    %fprintf('      %1.1fmin %2.1fs',floor(toc/60),rem(toc,60))
%     if dyn.kdet == 0.05
%         strkdet = '0p05';
%     elseif dyn.kdet == 0.1
%         strkdet = '0p1';
%     elseif dyn.kdet == 0.2
%         strkdet  = '0p2';
%     elseif dyn.kdet == 0.5
%         strkdet = '0p5';
%     elseif dyn.kdet == 1
%         strkdet = '1p0';
%     else
%         error('none');
%     end
%     str = sprintf('vel%g_kdet%s',8*(dyn.kfor - dyn.kback),strkdet);
%     save(str);
    end
  save('allvelandparam')  
end






