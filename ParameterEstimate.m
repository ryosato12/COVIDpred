%% 
clear;
close all;
rng('default');

%% 
load('TotalCases.mat')

%%
totalPopulation=1936259; 
tspan=20:length(TotalCases);  

%% Find parameter values 
noStartPoints=10;
lowerBound=[1E-8  0.0004 1E-2  1E-6];      
upperBound=[1E-4  0.010  1     1   ];     
xstart=.5*(lowerBound+upperBound);                            

problem=createOptimProblem('fmincon','objective', ...
                              @(z)SEIR_RUN_ODE45(z, ...
                              totalPopulation, ...
                              TotalCases, ...
                              tspan), ...
                              'x0',xstart, ...
                              'lb',lowerBound, ...
                              'ub',upperBound);

problem.options=optimoptions(problem.options, ...
                               'MaxFunEvals',9999, ...
                               'MaxIter',9999);
numstartpoints=noStartPoints;                              
ms=MultiStart('UseParallel',true,'Display','iter');      
[b,fval,exitflag,output,manymins]=run(ms,problem,numstartpoints);  

for i=1:length(manymins)
    SEIRParameters(i,:)=manymins(i).X;       
end

for i=1:length(manymins)
    fvalues(i)=manymins(i).Fval;           
end

for i=1:length(manymins)
    ExitFlags(i)=manymins(i).Exitflag;     
end

%% 
I0=TotalCases(20)-TotalCases(10);            
E0=.40*(TotalCases(20)-TotalCases(10));          
R0=.95*TotalCases(10);
CI0=TotalCases(20);

%%
S0=totalPopulation-I0-R0-E0;
initialValues=[S0,I0,E0,R0,CI0];

%%
[t,y]=ode45(@(t,y) SEIR_Model(t,y,SEIRParameters(1,:)),tspan,initialValues);

S=y(:,1);
I=y(:,2);
R=y(:,3);
E=y(:,4); 
CI=y(:,5);

%%
% figure(1); clf; 
% plot(tspan,CI./max(abs(CI)))
% hold on
% % scatter(tspan,TotalCases(20:end)./max(TotalCases(20:end)),'filled')
% plot(tspan,E./max(abs(E)))
% plot(tspan,S./max(abs(S)))
% plot(tspan,I./max(abs(I)))
% plot(tspan,R./max(abs(R)))
% xlabel('Days')
% ylabel('Normalized Daily Incidence')
% legend('Culmulative','Exposed','Susceptible','Infectious','Recovered')
% hold off 

%%
% figure(2); clf; 
% plot(tspan(2:end),diff(CI))
% hold on
% scatter(tspan(2:end),diff(TotalCases(20:end)),'filled')
% xlabel('Days')
% ylabel('Daily Incidence')
% xlim([20 140])
% title('Santa Clara Positive Cases')

%%
figure(3); clf; 
plot(tspan,CI)
hold on
scatter(tspan,TotalCases(20:end),'filled')
xlabel('Days')
ylabel('Daily Incidence')
xlim([20 140])
title('Santa Clara Positive Cases')


%%
save('ParameterEstimate.mat','SEIRParameters')

%% Compute the error  
function value=SEIR_RUN_ODE45(z,totalPopulation,TotalCases,tspan) 
    I0=TotalCases(20)-TotalCases(10);            
    E0=.40*(TotalCases(20)-TotalCases(10));          
    R0=.95*TotalCases(10);
    CI0=TotalCases(20);

    S0=totalPopulation-I0-R0-E0;
    
    initialValues=[S0,I0,E0,R0,CI0];
    
    [t,y]=ode45(@(t,y) SEIR_Model(t,y,z),tspan,initialValues);        
    CI=y(:,5);
    diff=CI-TotalCases(20:end);
    value=norm(diff,2);
end 

%%  
function dydt=SEIR_Model(t,y,z)
    beta=z(1);     
    k=1/4;
    h=1/10; 
    delta=z(4);
    sigma=z(2);
    epsilon=z(3); 
    
    % Initiate DE variables
    dydt=zeros(5,1);
    
    S=y(1);
    E=y(2);
    I=y(3);
    R=y(4);
    
    dydt(1)=-beta*S*I-sigma*E;
    dydt(2)=beta*S*I-k*E-h*E;
    dydt(3)=k*E-delta*I-epsilon*I;
    dydt(4)=delta*I+h*E-sigma*R;
    dydt(5)=beta*S*I;
end


