function [plasmid_counts,yend,types,ystore] = simulate_multi_plasmid(par)
%This function simulates the n-plasmid ODE model

if ~isfield(par,'additive_fitness')
    par.additive_fitness = 1;
end

%Generate matrices for conjugation and loss 
[N,K,types,type_combo_index] = generate_plasmid_matrix(par.n_plasmid);

%Vector of fitness losses
if par.additive_fitness
    Fv = (1 - sum(types,2).*par.f);
else
    Fv = (1 - par.f).^sum(types,2);
end

K = K.*repmat(Fv',size(K,1),1);

%Get death vector
Dv = zeros(size(Fv)) + par.delta1;

%Define wrapped rhs functions
rhs_fun = @(y) RHS_multi_plasmid(y,par.gammac,par.rl,N,K,Dv,Fv,par.alpha1,par.Gamma,type_combo_index);

%Simulation parameters
tn = 1; %Timestep tracking variable 
u = 0; %Loop exit variable

%Initial conditions and storage vectors
y0 = [zeros(size(Fv)) + 1/length(Fv);1];
ystore = zeros(length(y0),1e6);
ystore(:,1) = y0;
tstore = zeros(1,1e6); %Time tracking variable
dt = par.dt;

%Simulate until error tolerance is met
while u == 0
    
    tn = tn + 1; %Advance timestep number
    
    %COMPUTE K1
    yi = ystore(:,tn-1);
    dy1 = rhs_fun(yi);
    
%     relchange = max(abs(dy1./yi));
%     dt = par.maxchange/relchange; %Set time step
    
    %COMPUTE K2
    yii = yi + 0.5*dt.*dy1;
    dy2 = rhs_fun(yii);
    
    %COMPUTE K3
    yiii = yi + 0.5*dt*dy2;
    dy3 = rhs_fun(yiii);
    
    %COMPUTE K4
    yiv = yi + dt*dy3;
    dy4 = rhs_fun(yiv);
    
    %TRUE TIMESTEP
    y = yi + dt*((1/6)*dy1 + (1/3)*dy2 + (1/3)*dy3 + (1/6)*dy4);
    
    %Set small values to zero to speed up integrator
    y(y<1e-20) = 0;
    
    %STORAGE
    ystore(:,tn) = y;
    tstore(tn) = tstore(tn-1) + dt;

    %EXIT CHECKS
    if min(y) < 0 %Check for negative concentration
        disp('NEGATIVE CONCENTRATIONS. SIMULATION ENDED. REDUCE STEP SIZE');
        u = 1;
    end
    
    error = max(abs((ystore(:,tn) - ystore(:,tn-1))./ystore(:,tn)));
    if error < par.error_threshold %Check for steady-state
        u = 1; %Loop exit variable    
    end
    
end

%Remove zeros and extract final value
ystore = ystore(:,1:tn);
yend = ystore(:,tn);

plasmid_counts = zeros(par.n_plasmid+1,1);
type_sums = sum(types,2);
for i = 0:par.n_plasmid
    plasmid_counts(i+1,:) = sum(yend(type_sums == i,:),1);
end

end

