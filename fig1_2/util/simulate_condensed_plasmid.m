function [counts,yend] = simulate_condensed_plasmid(par)

%This function simulates the condensed plasmid system (equivalent to the
%full multiplasmid system where all plasmids are identical)

if ~isfield(par,'additive_fitness')
    par.additive_fitness = 1;
end

if ~isfield(par,'rand_init')
    par.rand_init = 0;
end

%Generate full matrices for conjugation and loss
[N,K,types,type_combo_index] = generate_plasmid_matrix(par.n_plasmid);

%Condense full matrices
[cN,cK,cond_combos,num_variants] = condense_plasmid_matrix(N,K,types,type_combo_index,par.n_plasmid);

%Vector of fitness losses
pop_plasmids = (0:par.n_plasmid)';
if par.additive_fitness
    Fv = (1 - pop_plasmids.*par.f);
else
    Fv = (1 - par.f).^pop_plasmids;
end

cK = cK.*repmat(Fv',size(cK,1),1);

disp('Finished constructing matrices!')

%Get death vector
Dv = zeros(size(Fv)) + par.delta1;

%Define wrapped rhs functions
rhs_fun = @(y) RHS_condensed_plasmid(y,par.gammac,par.rl,cN,cK,Dv,Fv,par.alpha1,par.Gamma,cond_combos,num_variants);

%Simulation parameters
tn = 1; %Timestep tracking variable
u = 0; %Loop exit variable

%Initial conditions and storage vectors
y0 = [zeros(par.n_plasmid + 1,1) + 1/(par.n_plasmid+1) ;1];

%Randomize ICs if option selected
if par.rand_init
    y0 = rand(size(y0));
    disp(y0)
end

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
    
    %relchange = abs(dy1./yi);
    %relchange = max(relchange(~isinf(relchange)));
    %dt = par.maxchange/relchange; %Set time step
    
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
    if min(y) < 0 || sum(isnan(y)) > 0 %Check for negative concentration
        disp('NEGATIVE CONCENTRATIONS OR NAN DETECTED. SIMULATION ENDED. REDUCE STEP SIZE');
        u = 1;
    end
    
    %error = max(abs((ystore(:,tn) - ystore(:,tn-1))./ystore(:,tn)));
    error = norm(dy1);
    if error < par.error_threshold %Check for steady-state
        u = 1; %Loop exit variable
    end
    
end

disp('Finished simulating!')

%Remove zeros and extract final value
ystore = ystore(:,1:tn);
yend = ystore(:,tn);

counts = yend(1:end-1);

end

