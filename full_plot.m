% Main plotting function
% Author: Liz Fedak
% Created: 10/13/19
% Updated: 10/20/20

function full_plot(varargin)
% FULL_PLOT() Generates all plots in "What p53 sees: ATM and ATR activation
% through crosstalk between DNA damage response pathways."
% Takes two vectors as optional input.
% It interprets the first vector as a list of figures to generate as
% follows:
% 1: Fig. 4
% 2: Fig. 5
% 3: Fig. 6
% 4: Fig. S7
% 5: Fig. S8
% 6: Fig. S9
% 7: Fig. S10
% This code does not include figures that rely on data sets from other
% papers.
%
% It interprets the second vector as a list of experimental conditions to
% use as follows:
% 0: Base model, no changes
% 1: ATM kinase activity does not activate ATR
% 2: ATR kinase activity does not activate ATM
% 3: ATM and ATR do not phosphorylate H2AX
% 4: Replication stress on MDS activates ATR
% 5: ssDNA does not break to form DSBs
% 6: ATM dissociates from end-resected DSBs
% 7: ATR dissociates from extensively end-resected DSBs
%
% If given no options, it will generate all figures for the base model.
    
    fig_list = 1:8; % if no input given, generate all figures
    exp_list = 0; % if no input given, use only base model
    
    if ~isempty(varargin) % but override figures & experimental conditions if given.
        fig_list = varargin{1};
        if length(varargin) > 1 % just in case someone wants to input 3 vectors
            exp_list = varargin{2};
        end
    end
    

    for i=exp_list
        
        [ode, ic_fun, par, DSBi, SSi, ATMi, ATRi] = choose_model(i); % Select model according to numerical value
        
        if ismember(1,fig_list) % If 1 is selected
           fig4(ode, ic_fun, par, DSBi, SSi, ATMi, ATRi) % Generate plots in Fig. 4
        end
        
        if ismember(2,fig_list) % If 2 is selected
           fig5(ode, ic_fun, par, DSBi, SSi, ATMi, ATRi) % Generate plots in Fig. 5
        end
        
        if ismember(3,fig_list) % If 3 is selected
           fig6(ode, ic_fun, par, ATMi, ATRi) % Generate plots in Fig. 6
        end
        
        if ismember(4,fig_list) % If 4 is selected
           figS7(ode, ic_fun, par, ATMi, ATRi) % Generate plots in Fig. S7
        end
        
        if ismember(5,fig_list) % If 5 is selected
           figS8(ode, ic_fun, par, ATRi) % Generate plots in Fig. S8
        end
        
        if ismember(6,fig_list) % If 6 is selected
           figS9(ode, ic_fun, par, ATMi, ATRi) % Generate plots in Fig. S9
        end
        
        if ismember(7,fig_list) % If 7 is selected
           figS10(ode, ic_fun, par, ATMi, ATRi) % Generate plots in Fig. S10
        end
        
    end
  
end

%% CHOOSING THE MODEL

function [ode, ic_fun, par, DSBi, SSi, ATMi, ATRi] = choose_model(n)
% Returns the ODE function, initial condition function, and parameter set for
% the base model or one of the seven experimental conditions. Can only
% handle one number at a time. The numbers refer to:
% 0: Base model, no changes
% 1: ATM kinase activity does not activate ATR
% 2: ATR kinase activity does not activate ATM
% 3: ATM and ATR do not phosphorylate H2AX
% 4: Replication stress on MDS activates ATR
% 5: ssDNA does not break to form DSBs
% 6: ATM dissociates from end-resected DSBs
% 7: ATR dissociates from extensively end-resected DSBs

params % 'par' is an array that is already defined here

% We need to know the indices of specific parameters within 'par'

rA_ind = 1;
kPL_ind = 5;
kaa_ind = 24;

switch n
    case 0
        ode = @full_ODEs;
        ic_fun = @s_init;
        DSBi  = 3:10; % Indices for classes of DSB
        SSi   = [11:15 28]; % Indices for classes of ssDNA. Final equation is directly solvable.
        ATMi  = 16:21; % Indices for classes of ATM
        ATRi  = 22:27; % Indices for classes of ATR
    case 1
        ode = @exp1_ODEs; % NOT FINISHED, NEED ODE
        ic_fun = @(s, D_IR, D_UV) [s_init(s, D_IR, D_UV) 0 0 0 0]; % four additional ODEs each with initial value 0
        DSBi  = [3:10 28:31]; % Indices for classes of DSB
        SSi   = [11:15 32]; % Indices for classes of ssDNA. Final equation is directly solvable.
        ATMi  = 16:21; % Indices for classes of ATM
        ATRi  = 22:27; % Indices for classes of ATR
    case 2
        ode = @full_ODEs;
        ic_fun = @s_init;
        par(kaa_ind) = 0;
        DSBi  = 3:10; % Indices for classes of DSB
        SSi   = [11:15 28]; % Indices for classes of ssDNA. Final equation is directly solvable.
        ATMi  = 16:21; % Indices for classes of ATM
        ATRi  = 22:27; % Indices for classes of ATR
    case 3
        ode = @full_ODEs;
        ic_fun = @s_init;
        par(rA_ind) = 1;
        DSBi  = 3:10; % Indices for classes of DSB
        SSi   = [11:15 28]; % Indices for classes of ssDNA. Final equation is directly solvable.
        ATMi  = 16:21; % Indices for classes of ATM
        ATRi  = 22:27; % Indices for classes of ATR
    case 4
        ode = @exp4_ODEs;
        MDS_per_Gy = 400; % Estimate of MDS created per Gy of gamma radiation from Ward 1988
        ic_fun = @(s, D_IR, D_UV) [s_init(s, D_IR, D_UV) MDS_per_Gy*D_IR*(1-s) MDS_per_Gy*D_IR*2*s 0 0]; % four additional ODEs; two are MDS classes split into replicated/non-replicated components
        % two additional parameters
        kCDSB = 1e-5;
        kMR   = log(2)/200/rA;
        par   = [par kCDSB kMR];
        DSBi  = 3:10; % Indices for classes of DSB
        SSi   = [11:15 30:32]; % Indices for classes of ssDNA. Final equation is directly solvable.
        ATMi  = 16:21; % Indices for classes of ATM
        ATRi  = 22:27; % Indices for classes of ATR
    case 5
        ode = @full_ODEs;
        ic_fun = @s_init;
        par(kPL_ind) = 0;
        DSBi  = 3:10; % Indices for classes of DSB
        SSi   = [11:15 28]; % Indices for classes of ssDNA. Final equation is directly solvable.
        ATMi  = 16:21; % Indices for classes of ATM
        ATRi  = 22:27; % Indices for classes of ATR
    case 6
        ode = @exp6_ODEs;
        ic_fun = @s_init;
        DSBi  = 3:10; % Indices for classes of DSB
        SSi   = [11:15 28]; % Indices for classes of ssDNA. Final equation is directly solvable.
        ATMi  = 16:21; % Indices for classes of ATM
        ATRi  = 22:27; % Indices for classes of ATR
    case 7
        ode = @exp7_ODEs;
        ic_fun = @s_init;
        DSBi  = 3:10; % Indices for classes of DSB
        SSi   = [11:15 28]; % Indices for classes of ssDNA. Final equation is directly solvable.
        ATMi  = 16:21; % Indices for classes of ATM
        ATRi  = 22:27; % Indices for classes of ATR
    otherwise
        disp('Invalid number entered for experimental condition.')
end
end

%% FIGURES

function fig4(ode, ic_fun, par, DSBi, SSi, ATMi, ATRi)
% Plots used in Figure 4
    
    % tiny bit of hardcoding here
    tmax = 48; % run for 48 hours
    
    % FIG. 4A: GAMMA RADIATION
    s = 0; % Percentage of S phase complete
    D_IR = 10; % IR damage induced, in Gy
    D_UV = 0; % UV damage induced, in J/m^2
    init = ic_fun(s, D_IR, D_UV); % generate initial conditions for above three parameters
    par_basic = [s_par(s, D_IR, D_UV) par]; % giving ICs for solvable state variables
    % init(1) = 0; % Uncomment to see activity outside of S phase
    dam = '10 Gy IR'; % a string specifying magnitude and units of damage
    
    % Show how all state variables evolve for given completed percentage of
    % S phase and damage state.
    basic_plot(ode, init, par_basic, dam, DSBi, SSi, ATMi, ATRi, tmax)
    
    % FIG. 4B: UV RADIATION
    
    D_IR = 0; % IR damage induced, in Gy
    D_UV = 10; % UV damage induced, in J/m^2
    init = ic_fun(s, D_IR, D_UV); % generate initial conditions for above three parameters
    par_basic = [s_par(s, D_IR, D_UV) par]; % giving ICs for solvable state variables
    dam = '10 J/m^2 UV'; % a string specifying magnitude and units of damage
     
    basic_plot(ode, init, par_basic, dam, DSBi, SSi, ATMi, ATRi, tmax)
end

% ----- %

function fig5(ode, ic_fun, par, DSBi, SSi, ATMi, ATRi)
% Plots used in Figure 5
    
    % Cell cycle effects: Plot ATM or ATR induction after IR or UV damage,
    % depending on what point in the cell cycle damage was induced.
    
    %FIG. 5A: DSB repair
     plot_cellcycle(ode, ic_fun, par, DSBi, 10, 0, ...
       'DSB repair rate relative to the start of S phase');
   
    % FIG. 5B: Photoproduct repair
     plot_cellcycle(ode, ic_fun, par, SSi, 0, 10, ...
       'Photoproduct repair rate relative to the start of S phase');
    
    % FIG. 5C: ATM after gamma radiation
     plot_cellcycle(ode, ic_fun, par, ATMi, 10, 0, ...
        'Level of ATM activation with changes in IR induction relative to the start of S phase');
    
    % FIG. 5D: ATR after gamma radiation
     plot_cellcycle(ode, ic_fun, par, ATRi, 10, 0, ...
         'Level of ATR activation with changes in IR induction relative to the start of S phase');
    
    % FIG. 5E: ATM after UV radiation
     plot_cellcycle(ode, ic_fun, par, ATMi, 0, 10, ...
       'Level of ATM activation with changes in UV induction relative to the start of S phase');
    
    % FIG. 5F: ATR after UV radiation
     plot_cellcycle(ode, ic_fun, par, ATRi, 0, 10, ...
         'ATR, 10 J/m^2 UV induction');
end

% ----- %

function fig6(ode, ic_fun, par, ATMi, ATRi)
% Plots used in Figure 6

    % Dose response curves; area under curve for 20 hours
    
    % FIG. 6A: Gamma radiation response
     dose_response(ode, ic_fun, par, 24, [20 0], 'Gy', ATMi, ATRi)
     
    % FIG. 6B: UV radiation response
     dose_response(ode, ic_fun, par, 24, [0 20], 'J/m^2', ATMi, ATRi)
     
    ATR_tot = 29; % position of ATR_tot
    par(ATR_tot) = 100; % decrease total ATR in cell to 100
    
    % FIG. 6C: Gamma radiation response, low ATR
     dose_response(ode, ic_fun, par, 24, [20 0], 'Gy', ATMi, ATRi)
     
    % FIG. 6D: UV radiation response, low ATR
     dose_response(ode, ic_fun, par, 24, [0 20], 'J/m^2', ATMi, ATRi)
end

% ----- %

function figS7(ode, ic_fun, par, ATMi, ATRi)
% Plots used in Figure S7

     % FIG. 7A: 10 Gy gamma radiation
     atmatr(ode, ic_fun, par, 10, 0, ATMi, ATRi, '10 Gy \gamma radiation')
     
     % FIG. 7B: 10 J/m^2 UV
     atmatr(ode, ic_fun, par, 0, 10, ATMi, ATRi, '10 J/m^2 UV')
     
     % FIG. 7C: 0.5 Gy gamma radiation
     atmatr(ode, ic_fun, par, 0.5, 0, ATMi, ATRi, '0.5 Gy \gamma radiation')
     
     % FIG. 7D: 1 J/m^2 UV
     atmatr(ode, ic_fun, par, 0, 1, ATMi, ATRi, '1 J/m^2 UV')

end

% ----- %

function figS8(ode, ic_fun, par, ATRi)
% Plots used in Figure S8
    kMRN = 7; % index of kMRN; check with ODE file
    kATR = 16;

    % FIG. 8A: End resection K/O in response to gamma radiation, like in Kousholt 2012
     knockout(ode, ic_fun, par, 10, 0, kMRN, ATRi, 'ATR', 0)
     
    % FIG. 8B: ATR kinase K/O in response to gamma radiation, like in Gamper et al. 2013
     knockout(ode, ic_fun, par, 10, 0, kATR, ATRi, 'ATR', 0.1)
     
end

% ----- %

function figS9(ode, ic_fun, par, ATMi, ATRi)
% Plots used in Figure S9
     ATM_tot = 29; % index of ATM_tot, check in ODE file
     ATR_tot = 30;

    % FIG. 9A: ATM response to gamma radiation, ATR K/O
     knockout(ode, ic_fun, par, 10, 0, ATR_tot, ATMi, 'ATM', 0) % gamma, ATR K/O
     
    % FIG. 9B: ATR response to gamma radiation, ATM K/O
     knockout(ode, ic_fun, par, 10, 0, ATM_tot, ATRi, 'ATR', 0) % gamma, ATM K/O
     
    % FIG. 9C: ATM response to UV radiation, ATR K/O
     knockout(ode, ic_fun, par, 0, 10, ATR_tot, ATMi, 'ATM',0) % UV, ATR K/O
     
    % FIG. 9D: ATR response to UV radiation, ATM K/O
     knockout(ode, ic_fun, par, 0, 10, ATM_tot, ATRi, 'ATR', 0) % UV, ATM K/O

end

% ----- %

function figS10(ode, ic_fun, par, ATMi, ATRi)
% Plots used in Figure S10

    kD = 4; % index of kD, same as ODE file

    % NER knockout in response to UV, like in Ray et al. 2016
     knockout(ode, ic_fun, par, 0, 20, kD, ATRi, 'ATR', 0)
     knockout(ode, ic_fun, par, 0, 20, kD, ATMi, 'ATM', 0)

end

%% HELPER FUNCTIONS %%

function basic_plot(ode, init, par, dam, DSBi, PPi, ATMi, ATRi, tmax)
% Basic plot showing how each quantity evolves.
% ode, init & par: same as above. 'dam' is a string w/ the name of the
% damage; used for the title. D is amount of damage in Gy for IR and J/m^2 
% for UV. Next 4 variables give the indices of all variables that belong 
% to the indicated compartment.

    colors = hsv(4); % make things look nice
     
    % run ODE
    [t,y] = ode45(ode, [0 tmax*60], init, [], par);
    
    PP = solvable_DEs(t,par);
    
    y = [y, PP]; % Append solvable DEs onto the end of the solution
    
    % Plot fits vs. data.
    figure
    gcf;
     
    set(gcf,'DefaultAxesFontSize',17) % make font big everywhere
    set(gcf, 'Position', [400, 400, 400, 400])
    lw = 2.5; % Change line width.
     
    % Plot DSBs, PPs, ATM, and ATR.
        plot(t/60, sum(y(:,DSBi),2), 'col', colors(1,:), 'LineWidth', lw)
        hold on
        plot(t/60, sum(y(:,PPi),2), 'col', colors(2,:), 'LineWidth', lw)
        hold on
        plot(t/60, sum(y(:,ATMi),2), 'col', colors(3,:), 'LineWidth', lw)
        hold on
        plot(t/60, sum(y(:,ATRi),2), 'col', colors(4,:), 'LineWidth', lw)
        title(['Model response to ' dam ' radiation'])
        legend('DSBs', 'PPs', 'ATMp', 'ATRp')
        ylim([0 5000])
        xlim([0 tmax])
        xlabel('Time (h)')
        ylabel('Proteins (ND)')

end

% ----- DOSE RESPONSE CURVES ----- %

function dose_response(ode, ic_fun, par, tmax, Dmax, units, ATMi, ATRi)
% Plot showing ATM and ATR dose response curves as the amount of UV or IR
% damage is varied (does both at once). The quantity being tracked is the
% area under the curve, since this is what most strongly affects the p53
% response.
% Dmax = [D_IR D_UV] specifies the maximum amount of damage to be induced, tmax the 
% maximum time in minutes. The parameter dam is again a string specifying 
% whether the damage induced is IR or UV. Really only going to work for one
% type of damage at a time, dimensionally. This is meant to handle cases
% where we expose the cell to only one form of damage.
    
    D = 1/sum(Dmax):0.25/sum(Dmax):1; % Initialize damage range, must start at 1/max(Dmax) for scaling
    n = length(D);
    ATM = zeros(1,n); % Initialize level of ATM/ATR activation
    ATR = zeros(1,n);
    % AT_tot = par(end); % total ATM or ATR is always last parameter
    % max_activation = AT_tot*tmax; % The maximum possible activation level; everything is always on
    s = 0;

    
    for i=1:n
        
        % run ODEs = i/s_tot; % percentage of S phase complete
        init_temp = ic_fun(s, D(i)*Dmax(1), D(i)*Dmax(2)); % s=0, Dmax(1) = D_IR, Dmax(2) = D_UV
        par_temp = [s_par(s, D(i)*Dmax(1), D(i)*Dmax(2)) par];
        [t,y] = ode45(ode, [0 tmax*60], init_temp, [], par_temp);
        ATM_tot = sum(y(:,ATMi),2);
        ATR_tot = sum(y(:,ATRi),2);
        ATM(i) = sum(diff(t).*ATM_tot(2:end)); % approximation of area under curve
        ATR(i) = sum(diff(t).*ATR_tot(2:end));
         
    end
    % Initialize plot and increase font size for legibility
    figure
     gcf;
     
     set(gcf,'DefaultAxesFontSize',17) % make font big everywhere
     set(gcf, 'Position', [400, 400, 400, 400])
     lw = 2.5; % Change line width.
     Atot = par(end-1) + par(end); % total number of ATM and ATR proteins

     
     plot(D*sum(Dmax), ATM, 'r', 'LineWidth', lw)
     hold on
     plot(D*sum(Dmax), ATR, 'b', 'LineWidth', lw)
     hold on
     plot(D*sum(Dmax), ATM+ ATR, 'm', 'LineWidth', lw)

     % title([dam ' Dose Response Curves'])
     legend('ATM', 'ATR', 'ATM + ATR')
     %ylim([1 sum(Atot)*tmax*60])
     xlabel(['Amount of Damage (' units ')'])
     ylabel('Area under curve')

end

% ----- ATM AND ATR K/Os ----- %

function knockout(ode, ic_fun, par, D_IR, D_UV, par_ind, ATi, AT_str, scale)
% Generates plots of ATM and ATR activity in and ouside of S phase for ATR
% or ATM K/O, respectively. Here, ATtot is the index of the total
% concentration of {ATM, ATR} and ATi gives the  of {ATR, ATM}. AT_str is a
% string that contains the name of the kinase whose activity is displayed,
% and scale is a scaling factor in case the knockout doesn't completely
% abrogate an interaction.
    
    % colors = hsv(4); % make things look nice
    
    s = 0;
    init = ic_fun(s, D_IR, D_UV);
    par1 = [s_par(s, D_IR, D_UV) par];
    tmax = 24;
    
    % WT, S phase
    [t1,y1] = ode45(ode, [0 tmax*60], init, [], par1);
    PP = solvable_DEs(t1,par1);
    y1 = [y1, PP]; % Append solvable DEs onto the end of the solution
    
    init(1) = 0; % Uncomment to see activity outside of S phase
    
    % WT, outside S phase
    [t2,y2] = ode45(ode, [0 tmax*60], init, [], par1);
    PP = solvable_DEs(t2,par1);
    y2 = [y2, PP];
    
    init = ic_fun(s, D_IR, D_UV);
    par2 = [s_par(s, D_IR, D_UV) par];
    par2(par_ind) = scale*par2(par_ind); % zeros(1,length(par_ind));
    
    % K/O, S phase
    [t3,y3] = ode45(ode, [0 tmax*60], init, [], par2);
    PP = solvable_DEs(t3,par2);
    y3 = [y3, PP];
    
    init(1) = 0;
    
    % K/O, outside S phase
    
    [t4,y4] = ode45(ode, [0 tmax*60], init, [], par2);
    PP = solvable_DEs(t4,par2);
    y4 = [y4, PP];
    
    
    % Plot fits vs. data.
    figure
     gcf;
     
     set(gcf,'DefaultAxesFontSize',17) % make font big everywhere
     set(gcf, 'Position', [400, 400, 400, 400])
     lw = 2.5; % Change line width.
     
     % Plot DSBs, PPs, ATM, and ATR.
        plot(t1/60, sum(y1(:,ATi),2), 'b', 'LineWidth', lw)
        hold on
        plot(t2/60, sum(y2(:,ATi),2), 'r', 'LineWidth', lw)
        hold on
        plot(t3/60, sum(y3(:,ATi),2), 'b--', 'LineWidth', lw)
        hold on
        plot(t4/60, sum(y4(:,ATi),2), 'r--', 'LineWidth', lw)
        legend('WT, S phase', 'WT, G1 phase', 'K/O, S phase', 'K/O, G1 phase')% 'k_{top} = 1, k_{MDSs} = 0.1; S phase', 'k_{top} = 1, k_{MDSs} = 0.1; G1 phase')
        % ylim([0 5000])
        % xlim([0 100])
        xlabel('Time (h)')
        ylabel([AT_str ' (Number of Proteins)'])

end

function atmatr(ode, ic_fun, par, D_IR, D_UV, ATMi, ATRi, dam)
% Generates plots of ATM and ATR activity in and ouside of S phase for ATR
% or ATM K/O, respectively. Here, ATtot is the index of the total
% concentration of {ATM, ATR} and ATi gives the  of {ATR, ATM}.
    
    s = 0;
    init = ic_fun(s, D_IR, D_UV);
    par1 = [s_par(s, D_IR, D_UV) par];
    
    % S phase
    [t1,y1] = ode45(ode, [0 24*60], init, [], par1);
    PP = solvable_DEs(t1,par1);
    y1 = [y1, PP]; % Append solvable DEs onto the end of the solution
    
    init(1) = 0; % Activity outside of S phase
    
    % Outside S phase
    [t2,y2] = ode45(ode, [0 24*60], init, [], par1);
    PP = solvable_DEs(t2,par1);
    y2 = [y2, PP];
    
    
    % Plot model.
    figure
     gcf;
     
     set(gcf,'DefaultAxesFontSize',17) % make font big everywhere
     set(gcf, 'Position', [400, 400, 400, 400])
     lw = 2.5; % Change line width.
     
     % Plot DSBs, PPs, ATM, and ATR.
        plot(t1/60, sum(y1(:,ATMi),2), 'b', 'LineWidth', lw)
        hold on
        plot(t1/60, sum(y1(:,ATRi),2), 'r', 'LineWidth', lw)
        hold on
        plot(t2/60, sum(y2(:,ATMi),2), 'b--', 'LineWidth', lw)
        hold on
        plot(t2/60, sum(y2(:,ATRi),2), 'r--', 'LineWidth', lw)
        title(['ATM and ATR response to ' dam ' radiation'])
        legend('ATM, S phase', 'ATR, S phase', 'ATM, G1 phase', 'ATR, G1 phase')
        ylim([0 5000])
        xlim([0 24])
        title(dam)
        xlabel('Time (h)')
        ylabel('Number of Proteins')

end

% ----- CELL CYCLE ----- %

function plot_cellcycle(ode, ic_fun, par, j, D_IR, D_UV, plot_title)
% Plots response to 10 J/m^2 UV damage induced in different parts of the
% cell cycle. Func is an ODE function handle, parfunc is either uv_par or
% ir_par, initfunc is either uv_init or ir_init, par is the par vector, j
% contains indices of the state variable to be plotted, plot_title is a string 
% to be used for the plot title. D_IR and D_UV are the amount of IR and UV
% damage induced in Gy and J/m^2, respectively.

opts = []; % odeset('MaxStep', 1e-1); % it likes jumping around the peak

% Change cell cycle duration parameters here. 
total_cycle = 22;
s_tot = 6; % number of hours cell spends in S phase, no-damage case
colors = hsv(s_tot); % rainbows
not_s_tot = total_cycle - s_tot; % number of hours cell spends outside S phase

figure

% Plot behavior of cells in S phase.
for i=0:s_tot-1 % do it in hours for now
    s = i/s_tot; % percentage of S phase complete
    par_temp = [s_par(s, D_IR, D_UV) par];
    init_temp = ic_fun(s, D_IR, D_UV);
    [t,y] = ode45(ode, [0 total_cycle*60], init_temp, opts, par_temp);
    PP = solvable_DEs(t,par); % Include solvable DEs
    y = [y, PP];
    p = plot(t/60, sum(y(:,j),2),'r','LineWidth', 3,'col', colors(i+1,:));
    hold on
end

% Plot behavior of cells where damage is induced before S phase.
for i=1:not_s_tot
    s = 0; % Outside of S phase, treat s as 0 and set G0 = 0.
    par_temp = [s_par(s, D_IR, D_UV) par];
    init_temp = ic_fun(s, D_IR, D_UV);
    init_temp(1) = 0; % Set G0 = 0 so replication doesn't initiate
    n = length(init_temp); % Needed to trim solvable quantities from ODE input
    % Pre-replication stage lasts for i hours
    [t1,y1] = ode45(ode, [0 i*60], init_temp, opts, par_temp);
    PP = solvable_DEs(t1,par); % Include solvable DEs
    y1 = [y1, PP];
    % When S phase starts, the initial condition is where the system is
    % after i hours.
    init_temp = y1(end,1:n);
    init_temp(1) = 3.3e9; % length of genome, set G0 to nonzero so that replication initiates
    [t2,y2] = ode45(ode, [i*60 total_cycle*60], init_temp, opts, par_temp);
    PP = solvable_DEs(t2,par); % Include solvable DEs
    y2 = [y2, PP];
    if i==1
        q = plot([t1; t2]/60,[sum(y1(:,j),2); sum(y2(:,j),2)],'k--','LineWidth', 3);
    elseif i==not_s_tot
        q = plot([t1; t2]/60,[sum(y1(:,j),2); sum(y2(:,j),2)],'k-.','LineWidth', 3);
    else
        q = plot([t1; t2]/60,[sum(y1(:,j),2); sum(y2(:,j),2)],'k','LineWidth', 3, 'HandleVisibility','off');
    end
    hold on
end

hold off

    title(plot_title)

    xlabel('Time after damage exposure (h)')
    ylabel('Net protein activation/damage level')
    legend('Damage induced at start of S phase', '1 hour into S phase', '2 hours into S phase', ...
        '3 hours into S phase', '4 hours into S phase', '5 hours into S phase', ...
        '1 hour before S phase', '16 hours before S phase')
    xlim([0 total_cycle])
    set(gca,'fontsize', 16);
end


%% Plot helper functions

% ----- %

function init = s_init(s, D_IR, D_UV)
% If damage is induced in the middle of S phase, correctly distribute the
% initial conditions according to how much DNA has already been replicated.
% D_IR or D_UV are only used the first time damage is induced.

% Consistent with full_ODEs
LDSB0      = D_IR*36*(1+s);
LDSBA0     = 0;
LCDSB0     = D_IR*4*(1+s); % total DNA = G_0(1 + s)
LCDSBA0    = 0;
LER0       = 0;
LERA0      = 0;
LEER0      = 0;
LEERA0     = 0;
LP0        = D_UV*4500*(1-s); % total untranscribed DNA = G_0(1 - s)
LPD0       = 0;
LPDA0      = 0;
FP0        = 0;
FPA0       = 0;
ATMpDSB0   = 0;
ATMpCDSB0  = 0;
ATMpLER0   = 0;
ATMpLEER0  = 0;
ATMpP0     = 0;
ATMpF0     = 0;
ATRpDSB0   = 0;
ATRpCDSB0  = 0;
ATRpLER0   = 0;
ATRpLEER0  = 0;
ATRpP0     = 0;
ATRpF0     = 0;

G0       = 3.3e9*(1-s); % 3.3e9 bps in genome
Poldf0   = 3.3e3; % if it's outside of S phase, all Pol-delta is free, but it doesn't actually matter

init = [G0, Poldf0, LDSB0, LDSBA0, LCDSB0, LCDSBA0, LER0, LERA0, LEER0, LEERA0, ...
    LP0, LPD0, LPDA0, FP0, FPA0, ...
    ATMpDSB0, ATMpCDSB0, ATMpLER0, ATMpLEER0, ATMpP0, ATMpF0, ...
    ATRpDSB0, ATRpCDSB0, ATRpLER0, ATRpLEER0, ATRpP0, ATRpF0];
end

% ----- %

function par = s_par(s, D_IR, D_UV)
% If damage is induced in the middle of S phase, correctly distribute the
% initial conditions according to how much DNA has already been replicated.
% D_IR or D_UV are only used the first time damage is induced. If S phase
% starts later, these four parameters are overwritten.

p_t      = 0.01; % hidden parameter representing percentage of actively transcribed genome

% LMDS0      = D_IR*400*(1+s);
PTR0     = D_UV*4500*p_t*2*s;

par  = PTR0;
end

% ----- %

function PP = solvable_DEs(t,par)
% Generates directly solvable quantities

    kDN = par(5);
    kDT = par(6);
    PNR = par(1)*exp(-kDN*t); % directly solvable
    PTR = par(2)*exp(-kDT*t);
    PP = PTR; % [PNR, PTR];
end

