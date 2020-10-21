% Crosstalk model parameters
% Author: Liz Fedak
% Created: 5/23/19
% Updated: 10/20/20


% Timestep and numerics control for ODE solver
% Timescale is minutes

rA       = 6; % Factor by which histone phosphorylation speeds binding rate of all proteins
kNHEJ    = log(2)/20/rA; % NHEJ, DSBs have a half-life of 20 minutes
kD       = 0.0015; % Detection rate of transcribed-strand CPDs. kDT > kDN
kSNHEJ   = log(2)/300/rA; % NHEJ requiring preprocessing. kSNHEJ < kNHEJ
kPL      = 1e-8; % RPA-bound ssDNA breakage rate, low
kMRN     = log(2)/200/rA; % MRN complex binding rate
MRNs     = 2; % Factor by which MRN binding rate increases during S phase
kBRCA    = log(2)/60/rA; % Rate of BRCA-mediated Rad51 outcompeting RPA
kSSA     = log(2)/200/rA; % Rate of single-strand annealing DSB repair, slower than kBRCA
kHR      = log(2)/600/rA; % HR rate, slow
kNER     = log(2)/20/rA; % NER rate, half-life of 20 min
kATM     = 0.5; % Maximal ATM activation rate
jATM     = 4000; % ATM half-maximal activation coefficient
XPAs     = 10; % increase in ATR binding efficacy during S phase
kATR     = 0.03; % Maximal ATR activation rate
jATR     = 40000; % ATR half-maximal activation cefficient
kdATM    = 0.005; % ATM basal dephosphorylation rate
kdATR    = 0.005; % ATR basal dephosphorylation rate

r        = 50*60; % rate at which DNAP traverses the DNA, bp/sec converted to minutes
kpd      = 0.1; % Pol delta binding rate to DNA
kdpd     = 0.005; % Basal Pol delta dissociation rate from DNA
kdpda    = 1e-5; % ATR-mediated blocking of Pol delta binding
kTLS     = log(2)/100/rA; % Rate of TLS
kaa      = 0.008; % Rate of ATM phosphorylation by ATR
ktop     = 4; % Enhanced binding rate of ATR to ATM-bound end-resected DSBs

Pold_tot  = 3.3e3; % total Pol delta, estimated from earlier model.
G_tot    = 3.3e9; % number of base pairs in the genome, fixed
ATM_tot  = 5000; % Total cellular ATM concentration, assumed constant
ATR_tot  = 1000; % Total cellular ATR concentration, assumed constant

par = [rA, kNHEJ, kD, kSNHEJ, kPL, kMRN, MRNs, kBRCA, kSSA, kHR, ...
       kNER, kATM, jATM, XPAs, kATR, jATR, kdATM, kdATR, ...
       r, kpd, kdpd, kdpda, kTLS, kaa, ktop, Pold_tot, G_tot, ATM_tot, ATR_tot];
   