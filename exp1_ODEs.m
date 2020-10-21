% ODEs for experimental condition 1
% Liz Fedak
% Created: 10/13/19
% Updated: 10/20/20


function ydot = exp1_ODEs(t,y,par)


%% Extract population values

G         = y(1); % number of unreplicated base pairs in genome
Poldf     = y(2); % free Pol delta available for replication
LDSB      = y(3); % Simple double-strand breaks (DSBs)
LDSBA     = y(4); % ATM- or ATR-modified DSBs
LCDSB     = y(5); % Complex DSBs (C-DSBs)
LCDSBA    = y(6); % ATM- or ATR-modified C-DSBs
LER       = y(7); % End-resected DSBs (ER)
LERA      = y(8); % ATM- or ATR-modified end-resected DSBs
LEER      = y(9); % Extensively end-resected DSBs (EER)
LEERA     = y(10); % ATM- or ATR-modified extensively end-resected DSBs
LP        = y(11); % Undetected photoproducts on DNA that has not already been replicated
LPD       = y(12); % Detected photoproducts
LPDA      = y(13); % ATR-modified detected photoproducts
FP        = y(14); % Stalled replication forks
FPA       = y(15); % ATR-bound ssDNA on stalled replication forks
ATMpDSB   = y(16); % Active ATM bound to simple DSBs
ATMpCDSB  = y(17); % Active ATM bound to complex DSBs
ATMpER    = y(18); % Active ATM bound to end-resected DSBs
ATMpEER   = y(19); % Active ATM bound to extensively end-resected DSBs
ATMpP     = y(20); % Active ATM bound to ATR on photoproducts
ATMpF     = y(21); % Active ATM bound to ATR on stalled replication forks
ATRpDSB   = y(22); % Active ATR bound to simple DSBs
ATRpCDSB  = y(23); % Active ATR bound to complex DSBs
ATRpER    = y(24); % Active ATR bound to end-resected DSBs
ATRpEER   = y(25); % Active ATR bound to extensively end-resected DSBs
ATRpP     = y(26); % Active ATR bound to ATR on photoproducts
ATRpF     = y(27); % Active ATR bound to ATR on stalled replication forks

% NEW ODES FOR EXPERIMENTAL CASE

% Classes of DSB that are bound to ATM but not ATR

LDSBATM   = y(28);
LCDSBATM  = y(29);
LERATM    = y(30);
LEERATM   = y(31);


%% Extract parameters 

% One initial condition for an easily solvable term

LPR0      = par(1);

rA        = par(2);
kNHEJ     = par(3);
kD        = par(4);
kSNHEJ    = par(5);
kPL       = par(6);
kMRN      = par(7);
MRNs      = par(8);
kBRCA     = par(9);
kSSA      = par(10);
kHR       = par(11);
kNER      = par(12);
kATM      = par(13);
jATM      = par(14);
XPAs      = par(15);
kATR      = par(16);
jATR      = par(17);
kdATM     = par(18);
kdATR     = par(19);

r         = par(20);
kpd       = par(21);
kdpd      = par(22);
kdpda     = par(23);
kTLS      = par(24);
kaa       = par(25);
ktop      = par(26);

Pold_tot  = par(27);
Gtot      = par(28);
ATM_tot   = par(29);
ATR_tot   = par(30);


%% S phase specific and conserved quantities 

Fa    = Pold_tot - Poldf; % conservation equation, active replication forks

if G < 1
    G = 0;
    dGdt = 0;
    rho = 0; % there should be no rate of replication if there's no DNA to replicate
    MRNs = 1; % and MRN complex is not upregulated
    XPAs = 1; % ATR binding is not upregulated
else
    rho = r*Fa./G; % 1/G is used to determine density of lesions
    dGdt = -r*Fa - rho.*LP;
end

LPR  = LPR0*exp(-kD*t);

ATMp  = ATMpDSB + ATMpCDSB + ATMpER + ATMpEER + ATMpP + ATMpF;
ATRp  = ATRpDSB + ATRpCDSB + ATRpER + ATRpEER + ATRpP + ATRpF;
ATM   = ATM_tot - ATMp;
ATR   = ATR_tot - ATRp;
ATR_g = kATR*XPAs*ATR./(jATR + FP + FPA + LPD + LPDA + LDSB + LDSBA + LDSBATM + LCDSB + LCDSBA + LCDSBATM + LER + LERA + LERATM + LEER + LEERA + LEERATM);
ATM_f = kATM*ATM./(jATM + LDSB + LDSBA + LDSBATM + LCDSB + LCDSBA + LCDSBATM + LER + LERA + LERATM + LEER + LEERA + LEERATM + ATRpP + ATRpF);


%% Core DEs

% Pol delta
dPoldfdt     = -kpd*Poldf.*(G/Gtot) + (kdpd + kdpda*ATRp)*Fa;

% DSBs and other IR lesions
dLDSBdt      = -(kNHEJ + kMRN*MRNs)*LDSB - LDSB.*(ATM_f + ATR_g);
dLDSBAdt     = ATR_g.*(LDSB + rA*LDSBATM) - rA*(kNHEJ + kMRN*MRNs)*LDSBA; % CHANGED
dLCDSBdt     = -(kSNHEJ + kMRN*MRNs)*LCDSB - LCDSB.*(ATM_f + ATR_g);
dLCDSBAdt    = ATR_g.*(LCDSB + rA*LCDSBATM) -rA*(kSNHEJ + kMRN*MRNs)*LCDSBA; % CHANGED
dLERdt       = kMRN*MRNs*(LDSB + LCDSB) + kPL*(LPD + FP) - LER.*(ATM_f + ATR_g) - (kBRCA + kSSA)*LER;
dLERAdt      = ATR_g.*(LER + rA*LERATM) + kPL*(LPDA + FPA) + rA*(kMRN*MRNs*(LCDSBA + LDSBA) - (kBRCA + kSSA)*LERA); % CHANGED
dLEERdt      = kBRCA*LER - LEER.*(ATM_f + ATR_g) - kHR*LEER;
dLEERAdt     = ATR_g.*(LEER + rA*LEERATM) + rA*(kBRCA*LERA - kHR*LEERA); % CHANGED

% UV photoproducts
dLPdt       = -kD*LP - rho.*LP;
dLPDdt       = kD*(LP + LPR) + kTLS*(FP + rA*FPA) - (kNER + kPL)*LPD;
dLPDAdt      = LPD.*ATR_g - (rA*kNER + kPL)*LPDA;

% S-phase-only DEs
dFPdt        = rho.*LP - FP.*ATR_g - (kTLS + kPL)*FP;
dFPAdt       = FP.*ATR_g - (rA*kTLS + kPL)*FPA;

% ATMp (FIRST FOUR CHANGED)
dATMpDSBdt   = (LDSB + rA*(LDSBA + LDSBATM)).*ATM_f - kdATM*ATMpDSB - rA*(kNHEJ + kMRN*MRNs).*ATMpDSB; 
dATMpCDSBdt  = (LCDSB + rA*(LCDSBA + LCDSBATM)).*ATM_f - kdATM*ATMpCDSB - rA*(kSNHEJ + kMRN*MRNs)*ATMpCDSB;
dATMpERdt    = (LER + rA*(LERA + LERATM)).*ATM_f + kPL*(ATMpP + ATMpF) - kdATM*ATMpER + rA*kMRN*MRNs*(ATMpDSB + ATMpCDSB) - rA*(kSSA + kBRCA).*ATMpER;
dATMpEERdt   = (LEER + rA*(LEERA + LEERATM)).*ATM_f - kdATM*ATMpEER + rA*(-kHR + kBRCA)*ATMpER;
dATMpPdt     = kaa*ATRpP.*ATM_f - (rA*kNER + kPL + kdATR + kdATM)*ATMpP;
dATMpFdt     = kaa*ATRpF.*ATM_f - (rA*kTLS + kPL + kdATR + kdATM)*ATMpF;

% ATRp (FIRST FOUR CHANGED)
dATRpDSBdt  = (LDSB + rA*LDSBATM + ktop*rA*LDSBA).*ATR_g - (kdATR + rA*kNHEJ + rA*kMRN*MRNs)*ATRpDSB;
dATRpCDSBdt = (LCDSB + rA*LCDSBATM + ktop*rA*LCDSBA).*ATR_g - (kdATR + rA*kSNHEJ + rA*kMRN*MRNs)*ATRpCDSB;
dATRpERdt   = (LER + rA*LERATM + ktop*rA*LERA).*ATR_g + rA*kMRN*MRNs*(ATRpDSB + ATRpCDSB) + kPL*(ATRpP + ATRpF) - (rA*(kSSA + kBRCA) + kdATR)*ATRpER;
dATRpEERdt  = (LEER + rA*LEERATM + ktop*rA*LEERA).*ATR_g + rA*kBRCA*ATRpER - (rA*kHR + kdATR)*ATRpEER;
dATRpPdt    = (LPD + ktop*rA*LPDA).*ATR_g - (rA*kNER + kPL + kdATR)*ATRpP;
dATRpFdt    = (FP + ktop*rA*FPA).*ATR_g - (rA*kTLS + kPL + kdATR)*ATRpF;

% new ODEs
dLDSBATMdt  = ATM_f.*LDSB - rA*ATR_g.*LDSBATM - rA*(kNHEJ + kMRN*MRNs)*LDSBATM;
dLCDSBATMdt = ATM_f.*LCDSB - rA*ATR_g.*LCDSBATM - rA*(kSNHEJ + kMRN*MRNs)*LCDSBATM;
dLERATMdt   = ATM_f.*LER - rA*ATR_g.*LERATM + rA*(kMRN*MRNs*(LCDSBATM + LDSBATM) - (kBRCA + kSSA)*LERATM);
dLEERATMdt  = ATM_f.*LEER - rA*ATR_g.*LEERATM + rA*(kBRCA*LERATM - kHR*LEERATM);
    
ydot = [dGdt;
        dPoldfdt;
        dLDSBdt;
        dLDSBAdt;
        dLCDSBdt;
        dLCDSBAdt;
        dLERdt;
        dLERAdt;
        dLEERdt;
        dLEERAdt;
        dLPdt;
        dLPDdt;
        dLPDAdt;
        dFPdt;
        dFPAdt;
        dATMpDSBdt;
        dATMpCDSBdt;
        dATMpERdt;
        dATMpEERdt;
        dATMpPdt;
        dATMpFdt;
        dATRpDSBdt;
        dATRpCDSBdt;
        dATRpERdt;
        dATRpEERdt;
        dATRpPdt;
        dATRpFdt;
        dLDSBATMdt;
        dLCDSBATMdt;
        dLERATMdt;
        dLEERATMdt];

