// EA_GNSS10 model

// Reference: Gerali, A., S. Neri, L. Sessa, F. Signoretti. (2010). 
// Credit and Banking in a DSGE Model of the Euro Area. Journal of Money Credit and Banking
// Supplement to Vol. 42, No. 6.

// The original model code has been provided by Stefano Neri.

// Last edited by: Sebastian Schmidt, January 2012
// Note: In order to define an output gap, I added a flexible price, flexible wage, 
// flexible loan and deposit rate case of the economy by eliminating adjustment costs
// in the respective sectors. The model block is added below the original model block.


// Authors' comments:

// model with: TWO WAGES; INVESTMENT ADJ. COSTS;  VARIABLE CAPITAL UTILIZATION;  CONSUMPTION HABITS
//            STICKY BANK RATES, PRICES & WAGES ? la Rotemberg with indexation to both past and st.st. inflation
// 9 blocks: 1) PATIENTs  2) IMPATIENTs  3) CAPITAL PRODUCERS  4) ENTREPRENEURS   5) BANKS   6) RETAILERS        
//           7) LABOR MKT WITH ONE UNION FOR EACH LABOR TYPE  8) AGGREGATION & EQUILIBRIUM   9) MONETARY POLICY  

//    AGGIUNTA DEL 16 APRILE: PROBLEMI INTERTEMPORALI TUTTI IN TERMINI REALE

// THIS VERSION: April 20th, 2009
// BBS has all items dated (t) - Banking capital in real terms is defined as: K_B(t)/p(t)
// Banking profits are defined in the model code at time (t) - All equations in exp form
// This is the full version (monop. competitive banking sector & BANK capital ? la Gerali with sticky rates) 



var 
c_p       // 1  PATIENT   HHs
h_p       // 2  PATIENT   HHs
d_p       // 3  PATIENT   HHs
l_p       // 4  PATIENT   HHs
lam_p     // 5  PATIENT   HHs
J_R       // 6  PATIENT   HHs
j_B       // 7  PATIENT   HHs
pie_wp    // 8  PATIENT   HHs  
c_i       // 9  IMPATIENT HHs
h_i       // 10 IMPATIENT HHs
b_i       // 11 IMPATIENT HHs 
l_i       // 12 IMPATIENT HHs
lam_i     // 13 IMPATIENT HHs
s_i       // 14 IMPATIENT HHs
pie_wi    // 15 IMPATIENT HHs
I         // 16 CAPITAL PRODUCERS
q_k       // 17 CAPITAL PRODUCERS
c_e       // 18 ENTREPRENEURS
k_e       // 19 ENTREPRENEURS
l_pd      // 20 ENTREPRENEURS
l_id      // 21 ENTREPRENEURS
b_ee      // 22 ENTREPRENEURS
y_e       // 23 ENTREPRENEURS
lam_e     // 24 ENTREPRENEURS
s_e       // 25 ENTREPRENEURS
u         // 26 ENTREPRENEURS Capital utilization rate
d_b       // 27 BANKS
b_h       // 28 BANKS
b_e       // 29 BANKS
r_d       // 30 BANKS
r_bh      // 31 BANKS
r_be      // 32 BANKS
R_b       // 33 BANKS
K_b       // 34 BANKS
pie       // 35 RETAILERS
x         // 36 RETAILERS
C         // 37 AGGREGATION & EQUILIBRIUM
Y         // 38 AGGREGATION & EQUILIBRIUM
D         // 39 AGGREGATION & EQUILIBRIUM
BE        // 40 AGGREGATION & EQUILIBRIUM
BH        // 41 AGGREGATION & EQUILIBRIUM
B         // 42 AGGREGATION & EQUILIBRIUM
w_p       // 43 AGGREGATION & EQUILIBRIUM
w_i       // 44 AGGREGATION & EQUILIBRIUM
J_B       // 45 AGGREGATION & EQUILIBRIUM
q_h       // 46 AGGREGATION & EQUILIBRIUM
K         // 47 AGGREGATION & EQUILIBRIUM
PIW       // 48 AGGREGATION & EQUILIBRIUM
r_ib      // 49 MONETARY POLICY
r_k       // 50 CAPITAL RENTAL RATE
ee_z      // 51 EXOGENOUS PROCESSES
A_e       // 52 EXOGENOUS PROCESSES
ee_j      // 53 EXOGENOUS PROCESSES
mk_d      // 54 EXOGENOUS PROCESSES
mk_be     // 55 EXOGENOUS PROCESSES
mk_bh     // 56 EXOGENOUS PROCESSES
ee_qk     // 57 EXOGENOUS PROCESSES
m_i       // 58 EXOGENOUS PROCESSES (IMPATIENT LTV)
m_e       // 59 EXOGENOUS PROCESSES (ENTREPRENEURS LTV)
eps_y     // 60 EXOGENOUS PROCESSES
eps_l     // 61 EXOGENOUS PROCESSES
eps_K_b   // 62 EXOGENOUS PROCESSES
Y1        // 63 output a prezzi di ss
rr_e      // 64 Entrep. Real Rate
aux1      // 65 auxiliary variable
bm        // 66 banks intermediation margins
spr_b     // 67 average bank spread (active-passive)
h
I_h
ee_qh
vi
//**************************************************************************
// Replication Variables                                                  //*
        interestPol interestH interestF inflation loansH loansF output consumption investment deposits interestDep bankcapital;
//**************************************************************************


varexo e_A_e e_eps_K_b e_j e_l e_me e_mi e_mk_be e_mk_bh e_mk_d e_r_ib e_qk e_y e_z e_qh e_vi;  //13 varexo



parameters  
            beta_p j phi beta_i m_i_ss beta_e m_e_ss alpha eksi_1 eksi_2    // HOUSEHOLDS & ENTREPRENEURS
            a_i a_p a_e gamma_p gamma_i gamma_e    ni                     // HOUSEHOLDS & ENTREPRENEURS
            eps_l_ss kappa_w                                                // HOUSEHOLDS (labor params)
            eps_d eps_bh eps_be                                             // BANKS 
            mk_d_ss mk_bh_ss mk_be_ss r_be_ss  r_bh_ss r_k_ss               // BANKS (SS)
            gamma_b beta_b delta_kb kappa_kb                             // BANKS 
            eps_y_ss kappa_p ind_p ind_w                                    // RETAILERS
            kappa_i kappa_d kappa_be kappa_bh deltak                        // OTHERS
            ind_d ind_be ind_bh                                             // OTHERS
            rho_ib phi_pie phi_y                                            // POLICY
            piss  r_ib_ss                                                   // STEADY STATE
            rho_ee_z rho_A_e rho_ee_j rho_mi rho_me rho_eps_y               // SHOCKS
            rho_mk_d rho_mk_be rho_mk_bh rho_ee_qk rho_eps_l rho_eps_K_b    // SHOCKS
            rho_ee_qh deltah kappa_ih
		rho_vi vi_ss;
            

% *********************			
% CALIBRATED PARAMETERS
% *********************


beta_p       = 0.9963;                                                     % discount factor patient households
beta_i       = 0.975;                                                      % discount factor impatient households     
beta_b       = beta_p;                                                     % discount factor bankers (not used in this version of the model)
beta_e       = beta_i;                                                     % discount factor entrepreneurs
j            = 0.2;                                                        % weight of housing in utility function
phi          = 7.822;                                                        % inverse Frisch elasticity of labor supply
m_i_ss       = 0.71  ;                                                      % loan-to-value ratio impatient households
m_e_ss       = 0.35 ;                                                      % loan-to-value ratio entrepreneurs
alpha        = 0.55;                                                      % capital share in the production function
eps_d        = -1.362400844;                                                   % elast. of subst. of deposits 
eps_bh       = 6.785728165;
eps_be       = 1.653421891; 
 
mk_d_ss      = eps_d   / (eps_d  - 1) ;                                    % steady state markdown on D (ok if eps_d<0; if eps_d>0 it should be eps_d/(eps_d+1) )
mk_bh_ss     = eps_bh  / (eps_bh - 1) ;                                    % steady state markup on loans to I
mk_be_ss     = eps_be  / (eps_be - 1) ;                                    % steady state markup on loans to E
book_ss      = 0; %-35                                                       % steady state value of (B-M)/D in the bank balance shhet
eps_y_ss     = 6;                                                          % 
eps_l_ss     = 5;                                                          % 
gamma_p      = 1;                                                          % shares of patient households
gamma_i      = 1; //1/3;                                                   % shares of impatient households
ni           = 0.9;                                                        % wage share of patient households
gamma_b      = 1; //0.10;												   % shares of bankers
gamma_e      = 1; //1 - gamma_p - gamma_i;								   % shares of entrepreneurs
deltak       = 0.049;                                                      % depreciation rate for physical capital
piss         = 1;                                                          % steady state gross inflation rate
r_ib_ss      = (piss/beta_p - 1) * (eps_d-1)/eps_d ;                       % steady state gross nominal interest rate 
r_be_ss      = r_ib_ss*eps_be/(eps_be-1) ;								   % steady state interest rate on loans to E
r_bh_ss      = r_ib_ss*eps_bh/(eps_bh-1) ;								   % steady state interest rate on loans to H
r_k_ss       = -(1-deltak)-m_e_ss*(1-deltak)*piss/beta_e*
               (1/(1+r_be_ss)-beta_e/piss)+1/beta_e;                       % steady state rental rate of capital
                                                         % fixed supply housing
eksi_1       = r_k_ss;                                                     % capital utilization cost parameter
eksi_2       = 0.1*r_k_ss;                                                 % capital utilization cost parameter
vi_ss=0.18;                                                                % Banking Capital ratio over Loans (BaselII)
eps_b        = mean([eps_bh,eps_be]);
% delta_kb     = r_ib_ss/vi * (eps_d - eps_b + vi*(eps_b-1))/((eps_b-1)*(eps_d-1));  % j_B with terms in r_ib                            
delta_kb     = r_ib_ss/vi_ss * (eps_d - eps_b + vi_ss*eps_d*(eps_b-1))/((eps_b-1)*(eps_d-1));                               

ind_d        = 0.0;                   % indexation deposit rates
ind_be       = 0.0;                   % indexation rates on loans to firms
ind_bh       = 0.0;                   % indexation rates on loans to households


deltah=0.01;

% *****************************************************************
% LOADING MEDIAN OF POSTERIOR: USES EXTRACT_MEDIAN_FROM_POSTERIOR.m (dummy way)
% *****************************************************************

load median_valuesROh.txt;
coeffs = median_valuesROh;

rho_ee_z	=	coeffs(1);  % 0.385953438168178	;
rho_ee_j     =	coeffs(2);  % 0.93816527333294	;
rho_me	=	coeffs(3);  % 0.921872719102206	;
rho_mi      =	coeffs(4);  % 0.90129485520182	;
rho_mk_d      =	coeffs(5);  % 0.922378382753078	;
rho_mk_bh	=	coeffs(6);  % 0.892731352899547	;
rho_mk_be	=	coeffs(7);  % 0.851229673864555	;
rho_A_e	=	coeffs(8);  % 0.873901213475799	;
rho_ee_qk	=	coeffs(9);  % 0.571692383714171	;
rho_eps_y	=	coeffs(10); % 0.294182239567384	;
rho_eps_l	=	coeffs(11); % 0.596186440884132	;
rho_eps_K_b	=	coeffs(12); % 0.813022758608552	;
rho_vi     =	coeffs(13); % 33.7705265016395	;
rho_ee_qh     =	coeffs(14); % 107.352040072465	;
kappa_p     =	coeffs(15); % 10.0305562248008	;
kappa_w     =	coeffs(16); % 2.77537377104213	;
kappa_i	=	coeffs(17); % 7.98005959044637	;
kappa_d	=	coeffs(18); % 9.04426718749482	;
kappa_be	=	coeffs(19); % 8.91481958034669	;
kappa_bh     =	coeffs(20); % 2.00384780180824  ;
kappa_kb      =	coeffs(21); % 0.750481873084311	;
kappa_ih       =	coeffs(22); % 0.303247771697294	;
phi_pie       =	coeffs(23); % 0.158112794106546	;
rho_ib       =	coeffs(24); % 0.300197804017489	;
phi_y	        =	coeffs(25); % 0.867003766306404	;
a_e         =   0.0     ;   % degree of habit formation: entrepreneurs
a_p         =   0.0     ;   % degree of habit formation: patient households
ind_p=coeffs(25);
ind_w=coeffs(26);
a_i=coeffs(27);     
//%------------------------------------------------------------
//% Model equations
//%------------------------------------------------------------

model;

//**************************************************************************
// Definition of Replication Variables in Terms of Original Model Variables //*

interestPol   = 400*exp(r_ib); 
interestH = 400*exp(r_bh);
interestF = 400*exp(r_be);                                         
inflation  = pie*100;    
loansH = BH*100;
loansF = BE*100;
output     = Y1*100; 
consumption = C*100;
investment = (I+I_h)*100;   
deposits = D*100;
interestDep = 400*(exp(r_d));
bankcapital = 100*K_b;                           
//**************************************************************************



// Model code:

////***********   1) PATIENT HHs ********************************************************6

  (1-a_i)*exp(ee_z)*(exp(c_p) - a_i*exp(c_p(-1)))^(-1) = exp(lam_p); // CON rescaling  (1) where a_p = a_e = a_i

j * (exp(e))  / exp(h_p) - exp(lam_p) * exp(q_h) + beta_p * exp(lam_p(+1)) * exp(q_h(+1))   = 0; // (2) 

exp(lam_p)  = beta_p * exp(lam_p(+1)) * (1+exp(r_d)) / exp(pie(+1)); // (3)


(1 - exp(eps_l)) * exp(l_p) + exp(l_p) ^(1+phi) / exp(w_p) * exp(eps_l)/exp(lam_p) 
                                         - kappa_w *( exp(pie_wp)     - exp(pie(-1)) ^ ind_w * piss ^ (1-ind_w) ) * exp(pie_wp)
   +  beta_p * exp(lam_p(+1))/exp(lam_p) * kappa_w *( exp(pie_wp(+1)) - exp(pie)     ^ ind_w * piss ^ (1-ind_w) ) * exp(pie_wp(+1)) ^2 / exp(pie(+1)) = 0 ; // Unions, labor supply (4)

exp(pie_wp) = exp(w_p) / exp(w_p(-1)) * exp(pie); // definition of wage inflation (5)

exp(c_p) + exp(q_h) * ( exp(h_p) - exp(h_p(-1)) ) + exp(d_p)  = exp(w_p) * exp(l_p)
   + (1+exp(r_d(-1)))*exp(d_p(-1))/exp(pie) + exp(J_R)/gamma_p ;  // patient household budget constraint (6), exp(J_R)/gamma_p = t_p


////***********   2) IMPATIENT HHs ********************************************************7

  (1-a_i)*exp(ee_z)*(exp(c_i) - a_i*exp(c_i(-1)))^(-1)  = exp(lam_i); // CON rescaling (7)

j * (exp(ee_j))  / exp(h_i) - exp(lam_i) * exp(q_h) + beta_i * exp(lam_i(+1)) * exp(q_h(+1))  + exp(s_i) * exp(m_i) *exp(q_h(+1))   * exp(pie(+1))  = 0;    // (8) s_i is multiplier on borrowing constraint
exp(lam_i) - beta_i * exp(lam_i(+1)) * (1+exp(r_bh)) / exp(pie(+1)) = exp(s_i) * (1+exp(r_bh)); // (9)

(1 - exp(eps_l)) * exp(l_i) + exp(l_i) ^(1+phi) / exp(w_i) * exp(eps_l)/exp(lam_i) 
                                         - kappa_w *( exp(pie_wi)     - exp(pie(-1))^ind_w * piss ^ (1-ind_w) ) * exp(pie_wi)
   +  beta_i * exp(lam_i(+1))/exp(lam_i) * kappa_w *( exp(pie_wi(+1)) - exp(pie)    ^ind_w * piss ^ (1-ind_w) ) * exp(pie_wi(+1)) ^2 / exp(pie(+1)) = 0; // (10)

exp(pie_wi) = exp(w_i) / exp(w_i(-1)) * exp(pie); // (11)

exp(c_i) + exp(q_h) * (exp(h_i) - exp(h_i(-1))) + (1+exp(r_bh(-1)))*exp(b_i(-1))/exp(pie) =  
   exp(w_i) * exp(l_i) + exp(b_i)  ;  // (12) budget constraint

(1+exp(r_bh)) * exp(b_i) = exp(m_i) * exp(q_h(+1))   *exp(h_i) * exp(pie(+1));     // (13) borrowing constraint impatient household 


////***********  3) CAPITAL PRODUCERS *****************************************************

exp(K) = (1-deltak) * exp(K(-1)) + ( 1 - kappa_i/2 * (exp(I)*exp(ee_qk)/exp(I(-1)) - 1)^2 ) * exp(I) ; // (14)  

1 = exp(q_k) * ( 1 -  kappa_i/2 * (exp(I)*exp(ee_qk)/exp(I(-1)) - 1)^2  - kappa_i * (exp(I)*exp(ee_qk)/exp(I(-1)) - 1) * exp(I)*exp(ee_qk)/exp(I(-1)) ) 
  + beta_e * exp(lam_e(+1)) / exp(lam_e) * exp(q_k(+1)) *   kappa_i * (exp(I(+1))*exp(ee_qk(+1))/exp(I) - 1) * exp(ee_qk(+1)) * (exp(I(+1))/exp(I))^2 ; // real price of capital (15)


////************  4) ENTREPRENEURS *********************************************************

  (1-a_i)*(exp(c_e) - a_i*exp(c_e(-1)))^(-1) = exp(lam_e);         // CON rescaling  (16) FOC consumption
//        (exp(c_e) - a_i*exp(c_e(-1)))^(-1) = exp(lam_e);         // SENZA rescaling

exp(s_e)  * exp(m_e) * exp(q_k(+1)) * exp(pie(+1)) * (1-deltak) 
     + beta_e * exp(lam_e(+1)) * ( exp(q_k(+1))*(1-deltak) + exp(r_k(+1))*exp(u(+1))
     - ( eksi_1*(exp(u(+1))-1)+eksi_2/2*( (exp(u(+1))-1)^2 ) ) )   = exp(lam_e) * exp(q_k) ;  // (17) FOC capital

exp(w_p) =    ni  * (1-alpha) * exp(y_e) / ( exp(l_pd) * exp(x) ); // (18) FOC labor patient households
exp(w_i) = (1-ni) * (1-alpha) * exp(y_e) / ( exp(l_id) * exp(x) ); // (19) FOC labor impatient households

exp(lam_e) - exp(s_e)  * (1+exp(r_be)) = beta_e * exp(lam_e(+1)) * (1+exp(r_be)) / exp(pie(+1));  // (20) FOC credit demand
exp(r_k)  = eksi_1 + eksi_2 * (exp(u)-1); // (21) FOC utilization rate

exp(c_e) + ((1+exp(r_be(-1))) * exp(b_ee(-1)) / exp(pie) ) +  (exp(w_p)*exp(l_pd) + exp(w_i)*exp(l_id)) + exp(q_k) * exp(k_e) 
   + ( eksi_1*(exp(u)-1)+eksi_2/2*(exp(u)-1)^2 ) * exp(k_e(-1)) = 
    exp(y_e) / exp(x) + exp(b_ee) + exp(q_k) * (1-deltak) * exp(k_e(-1))  ;   // budget constraint entrepreneurs (22)

exp(y_e) = exp(A_e) * (exp(u)*exp(k_e(-1)))^(alpha) * ( exp(l_pd)^ni * exp(l_id)^(1-ni) ) ^ (1-alpha); // production technology (23)

(1+exp(r_be)) * exp(b_ee) = exp(m_e) * exp(q_k(+1))  *exp(pie(+1)) * exp(k_e) * (1-deltak); // borrowing constraint entrepreneurs (24)

exp(r_k) = alpha * exp(A_e) * exp(u)^(alpha-1) * exp(k_e(-1))^(alpha-1) * ( exp(l_pd)^ni * exp(l_id)^(1-ni) ) ^ (1-alpha) /exp(x);  // definition (25)

////*************  5)BANKS ****************************************************************

exp(R_b) = - kappa_kb * ( exp(K_b) / exp(B) - exp(vi) ) * (exp(K_b)/exp(B)) ^2  + exp(r_ib) ; // (26) FOC wholesale branch, where r_ib is policy rate

exp(K_b) * exp(pie) = (1-delta_kb) * exp(K_b(-1)) / exp(eps_K_b) + exp(j_B(-1)) ; // (27) accumulation of bank capital

gamma_b * exp(d_b)  = gamma_p * exp(d_p) ;                                      
gamma_b * exp(b_h)  = gamma_i * exp(b_i) ; 
gamma_b * exp(b_e)  = gamma_e * exp(b_ee); 

exp(b_h) + exp(b_e)  =  exp(d_b) + exp(K_b) ;

/// PRICING in terms of MK ///

// (28) FOC deposit branch:
- 1 + exp(mk_d)/(exp(mk_d)-1)  - exp(mk_d)/(exp(mk_d)-1)  * exp(r_ib)/exp(r_d)  - kappa_d  * ( exp(r_d)/exp(r_d(-1)) - ( exp(r_d(-1)) / exp(r_d(-2)) )^ind_d  )  * exp(r_d)/exp(r_d(-1)) 
  + beta_p * ( exp(lam_p(+1))/exp(lam_p) ) * kappa_d  * ( exp(r_d(+1))/exp(r_d) - ( exp(r_d)/exp(r_d(-1)))^ind_d )   * ( (exp(r_d(+1))/exp(r_d))^2 )   * (exp(d_b(+1))/exp(d_b)) = 0;// 
  
// (29) FOC loan branch, entrepreneurs:
+ 1 - exp(mk_be)/(exp(mk_be)-1)  +  exp(mk_be)/(exp(mk_be)-1)  * exp(R_b)/exp(r_be) - kappa_be * (exp(r_be)/exp(r_be(-1)) - ( exp(r_be(-1)) / exp(r_be(-2)) )^ind_be ) * exp(r_be)/exp(r_be(-1)) 
  + beta_p * ( exp(lam_p(+1))/exp(lam_p) ) * kappa_be * ( exp(r_be(+1))/exp(r_be) - ( exp(r_be)/exp(r_be(-1)))^ind_be ) * ( (exp(r_be(+1))/exp(r_be))^2 ) * (exp(b_e(+1))/exp(b_e)) = 0;//   30
  
// (30) FOC loan brach, impatient households
+ 1 - exp(mk_bh)/(exp(mk_bh)-1)  +  exp(mk_bh)/(exp(mk_bh)-1)  * exp(R_b)/exp(r_bh) - kappa_bh * (exp(r_bh)/exp(r_bh(-1)) - ( exp(r_bh(-1)) / exp(r_bh(-2)))^ind_bh ) * exp(r_bh)/exp(r_bh(-1)) 
  + beta_p * ( exp(lam_p(+1))/exp(lam_p) ) * kappa_bh * ( exp(r_bh(+1))/exp(r_bh) - ( exp(r_bh)/exp(r_bh(-1)))^ind_bh ) * ( (exp(r_bh(+1))/exp(r_bh))^2 ) * (exp(b_h(+1))/exp(b_h)) = 0;//

// (31) overall bank profits:
exp(j_B) = + exp(r_bh)  *  exp(b_h)
           + exp(r_be)  *  exp(b_e) 
           - exp(r_d)   *  exp(d_b)           
           - kappa_d/2  * ( (exp(r_d)/exp(r_d(-1))-1)^2)   * exp(r_d) *exp(d_b) 
           - kappa_be/2 * ( (exp(r_be)/exp(r_be(-1))-1)^2) * exp(r_be)*exp(b_e) 
           - kappa_bh/2 * ( (exp(r_bh)/exp(r_bh(-1))-1)^2) * exp(r_bh)*exp(b_h)
           - kappa_kb/2 * ( (exp(K_b) / exp(B)  - exp(vi) ) ^2) * exp(K_b); //  32

////***********  6)RETAILERS **************************************************************

exp(J_R)  = exp(Y)*(1 - (1/exp(x))    - (kappa_p/2) * (exp(pie) - ( exp(pie(-1)) ^ ind_p * piss ^ (1-ind_p) ))^2 ) ; // (32) aggregate retail profits

// (33) New Keynesian Phillips curve:
1 - exp(eps_y) + exp(eps_y) / exp(x) -    kappa_p * (exp(pie)     - ( exp(pie(-1)) ^ ind_p * piss ^ (1-ind_p) )) * exp(pie) 
    + beta_p*(exp(lam_p(+1))/exp(lam_p))* kappa_p * (exp(pie(+1)) - ( exp(pie)     ^ ind_p * piss ^ (1-ind_p) )) * exp(pie(+1)) * (exp(Y(+1))/exp(Y)) = 0;  // 34

////************  7) AGGREGATION & EQUILIBRIUM  ************************************************

exp(C)              = gamma_p * exp(c_p) + gamma_i * exp(c_i) + gamma_e * exp(c_e);
exp(BH)             = gamma_b * exp(b_h);
exp(BE)             = gamma_b * exp(b_e);
exp(B)              = (exp(BH) + exp(BE)); //
exp(D)              = gamma_p * exp(d_p) ; // oppure: (gamma_b * exp(d_b))
exp(Y)              = gamma_e * exp(y_e); //
exp(J_B)            = gamma_b * exp(j_B);  // 
gamma_e * exp(l_pd) = gamma_p * exp(l_p);
gamma_e * exp(l_id) = gamma_i * exp(l_i);
exp(h)                   = gamma_p * exp(h_p) + gamma_i * exp(h_i); //
exp(K)              = gamma_e * exp(k_e); //
% exp(Y1)             = exp(C) +    1     * (exp(K)-(1-deltak)*exp(K(-1))) + delta_kb * exp(K_b(-1)); //   47
exp(Y1)             = exp(C) +    1     * (exp(K)-(1-deltak)*exp(K(-1))) + 1*(exp(h)-(1-deltah)*exp(h(-1))); //   47
//exp(Y)            = exp(C) + exp(q_k) * (exp(K)-(1-deltak)*exp(K(-1))) + delta_kb * exp(K_b(-1))
//                      + (eksi_1*(exp(u)-1) + eksi_2/2*((exp(u)-1)^2))
//                      + kappa_p/2  * (  exp(pie) - ( exp(pie(-1)) ^ ind_p * piss ^ (1-ind_p) ))^2 * exp(Y)
//                      + kappa_d/2  * ( (exp(r_d(-1))/exp(r_d(-2))-1)^2)   * exp(r_d(-1)) *exp(d_b(-1)) 
//                      + kappa_be/2 * ( (exp(r_be(-1))/exp(r_be(-2))-1)^2) * exp(r_be(-1))*exp(b_e(-1)) 
//                      + kappa_bh/2 * ( (exp(r_bh(-1))/exp(r_bh(-2))-1)^2) * exp(r_bh(-1))*exp(b_h(-1))
//                         ;  //
exp(PIW)            = ( exp(w_p) + exp(w_i) )  / ( exp(w_p(-1)) + exp(w_i(-1)) ) * exp(pie);

////***********  8) TAYLOR RULE & PROFITS CB *****************************************************                                                  

(1+exp(r_ib)) = (1+r_ib_ss)^(1 - rho_ib ) * (1+exp(r_ib(-1)))^rho_ib * (( exp(pie) / piss ) ^phi_pie *  
             (exp(Y1)/exp(Y1(-1)))^phi_y  ) ^ ( 1 - rho_ib ) * (1+e_r_ib) ;// (34) monetary policy

////***********  9) EXOGENOUS PROCESSES ****************************************************12

exp(ee_z)     = 1 - rho_ee_z   *    1          + rho_ee_z   * exp(ee_z(-1))    + e_z;
exp(A_e)      = 1 - rho_A_e    *    1          + rho_A_e    * exp(A_e(-1))     + e_A_e;
exp(ee_j)     = 1 - rho_ee_j   *    1          + rho_ee_j   * exp(ee_j(-1))    - e_j;
exp(m_i)      = (1-rho_mi)     *  m_i_ss       + rho_mi     * exp(m_i(-1))     + e_mi;
exp(m_e)      = (1-rho_me)     *  m_e_ss       + rho_me     * exp(m_e(-1))     + e_me;
exp(mk_d)     = (1-rho_mk_d)   * mk_d_ss       + rho_mk_d   * exp(mk_d(-1))    + e_mk_d;
exp(mk_be)    = (1-rho_mk_be)  * mk_be_ss      + rho_mk_be  * exp(mk_be(-1))   + e_mk_be;
exp(mk_bh)    = (1-rho_mk_bh)  * mk_bh_ss      + rho_mk_bh  * exp(mk_bh(-1))   + e_mk_bh;
exp(ee_qk)    =  1-rho_ee_qk   *    1          + rho_ee_qk  * exp(ee_qk(-1))   + e_qk;
exp(eps_y)    = (1-rho_eps_y)  * eps_y_ss      + rho_eps_y  * exp(eps_y(-1))   + e_y;
exp(eps_l)    = (1-rho_eps_l)  * eps_l_ss      + rho_eps_l  * exp(eps_l(-1))   + e_l;
exp(eps_K_b)  = (1-rho_eps_K_b)*    1          + rho_eps_K_b* exp(eps_K_b(-1)) + e_eps_K_b;
exp(ee_qh)=  1-rho_ee_qh   *    1          + rho_ee_qh  * exp(ee_qh(-1))   + e_qh;
exp(vi)       = (1- rho_vi)*vi_ss                    + e_vi;


////***********  10) AUXILIARY VARIABLES *****************************************************4

rr_e         = exp(lam_e) - beta_e*exp(lam_e(+1))*(1+exp(r_be))/exp(pie(+1));
aux1         = exp(K_b)/exp(B);
exp(bm)      = (exp(b_h(-1))/(exp(b_h(-1))+exp(b_e(-1))) * exp(r_bh(-1)) + exp(b_e(-1))/(exp(b_h(-1))+exp(b_e(-1))) * exp(r_be(-1))) - exp(r_d(-1));
exp(spr_b)   =  0.5*exp(r_bh) + 0.5*exp(r_be) - exp(r_d);                                 //64


////***********  Housing PRODUCERS *****************************************************

exp(h) = (1-deltah) * exp(h(-1)) + ( 1 - kappa_ih/2 * (exp(I_h)*exp(ee_qh)/exp(I_h(-1)) - 1)^2 ) * exp(I_h) ; // (14)  

1 = exp(q_h) * ( 1 -  kappa_ih/2 * (exp(I_h)*exp(ee_qh)/exp(I_h(-1)) - 1)^2  - kappa_ih * (exp(I_h)*exp(ee_qh)/exp(I_h(-1)) - 1) * exp(I_h)*exp(ee_qh)/exp(I_h(-1)) ) 
  + beta_p * exp(lam_p(+1)) / exp(lam_p) * exp(q_h(+1)) *   kappa_ih * (exp(I_h(+1))*exp(ee_qh(+1))/exp(I_h) - 1) * exp(ee_qh(+1)) * (exp(I_h(+1))/exp(I_h))^2 ; // real price of capital (15)


end;


initval;
c_p               =1.5979;
h_p               =-0.0245631;
d_p               =3.57254;
l_p               =-0.0738612;
lam_p             =-1.5979;
J_R               =0.463875;
j_B               =-1.52859;
pie_wp            =8.32735E-16;
c_i               =-1.0838;
h_i               =-3.71877;
b_i               =1.54706;
l_i               =-0.0189437;
lam_i             =1.0838;
s_i               =-2.78576;
pie_wi            =9.02975E-16;
I                 =1.71553;
q_k               =-3.03391E-17;
c_e               =0.667579;
k_e               =4.73146;
l_pd              =-0.0738612;
l_id              =-0.0189437;
b_ee              =3.62249;
y_e               =2.25563;
lam_e             =-0.667579;
s_e               =-4.79487;
u                 =-0.57437;
d_b               =3.57254;
b_h               =1.54706;
b_e               =3.62249;
r_d               =-5.59572;
r_bh              =-5.48472;
r_be              =-4.71577;
R_b               =-5.64415;
K_b               =1.87507;
pie               =8.9463E-16;
x                 =0.182322;
C                 =1.97831;
Y                 =2.25563;
D                 =3.57254;
BE                =3.62249;
BH                =1.54706;
B                 =3.74072;
w_p               =1.24331;
w_i               =-1.00884;
J_B               =-1.52859;
q_h               =5.61245;
K                 =4.73146;
PIW               =8.51348E-16;
r_ib              =-5.04529;
r_k               =-2.68162;
ee_z              =-4.57284E-16;
A_e               =-1.35501E-17;
ee_j              =4.21102E-16;
mk_d              =-0.55043;
mk_be             =0.928379;
mk_bh             =0.159427;
ee_qk             =-6.16294E-17;
m_i               =-0.34249;
m_e               =-1.04982;
eps_y             =1.79176;
eps_l             =1.60944;
eps_K_b           =3.75294E-16;
Y1                =2.54867;
rr_e              =0.00834612;
aux1              =0.154796;
bm                =-5.35942;
spr_b             =-5.8648;
vi                =-1.96053;
interestPol       =2.57585;
interestH         =1.65988;
interestF         =3.58119;
inflation         =0.000000000000089463;
loansH            =154.706;
loansF            =362.249;
output            =254.867;
consumption       =197.831;
investment        =171.553;
deposits          =357.254;
interestDep       =1.4855;
bankcapital       =187.507;
h                 =5.54525;
I_h               =0.940083;
ee_qh             =2.23792E-16;
end;

resid(1);
//// COMPUTE THE SS /////
steady(solve_algo=4);


    
shocks;
//estimated st.dev.
var e_z         = 0.1611^2;
var e_A_e       = 0.0964^2;
var e_j         = 0.4392^2;
var e_me        = 0.0331^2;
var e_mi        = 0.0757^2;
var e_mk_d      = 0.8979^2;
var e_mk_bh     = 0.0075^2;
var e_mk_be     = 0.3835^2;
var e_qk        = 0.2138^2;
var e_r_ib      = 0.0077^2;
var e_y         = 1.4149^2;
var e_l         = 0.3881^2;
var e_eps_K_b   = 0.1565^2;
var e_qh = 0.6077^2;
var e_vi = 0.1102^2;
end;

// stoch_simul(order=1, irf=20, noprint, nograph) interestPol interestH interestF inflation loansH loansF output consumption investment deposits interestDep bankcapital;
stoch_simul(order=1, irf=20, irf_shocks=(e_vi) );
%stoch_simul(order=1, irf=20, irf_shocks=(e_j) );
stoch_simul(order=1, irf=20, irf_shocks=(e_r_ib) );
stoch_simul(order=1, irf=20, irf_shocks=(e_A_e) );
