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
//            STICKY BANK RATES, PRICES & WAGES ?la Rotemberg with indexation to both past and st.st. inflation
// 9 blocks: 1) PATIENTs  2) IMPATIENTs  3) CAPITAL PRODUCERS  4) ENTREPRENEURS   5) BANKS   6) RETAILERS        
//           7) LABOR MKT WITH ONE UNION FOR EACH LABOR TYPE  8) AGGREGATION & EQUILIBRIUM   9) MONETARY POLICY  

//    AGGIUNTA DEL 16 APRILE: PROBLEMI INTERTEMPORALI TUTTI IN TERMINI REALE

// THIS VERSION: April 20th, 2009
// BBS has all items dated (t) - Banking capital in real terms is defined as: K_B(t)/p(t)
// Banking profits are defined in the model code at time (t) - All equations in exp form
// This is the full version (monop. competitive banking sector & BANK capital ?la Gerali with sticky rates) 



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

//**************************************************************************
// Replication Variables                                                  //*
//      68  69    70    71   72    73    74      75        76        77          78       79       80          81 
        dc dinv infobs wages hp loansH loansF deposits interestPol interestH interestF output interestDep bankcapital;
//**************************************************************************


varexo  e_z e_A_e e_j e_me e_mi e_mk_d e_mk_bh e_mk_be e_qk e_r_ib e_y e_l e_eps_K_b  ;  //13 varexo



parameters  
            beta_p j phi beta_i m_i_ss beta_e m_e_ss alpha eksi_1 eksi_2    // HOUSEHOLDS & ENTREPRENEURS
            h a_i a_p a_e gamma_p gamma_i gamma_e    ni                     // HOUSEHOLDS & ENTREPRENEURS
            eps_l_ss kappa_w                                                // HOUSEHOLDS (labor params)
            eps_d eps_bh eps_be                                             // BANKS 
            mk_d_ss mk_bh_ss mk_be_ss r_be_ss  r_bh_ss r_k_ss               // BANKS (SS)
            gamma_b beta_b delta_kb vi kappa_kb                             // BANKS 
            eps_y_ss kappa_p ind_p ind_w                                    // RETAILERS
            kappa_i kappa_d kappa_be kappa_bh deltak                        // OTHERS
            ind_d ind_be ind_bh                                             // OTHERS
            rho_ib phi_pie phi_y                                            // POLICY
            piss  r_ib_ss                                                   // STEADY STATE
            rho_ee_z rho_A_e rho_ee_j rho_mi rho_me rho_eps_y               // SHOCKS
            rho_mk_d rho_mk_be rho_mk_bh rho_ee_qk rho_eps_l rho_eps_K_b    // SHOCKS
            ;
            

% *********************			
% CALIBRATED PARAMETERS
% *********************

beta_p       = 0.9943;                                                     % discount factor patient households
beta_i       = 0.975;                                                      % discount factor impatient households     
beta_b       = beta_p;                                                     % discount factor bankers (not used in this version of the model)
beta_e       = beta_i;                                                     % discount factor entrepreneurs
j            = 0.2;                                                        % weight of housing in utility function
phi          = 1.0;                                                        % inverse Frisch elasticity of labor supply
m_i_ss       = 0.7  ;                                                      % loan-to-value ratio impatient households
m_e_ss       = 0.35 ;                                                      % loan-to-value ratio entrepreneurs
alpha        = 0.250;                                                      % capital share in the production function
eps_d        = -1.46025;                                                   % elast. of subst. of deposits 
eps_bh       = 2.932806;
eps_be       = 2.932806; 
 
mk_d_ss      = eps_d   / (eps_d  - 1) ;                                    % steady state markdown on D (ok if eps_d<0; if eps_d>0 it should be eps_d/(eps_d+1) )
mk_bh_ss     = eps_bh  / (eps_bh - 1) ;                                    % steady state markup on loans to I
mk_be_ss     = eps_be  / (eps_be - 1) ;                                    % steady state markup on loans to E
book_ss      = 0; %-35                                                       % steady state value of (B-M)/D in the bank balance shhet
eps_y_ss     = 6;                                                          % 
eps_l_ss     = 5;                                                          % 
gamma_p      = 1;                                                          % shares of patient households
gamma_i      = 1; //1/3;                                                   % shares of impatient households
ni           = 0.8;                                                        % wage share of patient households
gamma_b      = 1; //0.10;												   % shares of bankers
gamma_e      = 1; //1 - gamma_p - gamma_i;								   % shares of entrepreneurs
deltak       = 0.025;                                                      % depreciation rate for physical capital
piss         = 1;                                                          % steady state gross inflation rate
r_ib_ss      = (piss/beta_p - 1) * (eps_d-1)/eps_d ;                       % steady state gross nominal interest rate 
r_be_ss      = r_ib_ss*eps_be/(eps_be-1) ;								   % steady state interest rate on loans to E
r_bh_ss      = r_ib_ss*eps_bh/(eps_bh-1) ;								   % steady state interest rate on loans to H
r_k_ss       = -(1-deltak)-m_e_ss*(1-deltak)*piss/beta_e*(1/(1+r_be_ss)-beta_e/piss)+1/beta_e;                       % steady state rental rate of capital
h            = 1;                                                          % fixed supply housing
eksi_1       = r_k_ss;                                                     % capital utilization cost parameter
eksi_2       = 0.1*r_k_ss;                                                 % capital utilization cost parameter

vi           = 0.12;                                                      % Banking Capital ratio over Loans (Basel II)
eps_b        = mean([eps_bh,eps_be]);
% delta_kb     = r_ib_ss/vi * (eps_d - eps_b + vi*(eps_b-1))/((eps_b-1)*(eps_d-1));  % j_B with terms in r_ib                            
delta_kb     = r_ib_ss/vi * (eps_d - eps_b + vi*eps_d*(eps_b-1))/((eps_b-1)*(eps_d-1));                               

ind_d        = 0.0;                   % indexation deposit rates
ind_be       = 0.0;                   % indexation rates on loans to firms
ind_bh       = 0.0;                   % indexation rates on loans to households


            
//%------------------------------------------------------------
//% Model equations
//%------------------------------------------------------------

model;


// Model code:

////***********   1) PATIENT HHs ********************************************************6

  (1-a_i)*exp(ee_z)*(exp(c_p) - a_i*exp(c_p(-1)))^(-1) = exp(lam_p); // CON rescaling  (1) where a_p = a_e = a_i

j * (exp(ee_j))  / exp(h_p) - exp(lam_p) * exp(q_h) + beta_p * exp(lam_p(+1)) * exp(q_h(+1))   = 0; // (2) 

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

  (exp(c_e) - a_i*exp(c_e(-1)))^(-1) = exp(lam_e);         // CON rescaling  (16) FOC consumption
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

exp(R_b) = - kappa_kb * ( exp(K_b) / exp(B) - vi ) * (exp(K_b)/exp(B)) ^2  + exp(r_ib) ; // (26) FOC wholesale branch, where r_ib is policy rate

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
           - kappa_kb/2 * ( (exp(K_b) / exp(B)  - vi ) ^2) * exp(K_b); //  32

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
h                   = gamma_p * exp(h_p) + gamma_i * exp(h_i); //
exp(K)              = gamma_e * exp(k_e); //
% exp(Y1)             = exp(C) +    1     * (exp(K)-(1-deltak)*exp(K(-1))) + delta_kb * exp(K_b(-1)); //   47
exp(Y1)             = exp(C) +    1     * (exp(K)-(1-deltak)*exp(K(-1))) ; //   47
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
exp(ee_j)     = 1 - rho_ee_j   *    1          + rho_ee_j   * exp(ee_j(-1))    + e_j;
exp(m_i)      = (1-rho_mi)     *  m_i_ss       + rho_mi     * exp(m_i(-1))     + e_mi;
exp(m_e)      = (1-rho_me)     *  m_e_ss       + rho_me     * exp(m_e(-1))     + e_me;
exp(mk_d)     = (1-rho_mk_d)   * mk_d_ss       + rho_mk_d   * exp(mk_d(-1))    + e_mk_d;
exp(mk_be)    = (1-rho_mk_be)  * mk_be_ss      + rho_mk_be  * exp(mk_be(-1))   + e_mk_be;
exp(mk_bh)    = (1-rho_mk_bh)  * mk_bh_ss      + rho_mk_bh  * exp(mk_bh(-1))   + e_mk_bh;
exp(ee_qk)    =  1-rho_ee_qk   *    1          + rho_ee_qk  * exp(ee_qk(-1))   + e_qk;
exp(eps_y)    = (1-rho_eps_y)  * eps_y_ss      + rho_eps_y  * exp(eps_y(-1))   + e_y;
exp(eps_l)    = (1-rho_eps_l)  * eps_l_ss      + rho_eps_l  * exp(eps_l(-1))   + e_l;
exp(eps_K_b)  = (1-rho_eps_K_b)*    1          + rho_eps_K_b* exp(eps_K_b(-1)) + e_eps_K_b;


////***********  10) AUXILIARY VARIABLES *****************************************************4

rr_e         = exp(lam_e) - beta_e*exp(lam_e(+1))*(1+exp(r_be))/exp(pie(+1));
aux1         = exp(K_b)/exp(B);
exp(bm)      = (exp(b_h(-1))/(exp(b_h(-1))+exp(b_e(-1))) * exp(r_bh(-1)) + exp(b_e(-1))/(exp(b_h(-1))+exp(b_e(-1))) * exp(r_be(-1))) - exp(r_d(-1));
exp(spr_b)   =  0.5*exp(r_bh) + 0.5*exp(r_be) - exp(r_d);                                 //64


//**************************************************************************
// Measurement(observation) equations //
dc = 100*(0.008522411 + C - C(-1));
dinv = 100*(0.006191921 + I - I(-1));
infobs = 100*(0.008813708 + pie);
wages = 100*(0.117237288 + PIW);
hp = 100*(-0.000052664 + q_h - q_h(-1));
loansH = 100*(0.038272095 + BH - BH(-1));
loansF = 100*(0.013037611 + BE - BE(-1)); 
deposits = 100*(0.024311868 + D - D(-1));   
interestPol   = 400*exp(r_ib); 
interestH = 400*exp(r_bh);
interestF = 400*exp(r_be);  
interestDep = 400*exp(r_d);                                       

//**************************************************************************

output     = Y1*100; 
bankcapital = 100*K_b;                           


end;


varobs dc dinv infobs wages hp loansH loansF deposits interestPol interestH interestF interestDep;

initval;

c_p      	=	-0.122029	;
h_p      	=	-0.0535366	;
d_p      	=	0.840797	;
l_p      	=	-0.262977	;
lam_p    	=	0.122029	;
J_R      	=	-1.52345	;
j_B      	=	-3.79369	;
c_i      	=	-1.91801	;
h_i      	=	-2.95404	;
b_i      	=	0.164689	;
l_i      	=	-0.0581332	;
lam_i    	=	1.91801	;
s_i      	=	-2.57912	;
I        	=	-1.94962	;
c_e      	=	-2.19923	;
k_e      	=	1.73925	;
l_pd     	=	-0.262977	;
l_id     	=	-0.0581332	;
b_ee     	=	0.313687	;
y_e      	=	0.268308	;
lam_e    	=	2.19923	;
s_e      	=	-2.29789	;
d_b      	=	0.940797	;
b_h      	=	0.164689	;
b_e      	=	0.313687	;
r_d      	=	-5.16157	;
r_bh     	=	-4.26486	;
r_be     	=	-4.26486	;
R_b      	=	-4.43360773	;
K_b      	=	-1.57284	;
pie      	=	1.16E-14	;
x        	=	0.182322	;
C        	=	0.133577	;
Y        	=	0.268308	;
D        	=	0.840797	;
BE       	=	0.313687	;
BH       	=	0.164689	;
B        	=	0.935107	;
w_p      	=	-0.161863	;
w_i      	=	-1.753	;
J_B      	=	-3.79369	;
q_h      	=	3.48936	;
K        	=	1.73925	;
r_ib     	=	log(r_ib_ss)   	;
r_k      	=	-3.03956	;
mk_d    	=	log(mk_d_ss)   	;
mk_be   	=	log(mk_be_ss)  	;
mk_bh   	=	log(mk_bh_ss)	;
m_i     	=	log(m_i_ss)	;
m_e     	=	log(m_e_ss)	;
Y1       	=	Y	;
rr_e     	=	exp(s_e)	;
aux1     	=	exp(s_i)    ;
bm       	=	-4.88329	;
spr_b    	=	-4.99645	;
eps_y    	=	log(eps_y_ss);
eps_l  		=   log(eps_l_ss);
end;


shocks;

var e_z;  stderr 0.0144;
var e_A_e;  stderr 0.0062;
var e_j;  stderr 0.0158;
var e_me;  stderr 0.0034;
var e_mi;  stderr 0.0023;
var e_mk_d;  stderr 0.0488;
var e_mk_bh;  stderr 0.0051;
var e_mk_be;  stderr 0.1454;
var e_qk;  stderr 0.0128;
var e_r_ib;  stderr 0.0018;
var e_y;  stderr 1.0099;
var e_l;  stderr 0.3721;
var e_eps_K_b;  stderr 0.050;

end;

estimated_params;

stderr e_z, inv_gamma_pdf, 0.01, 0.05;
stderr e_j, inv_gamma_pdf, 0.01, 0.05;
stderr e_me, inv_gamma_pdf, 0.01, 0.05;
stderr e_mi, inv_gamma_pdf, 0.01, 0.05;
stderr e_mk_d, inv_gamma_pdf, 0.01, 0.05;
stderr e_mk_bh, inv_gamma_pdf, 0.01, 0.05;
stderr e_mk_be, inv_gamma_pdf, 0.01, 0.05;
stderr e_A_e, inv_gamma_pdf, 0.01, 0.05;
stderr e_qk, inv_gamma_pdf, 0.01, 0.05;
stderr e_r_ib, inv_gamma_pdf, 0.01, 0.05;
stderr e_y, inv_gamma_pdf, 1.0, 0.05;
stderr e_l, inv_gamma_pdf, 0.40, 0.05;
stderr e_eps_K_b, inv_gamma_pdf, 0.05, 0.05;

rho_ee_z, beta_pdf, 0.8, 0.10;	    // consumption preference
rho_ee_j, beta_pdf, 0.8, 0.10; 	    // housing preference
rho_me, beta_pdf, 0.8, 0.10;        // Firms' LTV
rho_mi, beta_pdf, 0.8, 0.10;        // HHs' LTV
rho_mk_d, beta_pdf, 0.8, 0.10;	    // deposits markup
rho_mk_bh, beta_pdf, 0.8, 0.10;	    // HHs loans markup
rho_mk_be, beta_pdf, 0.8, 0.10;	    // firms loans markup
rho_A_e, beta_pdf, 0.8, 0.10;       // technology
rho_ee_qk, beta_pdf, 0.8, 0.10;	    // investment efficiency
rho_eps_y, beta_pdf, 0.8, 0.10;	    // p mark-up 
rho_eps_l, beta_pdf, 0.8, 0.10;	    // w mark-up
rho_eps_K_b, beta_pdf, 0.8, 0.10;	// balance sheet

kappa_p, gamma_pdf, 50, 20;    
kappa_w, gamma_pdf, 50, 20;      
kappa_i, gamma_pdf, 2.5, 1.0;     
kappa_d, gamma_pdf, 10, 2.5;     
kappa_be, gamma_pdf, 3, 2.5;	
kappa_bh, gamma_pdf, 6, 2.5;	
kappa_kb, gamma_pdf, 10, 5.0;	
phi_pie, gamma_pdf, 2.0, 0.5;     
rho_ib, beta_pdf, 0.50, 0.15;      
phi_y, normal_pdf, 0.10, 0.15;       
ind_p, beta_pdf, 0.50, 0.15;       
ind_w, beta_pdf, 0.50, 0.15;      
a_i, beta_pdf, 0.50, 0.10;	   // degree of habit formation: p HHs, a_p = a_i     
//a_e=0; // degree of habit formation: Es         
// a_p=0;         
end;

    
estimation(datafile=dataGerali,mh_replic=100000,mh_nblocks=10,mh_drop=0.5,mh_jscale=0.0030,load_mh_file,filtered_vars,mode_compute=6,smoother) 
output C I pie wages loansH loansF deposits interestPol interestH interestF interestDep bankcapital; 

//stoch_simul(order=1, irf=20, noprint, nograph) interestPol interestH interestF inflation loansH loansF output consumption investment deposits interestDep bankcapital;
stoch_simul(order=1, irf=20, irf_shocks=(e_eps_K_b) );
stoch_simul(order=1, irf=20, irf_shocks=(e_j) );
stoch_simul(order=1, irf=20, irf_shocks=(e_r_ib) );