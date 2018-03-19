/* **************************************************** */
//      REPLICATION SOMMER ET SULLIVAN (AER, 2017)      //
/*                  February 2018                       */
//               Alexandre Gaillard                     //
/* **************************************************** */




/****************************************/
// FUNCTION CONSUMPTION (OR MAX WEALTH) //
/****************************************/

double progtaxfun(double income_deduct){

int tt;
double taxes_paid, stop;

taxes_paid = 0;
stop = 0;
tt = 0;

while(stop == 0){
    if(income_deduct >= threshold[tt] && income_deduct < threshold[tt+1]){
        taxes_paid = taxes_paid + prog_rate[tt]*(income_deduct - threshold[tt]);
        stop = 1;
    }else{
        taxes_paid = taxes_paid + prog_rate[tt]*(threshold[tt+1] - threshold[tt]);
    }
    
    tt++;
    
    if(tt > 5){printf("mistake %f %d %f %f %f", taxes_paid, tt-1, income_deduct, threshold[tt], threshold[tt-1]); getchar();}
}

return(taxes_paid);
}





double COHfun(const int k, const int y, const int h, const int hp, const int s){
double cashonhand, itemized, payroll_taxes, fixed_cost, M, frac_rented, frac_occupied, Idummy, Ik, standard, renthousing, housing, income, income_deduct, taxes_paid;

Idummy = 0;
Ik = 0;
income_deduct= 0;
taxes_paid = 0;
income = 0;

if(K[inxKH(k,h)] >= 0){Ik = 1;}else{Ik = 0;}
if(H[h] != H[hp]){Idummy = 1;}else{Idummy = 0;}


// COMPUTE DISPOSAL INCOME //

// renter //
if(H[hp] == 0){
    income = wage*Y[y] + Ik*irateS*K[inxKH(k,h)];
    income_deduct = income - xipar - epar;
    
    M = 0.0;
    fixed_cost = 0.0;
}

// landlord //
if(H[hp] > 0 && H[hp] > S[s]){
    frac_rented = ((H[hp] - S[s])/H[hp]);
    frac_occupied = 1.0 - frac_rented;
    
    income = wage*Y[y] + Ik*irateS*K[inxKH(k,h)] + rent*(H[hp] - S[s]);
    
    itemized = - epar  // personal exemption
        + policy_mortgage*frac_occupied*(tau_m*irateB)*K[inxKH(k,h)]*(1 - Ik) // ownocc_int_deduct*frac_occupied*( tauM*(intrate(ir)+kappa)*min_dm0 )
        - frac_occupied*tau_h*price*H[hp]
        + frac_rented*((tau_m*irateB)*K[inxKH(k,h)]*(1 - Ik))    // frac_rented*( tauM*(intrate(ir)+kappa)*min_dm0 )
        - frac_rented*tau_h*price*H[hp];

    itemized = itemized - (tau_LL + deltaH)*price*(H[hp] - S[s]);
    
    if(fabs(itemized) > epar + xipar){
        income_deduct = income + itemized;
    }else{
        income_deduct = income - epar - xipar;
    }
    
    M = deltaH*H[hp];
    fixed_cost = phipar;
}

// homeowner //
if(H[hp] > 0 && H[hp] <= S[s]){

    //if(S[s] > H[hp]){printf("MISTAKE, s > hp, %d %d", s, hp);getchar();}
    
    income = wage*Y[y] + Ik*irateS*K[inxKH(k,h)] + rent*(H[hp] - S[s]);
    
    itemized = - epar  // personal exemption
        + policy_mortgage*(tau_m*irateB)*K[inxKH(k,h)]*(1 - Ik) // ownocc_int_deduct*tauM*(intrate(ir)+kappa)*min_dm0
        - tau_h*price*H[hp];

    if(fabs(itemized) > epar + xipar){
        income_deduct = income + itemized;
    }else{
        income_deduct = income - epar - xipar;
    }
    
    M = deltaH*H[hp];
    fixed_cost = 0.0;
}


income_deduct = max(0,income_deduct);

// taxes (piecewise function) //
taxes_paid = progtaxfun(income_deduct);
if(taxes_paid < 0){printf("mistake taxes_paid < 0");getchar();}

payroll_taxes = tau_p*wage*Y[y];


cashonhand = wage*Y[y]
            - price*H[hp]
            + price*H[h]
            + (1+irateB)*K[inxKH(k,h)]*(1 - Ik)
            + (1+irateS)*K[inxKH(k,h)]*Ik
            - Idummy*tau_s*price*H[h]
            - Idummy*tau_b*price*H[hp]
            - rent*(S[s]-H[hp])
            - taxes_paid
            - payroll_taxes
            - M*price
            - fixed_cost
            - tau_h*price*H[hp];



return(cashonhand); // may be negative for some values

}


