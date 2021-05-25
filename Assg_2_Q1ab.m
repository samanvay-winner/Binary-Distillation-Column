% This will use the code from Assg 1 & build up on that.
clc
clear

x_top = input("The top composition value (zd): ");
x_bottom = input("The bottom composition value (zw): ");
x_feed = input("The feed composition value (zf): ");
reflux_ratio = input("The reflux ratio (R): ");
boilup_ratio = input("The boilup ratio (S): ");

% Will assume number of moles of feed = 1, since our results aren't
% dependent on that.

Tc_prop=537.32; %K
Pc_prop= 51.78;  %bar

Tc_water = 647.3; %K
Pc_water = 220.9; %bar

global R;
R = 0.0832;     %L.bar/mol.K  (0.0821(in L.atm/mol.K)*1.01325(bar/atm))
global a_water;
a_water = ((27/64)*(R*Tc_water)^2)/(Pc_water);  %L^2.bar/mol^2
global b_water;
b_water = (R*Tc_water)/(8*Pc_water); %L/mol
global a_prop;
a_prop = ((27/64)*(R*Tc_prop)^2)/(Pc_prop);  %L^2.bar/mol^2 
global b_prop;
b_prop = (R*Tc_prop)/(8*Pc_prop);  %L/mol
global A12;
A12 = 2.576;  %ln(gamma_inf) for propane
global A21;
A21 = 1.201;  %ln(gamma_inf) for water

B=1441.629;
C=-74.299;
A=4.87601;

P = 1; %bar

x_vec = [0: 0.02: 1];
y_vec = [];
    
for i = 1:length(x_vec)
    x = x_vec(i);
    y = calc_y_from_x(x, 'p', P);
    y_vec = [y_vec, y];
end;

ylim([0, 1]);

hold on
plot(x_vec, y_vec);
plot(x_vec, x_vec, 'r');


%Point of intersection of the 3 lines - feel line & the 2 operating lines
x_intersection = (boilup_ratio*x_top + (reflux_ratio +1)*x_bottom)/(boilup_ratio + reflux_ratio + 1);
y_intersection = ((boilup_ratio + 1)*x_top + reflux_ratio*x_bottom)/(boilup_ratio + reflux_ratio + 1);


top_line = @(x) (reflux_ratio*x + x_top)/(reflux_ratio + 1);
bottom_line = @(x) ((boilup_ratio+1)*x - x_bottom)/(boilup_ratio);
%feed_line = @(x) (x_feed*x + 1)/(x_feed-1);

%m_feed = x_feed/(x_feed- 1);
%x_intersection = fsolve(@(x) (P*(m_feed*x - x_feed*(m_feed-1))) - x*exp(A12/(1+ (A12*x/(A21*(1-x))))^2), x_feed)
%y_intersection = feed_line(x)

x_top_range = [x_intersection, x_top];
x_bottom_range = [x_bottom, x_intersection];
%feed_range = [min(x_intersection, x_feed), max(x_intersection, x_feed)];


fplot(top_line, x_top_range);
fplot(bottom_line, x_bottom_range);
%fplot(feed_line, feed_range)


%Starting point of the iterations (McAbe Thiele)
count = 0;
x_val = x_top;
y_val = x_top; %It'll be x_top only, btw!! Can use that itself directly!!

%y_from_x = @(x) exp(A12/(1+ (A12*x/(A21*(1-x)))^2))*x/P ;
%x_from_y = @(x, y) P*y/exp(A12/(1+ (A12*x/(A21*(1-x)))^2));

%solve_for_hori_inter = @(x) P*y_val - x*exp(A12/(1+ (A12*x/(A21*(1-x))))^2);
%solve_for_vert_inter = @(y) y - (reflux_ratio*x_new + x_top)/(reflux_ratio + 1)

%Horizontal line => y = fixed, x = changing
x_intersections = [];
y_intersections = [];

x_new = x_val;
y_new = y_val;

while(x_new >= x_bottom)
    solve_for_hori_inter = @(x) P*y_val - x*exp(A12/(1+ (A12*x/(A21*(1-x))))^2);
    x_new = fsolve(solve_for_hori_inter, x_val)
    line([x_new, x_val], [y_new, y_new]);

    x_val = x_new;
    
    if(x_new > x_intersection)
        y_new = (reflux_ratio*x_new + x_top)/(reflux_ratio + 1);
    else
        y_new = ((boilup_ratio+1)*x_new - x_bottom)/(boilup_ratio);
    end
    
    count = count + 1;
    line([x_new, x_new], [y_val, y_new])
    y_val = y_new;
    
%    plot(x_new, y_val, '*')
%    plot(x_new, y_new, '*')

%Ignore the below line...
%x_new = fsolve(x_from_y(x_val, y_val), x_val)  %This gave a velue of x_new = -1.7114e-09  :(( Pretty wrong

end

display(count);

%figure; % To generate the next diagram in a new plot window

min_trays = 0;
x_val = x_top;
y_val = x_top;
x_new = x_top;
y_new = x_top;

while(x_new >= x_bottom)
    solve_for_hori_inter = @(x) P*y_val - x*exp(A12/(1+ (A12*x/(A21*(1-x))))^2);
    x_new = fsolve(solve_for_hori_inter, x_val)
    line([x_new, x_val], [y_new, y_new], 'LineStyle','-.', 'Color','g');

    x_val = x_new;
    
    y_new = x_new;
    
    min_trays = min_trays + 1;
    line([x_new, x_new], [y_val, y_new], 'LineStyle','-.', 'Color','g')
    y_val = y_new;
    
%    plot(x_new, y_val, '*')
%    plot(x_new, y_new, '*')

%Ignore the below line...
%x_new = fsolve(x_from_y(x_val, y_val), x_val)  %This gave a velue of x_new = -1.7114e-09  :(( Pretty wrong

end

display(min_trays)

function y = calc_y_from_x(x, subs, P) %subs = 'p' or 'w'. P = Total pressure
    %fuga_liq_array = find_fuga_liq(P, x, 'p')
    [gamma1, gamma2] = find_gamma(x, subs);
        
%   y = fuga_liq_array*gamma1*x/P;
    if(subs == 'p')
        y = gamma1*x/P;
    elseif(subs == 'w')
        y = gamma2*x/P;
    else
        error("The substance isn't 'p' or 'w'. Kindly correct it.");
    end
    
end

function [fuga_liq_array] = find_fuga_liq(P, n, subs)
    global R a_prop b_prop a_water b_water

    %A, B, C = Antoine's coeff. parameters
    B=1441.629;
    C=-74.299;
    A=4.87601;
    fuga_coeff_array = [];
    fuga_liq_array = [];
    
%     temp_range = 350:1:353;
        
%     for i = 1:length(temp_range)
%         T = temp_range(i)
        
        T = 400;
        display(size(T));
        
        if(subs == 'p')
            a = a_prop;
            b = b_prop;
        elseif(subs == 'w')
            a = a_water;
            b = b_water;
        else
            error("The substance isn't 'p' or 'w'. Kindly correct it.");
        end

               
            P_sat = 10^(A - (B./(T + C)))
            express = [1, -n.*(b+R.*T/P_sat), (n.^2)*(a/P_sat), -(n.^3)*(a*b)/P_sat];
            V = roots(express)
            vol = min(V(V>0)) %Liquid volume at given conds.
                       
            Z = P_sat*vol/(n*R.*T)
            fuga_coeff = (vol/(Z*(vol-n*b)))*exp(n*b/(vol - n*b) - 2*a*n/(vol*R.*T))
            fuga_coeff_array = [fuga_coeff_array, fuga_coeff];
    % This is the fugacity coefficient (phi_i_sat) at T, P_sat
            Poyntingfactor = exp((vol*(P-P_sat))/(n*R.*T));
            
            fuga_liq = P_sat*fuga_coeff*Poyntingfactor
            fuga_liq_array = [fuga_liq_array, fuga_liq];
%     end
end

function [gamma1, gamma2] = find_gamma(x, subs)  %Code: 'p' for propane (1) and 'w' for water (2)
    global A12 A21;
    
    if(subs == 'p')
        x1 = x;
        x2 = 1-x;
    elseif(subs == 'w')
        x1 = 1-x;
        x2 = x;
    else
        error("The substance isn't 'p' or 'w'. Kindly correct it.");
    end
    
    ln_gamma1 = A12/(1+ (A12*x1/(A21*x2)))^2;
    ln_gamma2 = A21/(1+ (A21*x2/(A12*x1)))^2;
    
    %Both of the next 2 lines should be uncommented only if gamma_i is reqd
    %not ln(gamma_i)
     gamma1 = exp(ln_gamma1);
     gamma2 = exp(ln_gamma2);

end

function find_phi(T, P, n, subs) %n is the number of moles, it could be one of xi's or yi's
%Note that a, b are not values, but function handles    
    global R a_prop b_prop a_water b_water;
    
    if(subs == 'p')
        ai = a_prop;
        bi = b_prop;
        y_prop = n; %Remember, total moles = 1 
    elseif(subs == 'w')
        ai = a_water;
        bi = b_water;
        y_prop = 1 - n;  %Remember, total moles = 1 
    else
        error("The substance isn't 'p' or 'w'. Kindly correct it.");
    end
    
    eqn = [P, -n*(P*bi+R*T), ai, -ai*bi];
    V = roots(eqn);
    Vg = max(V);
    
    %f = @(v) ((1/(v-n*b) + b*n/(v-n*b)^2 - 2*a*n/(v*R*T)) - 1/v);
    %ln_phi = int(f, vg, inf) - ln(P*vg/(n*R*T));  %Note that the latter term has only vg, and not v.
    
    b = @(y_prop) y_prop*b_prop + (1-y_prop)*b_water;
    a = @(y_prop) (y_prop*sqrt(a_prop) + (1-y_prop)*sqrt(a_water))^2;

    ln_phi = @(y_prop) log(Vg/(Vg-n*b(y_prop))) + bi/(Vg-n*b(y_prop)) - log(P*Vg/(n*R*T)) -2*sqrt(ai)*(sqrt(a(y_prop)))/(R*T);
    %Note: if not converging well, use (y_prop*sqrt(a_prop) +
    %(1-y_prop)*sqrt(a_water)) in place of the last term sqrt(a).
    
%    soln = fzero(ln_phi, 0);
    %fzero wouldn't work at all, actually (we don't need to drive a, b, ln_phi to 0). 
    %Something else is reqd. Tried fsolve, but didn't quite work... Also, constraint on y_prop to be
    %between 0 and 1 is reqd, I think...
   
end