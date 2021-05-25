% This file finds the Minimum Reflux ratio
% This will use the code from Assg 1 & build up on that.
clc
clear

%Replace these with input commands later!

x_top = input("The top composition value (zd): ");
x_bottom = input("The bottom composition value (zw): ");
x_feed = input("The feed composition value (zf): ");


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
end

ylim([0, 1]);

hold on
plot(x_vec, y_vec);
plot(x_vec, x_vec, 'r');


%Point of intersection of the 3 lines - feel line & the 2 operating lines
%x_intersection = (boilup_ratio*x_top + (reflux_ratio +1)*x_bottom)/(boilup_ratio + reflux_ratio + 1);
%y_intersection = ((boilup_ratio + 1)*x_top + reflux_ratio*x_bottom)/(boilup_ratio + reflux_ratio + 1);

m_feed = x_feed/(x_feed- 1);
feed_line = @(x) m_feed*x - x_feed*(m_feed-1);

x_intersection = fsolve(@(x) (P*(m_feed*x - x_feed*(m_feed-1))) - x*exp(A12/(1+ (A12*x/(A21*(1-x))))^2), x_feed);
y_intersection = feed_line(x_intersection);

x_top_range = [x_intersection, x_top];
x_bottom_range = [x_bottom, x_intersection];
feed_range = [min(x_intersection, x_feed), max(x_intersection, x_feed)];

slope = (y_intersection - x_top)/(x_intersection - x_top); % Of the Enriching section line; = y2-y1/(x2-x1)

min_reflux_ratio = slope/(1-slope);

display(min_reflux_ratio)

% The boilup line can be made pretty easily too!
% We know the two points it passes through => Slope = known
% And the slope = (B+1)/B. Hence, B = slope/(slope-1)

slope1 = (y_intersection - x_bottom)/(x_intersection - x_bottom);
corr_boilup_ratio = slope1/(slope1-1);


top_line = @(x) (min_reflux_ratio*x + x_top)/(min_reflux_ratio + 1);
bottom_line = @(x) slope1*(x - x_bottom) + x_bottom;
feed_line = @(x) m_feed*x - x_feed*(m_feed-1);



fplot(top_line, x_top_range);
fplot(bottom_line, x_bottom_range);
fplot(feed_line, feed_range)


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