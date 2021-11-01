 function [ xt, diff, fn, flagin] = nonmonotoneArmijocriterion( fn, fmax, t, gamma, gtd, x,d,p, xt,diff, A)
%% Nonmonotone Armijo Criterion
%
% Syntax: 
%        [ xt, diff, fn, flagin] = nonmonotoneArmijocriterion( fn, fmax, t, gamma, gtd, x,d,p, xt,diff, A)
%
% Reference: G. Birgin, J. M. Martinez, and M. Raydan. Journal of Statistical Software, 2014 

    flagin=1;
    
    while fn > fmax + t*gamma*gtd   %Test    -like criterion
        t = t/2;    
        if (t< 1e-8)                  %Will hardly change x on update
            flagin = -2;
            return
        end       
        xt = x+t*d;    
        diff = A*xt-p;                                
        normdiff = norm(diff,2);                                                 
        fn = 0.5*normdiff^2;
    end   
end