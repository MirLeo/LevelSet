function L = Lim(x,y)
% Superbee limiter
% R. Leveque, Finite Volume Methods for Hyperbolic Problems,...
% Cambridge University Press, Cambridge, 2002.
if x*y > 0
    
    if (abs(y) > 2*abs(x) || abs(y) < abs(x)/2)
        
        L = 2 * sign(x) * min(abs(x),abs(y));
    else
        
        L = sign(x) * min(abs(x),abs(y));
    end
    
else
    
    L = 0;
    
end