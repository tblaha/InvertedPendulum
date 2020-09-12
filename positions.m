function [cgpos, endpoints] = positions(s, th, l, n)
    
    cgpos = [s; 0] + [sin(th); cos(th)] .* l * (0.5*eye(n) + triu(ones(n), 1));
    endpoints = [s; 0] + [[0; 0], [sin(th); cos(th)].* l * triu(ones(n))];
    
end