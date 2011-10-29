function O = o_zajlis(X,t)
    alpha = 1;
    O = zeros(size(X));
    O(1) = X(1) * (1 - alpha * X(2));
    O(2) = X(2) * (X(1) - 1);
end

function O = o_bolniki(X,t)
    mu = 1;
    n = max(size(X));
    nu = 1/(n-2);
    O(1) = - mu * X(1) * sum(X(2:n-1));
    O(2) = mu * X(1) * sum(X(2:n-1)) - nu * X(2);
    for i = 3:n-1
        O(i) = nu * ( X(i-1) - X(i) );
    endfor
    O(n) = nu * X(n-1);
end

function O = o_laser(X, t)
    O = zeros(size(X));
    A = 1;
    B = 1;
    C = 1;
    R = 1.5;
    O(1) = -A * X(1) * X(2) - B * X(1) + R;
    O(2) = A * X(1) * X(2) - C * X(2);
end

function grafi()
    t = linspace(0,100,100)';
    
    X = lsode("o_bolniki", [0.7; 0.01; 0.29], t);
    save bolniki_1.dat X
    
    X = lsode("o_bolniki", [0.7; 0.01; 0; 0; 0; 0; 0; 0; 0.29], t);
    save bolniki_7.dat X
    
    X = lsode("o_laser", [2; 0.1], t);
    save laser_0.dat X
    
    X = lsode("o_laser", [0.5; 1], t);
    save laser_1.dat X
endfunction

grafi()
    