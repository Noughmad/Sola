global A

function O = o_zajlis(X,t)
    global A
    gamma = A;
    O = zeros(size(X));
    O(1) = X(1) * (1 - X(2));
    O(2) = gamma * X(2) * (X(1) - 1);
end

function O = o_bolniki(X,t)
    global A
    mu = A/2;
    n = max(size(X));
    nu = (n-2)/2;
    O(1) = - mu * X(1) * sum(X(2:n-1));
    O(2) = mu * X(1) * sum(X(2:n-1)) - nu * X(2);
    for i = 3:n-1
        O(i) = nu * ( X(i-1) - X(i) );
    endfor
    O(n) = nu * X(n-1);
end

function O = o_laser(X, t)
    global A
    R = A;
    O = zeros(size(X));
    O(1) = -X(1) * X(2) - X(1) + R;
    O(2) = X(1) * X(2) - X(2);
end

function grafi()
    global A
    t = linspace(0,20,500)';
    
    b = 0.01;
    i = 0.09;
    z = 1-b-i;
    
    
    A = 2;
    X = lsode("o_bolniki", [z; b; i], t);
    X = [X(:,1) sum(X(:,1:2),2), ones(size(X(:,1)))];
    save bolniki_1_1.dat X
    
    A = 5;
    X = lsode("o_bolniki", [z; b; i], t);
    X = [X(:,1) sum(X(:,1:2),2), ones(size(X(:,1)))];
    save bolniki_1_2.dat X
    
    A = 2;
    X = lsode("o_bolniki", [z; b; 0; 0; 0; 0; 0; 0; 0; 0; 0; i], t);
    X = [X(:,1) sum(X(:,2:11),2) X(:,12)];
    X = [X(:,1) sum(X(:,1:2),2), ones(size(X(:,1)))];
    save bolniki_7_1.dat X
    
    A = 5;
    X = lsode("o_bolniki", [z; b; 0; 0; 0; 0; 0; 0;0;0;0; i], t);
    X = [X(:,1) sum(X(:,2:11),2) X(:,12)];
    X = [X(:,1) sum(X(:,1:2),2), ones(size(X(:,1)))];
    save bolniki_7_2.dat X
    
    for f = 1:20
        A = 2;
        X = lsode("o_laser", [0; 0.1*f], t);
        save(["laser_" int2str(f) "_1.dat"], "X")
        A = 5;
        X = lsode("o_laser", [0; 0.1*f], t);
        save(["laser_" int2str(f) "_2.dat"], "X")
    endfor
    
    for l = 1:20
        A = 1;
        X = lsode("o_zajlis", [1; 0.05*l], t) .- 1;
        save(["zajlis_" int2str(l) "_1.dat"], "X")
        A = 0.2;
        X = lsode("o_zajlis", [1; 0.05*l], t) .- 1;
        save(["zajlis_" int2str(l) "_2.dat"], "X")
    endfor
endfunction

function vektorji()
    V = [];
    L = [];
    for x = linspace(0,3,30)
        for y = linspace(0,3,30)
            o = o_zajlis([x;y],0)';
            V = [V; x y o ./ sqrt(sumsq(o)) ./ 20];
            o = o_laser([x;y],0)';
            L = [L; x y o ./ sqrt(sumsq(o)) ./ 20];
        endfor
    endfor
    save zajlis_f.dat V
    save laser_f.dat L
endfunction

function bolezni(mu)
    global A
    A = mu;
    b = 0.01;
    t = linspace(0,20,500)';
    Y = [];
    Z = [];
    for i = 0:0.09:1
        X = lsode("o_bolniki", [1-b-i; b; i], t);
        Y = [Y; i max(X(:,2)), (X(500,3) - i)];
        X = lsode("o_bolniki", [1-b-i; b; 0; 0; 0; 0; 0; 0;0;0;0; i], t);
        X = [X(:,1) sum(X(:,2:11),2) X(:,12)];
        Z = [Z; i max(X(:,2)), (X(500,3) - i)];
    endfor
    
    save(["bolni_1_" int2str(A) ".dat"], "Y");
    save(["bolni_7_" int2str(A) ".dat"], "Z");
endfunction

grafi()
vektorji()  
bolezni(1)
bolezni(2)