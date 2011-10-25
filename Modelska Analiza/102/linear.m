function semafor_lin(N, vmax, amax, amin, y0)
  c = linspace(0,1,N)';

  O = ones(N,1);

  Aa = diag(ones(N,1)) - diag(ones(N-1,1),-1);
  bvu = amax * O;
  bvu(1) = bvu(1) + y0;
  bvl = amin * O;
  bvl(1) = bvl(1) + y0;
  
  As = ones(1,N);
  bs = N;

  A = [ Aa; Aa; As ];
  b = [ bvu; bvl; bs ];

  ctype = [];
  for i = 1:N
    ctype = [ctype "U"];
  end

  for i = 1:N
    ctype = [ctype "L"];
  end

  ctype = [ctype "S"];

  vartype = [];
  for i = 1:N
    vartype = [ vartype "C" ];
  endfor

  [X, F, status, extra] = glpk(c, A, b, zeros(100,1), vmax * ones(100,1), ctype, vartype, -1);

  status
  F
  c' * X

  saveTo = [ "lin_" int2str(10*y0) ".dat" ]
  save(saveTo, "X");
end

semafor_lin(100, 2, 0.02, -0.05, 1.5)