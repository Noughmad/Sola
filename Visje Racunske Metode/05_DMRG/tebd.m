source mpa.m
source schmidt.m

function S = didx(s1, s2, d)
  S = (s1-1) * d + s2;
endfunction

function U = dvodelcni(z)
  U = spalloc(4, 4, 6);
  U(1,1) = exp(2*z);
  U(4,4) = exp(2*z);
  U(2,2) = U(3,3) = cosh(2*z);
  U(2,3) = U(3,2) = sinh(2*z);
  U = exp(-z) * U;
endfunction

function propagate(A, j, U)
  BB = {};
  
  d = sqrt(length(U));
  
  Mjm = size(A.B{j,1}, 1);
  Mjp = size(A.B{j+1,1}, 2);
  
  for sj = 1:d
    for sjp = 1:d
      BB{sj, sjp} = zeros(Mjm, Mjp);
      for sj2 = 1:d
        for sjp2 = 1:d
          BB{sj, sjp} += U(didx(sj, sjp, d), didx(sj2, sjp2, d)) * A.B{j, sj2} * A.L{j} * A.B{j+1, sjp2};
        endfor
      endfor
    endfor
  endfor
  
  Q = zeros(Mjm*d, Mjp*d);

  if j == 1
    lm = [1];
  else
    lm = A.L{j-1};
  endif

  for km = 1:Mjm
    for kp = 1:Mjp
      for sj = 1:d
        for sjp = 1:d
          Q(didx(km, sj, d), didx(kp, sjp, d)) = lm(km,km) * BB{sj, sjp}(km, kp);
        endfor
      endfor
    endfor
  endfor
  
  [U, S, V, e] = svd_limited(Q);
  
  A.L{j} = S;
  for s = 1:d
    A.B{j, s} = lm \ U(Mjm*(s-1)+1:Mjm*s,:);
    A.B{j+1, s} = V(Mjp*(s-1)+1:Mjp*s,:)';
  endfor
endfunction

function step(A, j0, U, n)
  for j = j0:2:n-1
    propagate(A, j, U);
  endfor
endfunction

function A = osnovno_stanje_heisenberg(n, beta)
  d = 2;
  N = d^n;
  psi = rand(N, 1);
  
  [A, e] = mpa(psi, d);
  
  P = 100;
  z = beta / P;
  U = dvodelcni(z);
  U2 = dvodelcni(z/2);
  
  step(A, 1, U2, n);
  step(A, 2, U, n);
  for s = 1:P-1
    step(A, 1, U, n);
    step(A, 2, U, n);
  endfor
  step(A, 1, U2, n);
endfunction