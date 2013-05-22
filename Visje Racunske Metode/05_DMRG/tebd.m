source mpa.m
source schmidt.m

function S = didx(s1, s2, d)
  S = (s1-1) * d + s2;
endfunction

function U = dvodelcni(z)
  z *= -1;
  U = spalloc(4, 4, 6);
  U(1,1) = exp(2*z);
  U(4,4) = exp(2*z);
  U(2,2) = U(3,3) = cosh(2*z);
  U(2,3) = U(3,2) = sinh(2*z);
  U = exp(-z) * U;
endfunction

function Ap = propagate(A, j, U)
  d = A.d;
  
  Mjm = size(A.B{j,1}, 1);
  Mjp = size(A.B{j+1,1}, 2);
  
  BB = {};
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
  
  Ap = A;
endfunction

function Ap = step(A, j0, U, n)
  for j = j0:2:n-1
    A = propagate(A, j, U);
  endfor
  Ap = A;
endfunction

function A = osnovno_stanje_heisenberg(n, beta)
  d = 2;
  N = d^n;
  psi = rand(N, 1);
  psi /= norm(psi);
  
  [A, e] = mpa(psi, d);
  
  P = 20;
  B = 20;
  z = beta / P;
  U = dvodelcni(z);
  U2 = dvodelcni(z/2);
  
  Norm = [0, mpa_norm(A, d)];
  
  for s = 1:P-1
    A = step(A, 1, U2, n);
    A = step(A, 2, U, n);
    for b = 1:B
      A = step(A, 1, U, n);
      A = step(A, 2, U, n);
    endfor
    A = step(A, 1, U2, n);
    Norm = [Norm; s, mpa_norm(A, d)];
  endfor
  
  dlmwrite("g_tebd_norm.dat", Norm, " ")
endfunction

function n = mpa_norm(A, d)
  n = 1;
  for j = 1:A.n
    t = 0;
    for s = 1:d
      t += kron(A.B{j,s}, A.B{j,s});
    endfor
    n *= t;
  endfor
  assert(size(n) == [1 1]);
endfunction
