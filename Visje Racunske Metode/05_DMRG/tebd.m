source mpa.m
source schmidt.m

function S = didx(s1, s2, d)
  S = (s1-1) * d + s2;
endfunction

function U = dvodelcni(z)
  U = zeros(4, 4);
  U(1,1) = U(4,4) = exp(2*z);
  U(2,2) = U(3,3) = cosh(2*z);
  U(2,3) = U(3,2) = sinh(2*z);
  U = sparse(U * exp(-z));
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

function [MA, factor] = mpa_normalize(A)
  psi = mpa_reverse(A);
  factor = norm(psi);
  MA = mpa(psi / norm(psi), A.d);
endfunction

function Norm = osnovno_stanje_heisenberg(n, beta)
  d = 2;
  N = d^n;
  psi = rand(N, 1);
  psi /= norm(psi);
  
  [A, e] = mpa(psi, d);
  
  P = 20;
  B = 20;
  z = -beta / P / B;
  U = dvodelcni(z);
  U2 = dvodelcni(z/2);
  logfactor = 0;
  
  SpinOp = [1, 0; 0, -1];
  Norm = [0, 0, log(mpa_norm(A, d)), exp_value_local(A, SpinOp, n/2)];
  
  for s = 1:P-1
    A = step(A, 1, U2, n);
    A = step(A, 2, U, n);
    for b = 1:B
      A = step(A, 1, U, n);
      A = step(A, 2, U, n);
    endfor
    A = step(A, 1, U2, n);
    Norm = [Norm; s, s*B*z, logfactor + log(mpa_norm(A, d)), exp_value_local(A, SpinOp, n/2)];

    [A, f] = mpa_normalize(A);
    logfactor += log(f);
  endfor
  
  dlmwrite("g_tebd_norm.dat", Norm, " ")
endfunction

function n = mpa_norm(A, d)
  n = 1;

  t = 0;
  for s = 1:d
    t += kron(A.B{1, s}, A.B{1, s});
  endfor
  n = t;

  for j = 2:A.n
    t = 0;
    for s = 1:d
      t += kron(A.L{j-1} * A.B{j,s}, A.L{j-1} * A.B{j,s});
    endfor
    n *= t;
  endfor
  
  n = sqrt(n);
endfunction

function test_mpa_norm(d, n)
  for i = 1:10
    psi = rand(d^n, 1);
    [A, e] = mpa(psi, d);

    assert(abs(Norm - norm(psi)) <= (1.0001 * e + 1e-12));
  endfor
endfunction

function o = exp_value_local(A, O, k)
  d = A.d;

  t = 0;
  for s = 1:d
    t += kron(A.B{1, s}, A.B{1, s});
  endfor
  o = t;

  for j = 2:k-1
    t = 0;
    for s = 1:d
      t += kron(A.L{j-1} * A.B{j,s}, A.L{j-1} * A.B{j,s});
    endfor
    o *= t;
  endfor

  n = o;

  t = 0;
  tn = 0;
  j = k;
  for s = 1:d
    for sp = 1:d
      t += O(s,sp) * kron(A.L{j-1} * A.B{j,s}, A.L{j-1} * A.B{j,s});
    endfor
    tn += kron(A.L{j-1} * A.B{j,s}, A.L{j-1} * A.B{j,s});
  endfor
  o *= t;
  n *= tn;
  
  for j = k+1:A.n
    t = 0;
    for s = 1:d
      t += kron(A.L{j-1} * A.B{j,s}, A.L{j-1} * A.B{j,s});
    endfor
    o *= t;
    n *= t;
  endfor

  o /= n;
endfunction

function o = exp_value_global(A, O, k)
  d = A.d;
  t = 0;
  for s = 1:d
    t += kron(A.B{1, s}, A.B{1, s});
  endfor
  o = t;

  for j = 2:A.n
    t = 0;
    for s = 1:d
      for sp = 1:d
        t += O(s,sp) * kron(A.L{j-1} * A.B{j,s}, A.L{j-1} * A.B{j,s});
      endfor
    endfor
    o *= t;
  endfor
  
  o = sqrt(n);  
endfunction

function casovni_razvoj(n, T)
  d = 2;
  N = d^n;
  psi = zeros(N,1);
  b = "";
  for i=1:n/2
    b = ["0" b "1"];
  endfor
  psi(bin2dec(b)) = 1;

  NormPsi = norm(psi)
  A = mpa(psi, 2);
  NormMpa = mpa_norm(A, 2)

  P = 50;

  clear i;
  z = -i*T / P;

  U = dvodelcni(z);
  U2 = dvodelcni(z/2);
  
  SpinOp = [1, 0; 0, -1];
  Norm = [0, 0, mpa_norm(A, d), exp_value_local(A, SpinOp, n/2)];
  
  for s = 1:P-1
    A = step(A, 1, U2, n);
    A = step(A, 2, U, n);
    A = step(A, 1, U2, n);
    Norm = [Norm; s, abs(s*z), mpa_norm(A, d), exp_value_local(A, SpinOp, n/2)];
  endfor
  
  dlmwrite("g_tebd_time.dat", Norm, " ")

endfunction