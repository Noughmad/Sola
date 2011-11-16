global gp;

function y = farm(x, k)
  # [y0, a, p] = k;
  y0 = k(1);
  a = k(2);
  p = k(3);
  y = y0 * x.^p ./ (x.^p + a^p);
endfunction

function y = farm2(x,k)
  # [y0, a] = k;
  y0 = k(1);
  a = k(2);
  y = k(1) * x ./ (x + a);
endfunction 

function y = farm_p(x, k)
  global gp
  y = k(1) * x.^gp ./ (x.^gp + k(2).^gp);
endfunction

function y = ledv1(t, p)
  # [A, lambda] = p;
  y = p(1)*exp(-p(2)*t);
endfunction

function y = ledv2(t, p)
  # [A, B, l1, l2] = p;
  y = p(1)*exp(-p(3).*t) + p(2)*exp(-p(4).*t);
endfunction

function y = ledvk(t, p)
  # [A, lambda, t0] = p;
  y = p(1)*exp(-p(2)*sqrt(abs(t-p(3))));
endfunction

function I = kor_u(U, P)
  I = P(1) * ( exp((U-P(4))/P(2)) - exp(-(U-P(4))/P(3)) );
endfunction

function I = kor(U,P)
  I = kor_u(U, [P; 0]);
endfunction

function farmacija()
  O = dlmread("receptorji.dat");
  x = O(:,1);
  y = O(:,2);
  dy = 3 * ones(size(y));

  [f,p,cvg,iter,corp,covp,covr,stdresid,Z,r2] = leasqr(x, y, [100,24,1], "farm", 0.0001, 20, dy .^-1);
  K = [p, sqrt(diag(covp))]
  chi2red = sumsq((f-y)./dy)/(max(size(y)) - 3)

  [f,p,cvg,iter,corp,covp,covr,stdresid,Z,r2] = leasqr(x, y, [100,24], "farm2", 0.0001, 20, dy .^-1);
  K = [p, sqrt(diag(covp))]
  chi2red = sumsq((f-y)./dy)/(max(size(y)) - 2)
endfunction

function farmacija_hi()
  O = dlmread("receptorji.dat");
  x = O(:,1);
  y = O(:,2);
  R = [];
  dy = 3 * ones(size(y));

  global gp;
  for gp = linspace(0.1,4)
    [f,p,cvg,iter,corp,covp,covr,stdresid,Z,r2] = leasqr(x, y, [100,24], "farm_p", 0.0001, 20, dy .^-1 );
    if (cvg > 0)
      R = [R; gp sumsq((f-y)./ dy)/(max(size(y)) - 2)];
    endif
  endfor

  save farma_hi.dat R
endfunction

function ledvice()
  K = dlmread("ledvice.dat");
  dK = sqrt(K);
  T = (0:80:80*27)';

  N = max(size(K));
  wt = dK .^ -1;

  [f,p,cvg,iter,corp,covp,covr,stdresid,Z,r2] = leasqr(T, K, [K(1), 1/400], "ledv1", 0.0001, 20, wt);
  parametri = [p, sqrt(diag(covp))]
  chi2red = sumsq((f-K)./dK)/(N-2)

  [f,p,cvg,iter,corp,covp,covr,stdresid,Z,r2] = leasqr(T, K, [K(1)/2, K(1)/2, 1/400, 1/4000], "ledv2", 0.0001, 50, wt);
  parametri = [p, sqrt(diag(covp))]
  chi2red = sumsq((f-K)./dK)/(N-4)

  [f,p,cvg,iter,corp,covp,covr,stdresid,Z,r2] = leasqr(T, K, [K(1), 1/30, -40], "ledvk", 0.00001, 20, wt);
  parametri = [p, sqrt(diag(covp))]
  chi2red = sumsq((f-K)./dK)/(N-3)
endfunction

function korozija()
  K = dlmread("korozija.txt");
  U = K(:,1);
  I = K(:,2);
  N = max(size(U));

  [f,p,cvg,iter,corp,covp,covr,stdresid,Z,r2] = leasqr(U, I, [1e-3, 100, 100], "kor");
  parametri = [p, sqrt(diag(covp))]
  chi2red = sumsq((f-I)*1e5)/(N-3)

  [f,p,cvg,iter,corp,covp,covr,stdresid,Z,r2] = leasqr(U, I, [1e-3, 100, 100, 1], "kor_u");
  parametri = [p, sqrt(diag(covp))]
  chi2red = sumsq((f-I)*1e5)/(N-4)
endfunction

# farmacija()
# farmacija_hi()
# korozija()
ledvice()