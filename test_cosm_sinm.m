%TEST_COSM_SINM Brief test of new matrix cosine and sine functions.

n = 100;
cosmdef = @(x) (expm(1i*x) + expm(-1i*x))/2;
sinmdef = @(x) (expm(1i*x) - expm(-1i*x))/(2*1i);

% % If the Symbolic Math Toolbox is available uncomment the next three
% % lines in order to compare with "exact" results.
% n = 10;
% cosmdef = @(x) double( (expm(1i*vpa(x)) + expm(-1i*vpa(x)))/2 );
% sinmdef = @(x) double( (expm(1i*vpa(x)) - expm(-1i*vpa(x)))/(2*1i) );

disp('Relative differences between cos and sin computed from the definition') 
disp('with EXPM (complex arithmetic) and from the new algorithms.') 
fprintf('Random A of size %g.\n', n) 

for use_schur = 0:1

  switch use_schur
     case 0, fprintf('\nWith Schur factorization.\n')
     case 1, fprintf('\nWithout Schur factorization.\n')
  end   

  disp('Normwise relative differences.')
  disp(' cos        sin        cossin')
  for k = 1:4
      A = randn(n);
      rescos = norm(cosmdef(A) - cosm(A, use_schur))/norm(cosmdef(A));
      ressin = norm(sinmdef(A) - sinm(A, use_schur))/norm(sinmdef(A));
      [C, S] = cosmsinm(A,use_schur);
      rescs = max(norm(cosmdef(A) - C)/norm(cosmdef(A)),...
                  norm(sinmdef(A) - S)/norm(sinmdef(A)));
      fprintf('%9.2e  %9.2e  %9.2e\n', rescos, ressin, rescs)
  end

end