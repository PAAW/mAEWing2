function regression_check( result, correct, tol )
% Verify that result is correct for sqp/slp_trust regression test
if (size(result)~=size(correct))
   if length(result)==length(correct)
      disp('regression_check: reshaping result and answer column vectors')
      result = result(:);
      correct = correct(:);
   else
      error('regression_check:size','Size result and correct answer differ')
   end
end
if nargin<3, tol=1e-6; end
if any(abs(result-correct)>tol)
   warning('Regression failure for sqp/slp_trust')
   disp(['Max difference = ',num2str(max(abs(result-correct))),...
        ' > tolerance = ',num2str(tol)])
   disp(' ')
   disp('Result vs. Correct Answer')
   disp([result(:), correct(:)])
end
end