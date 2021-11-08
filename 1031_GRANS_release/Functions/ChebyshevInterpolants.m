% ======================================================================= %
% acquaria function
% M.A. Habisreutinger 2013
% ======================================================================= %
 
 
function PHI = ChebyshevInterpolants(np,L,xq)
  % M.A. Habisreutinger 2011
  
  nq   = length(xq) ;
      
  PHI  = zeros(nq,np) ; PHI(:,1)  = 1 ;
  
  if np > 1
    PHI(:,2) = 2*xq/L ;
  end
  
  for ip = 3 : np
    PHI(:,ip)   = 2*PHI(:,2).*PHI(:,ip-1)-PHI(:,ip-2) ;
  end
    
  
end