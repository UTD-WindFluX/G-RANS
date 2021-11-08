function [Eder]=interderivaEtoL(N,xGLC,xGC)


T = toeplitz((-1).^(0:N)); 
T=T(1:N-1,1:N);

 A=repmat(sqrt(1-xGC'.^2) , [N-1,1] );
 B=repmat( xGLC(2:end-1) , [1,N] );
 C=repmat( xGC' , [N-1,1] );
 
 Eder=T.*A./(B-C).^2/N;
 
 E0=((-1).^(0:N-1)).*sqrt(1-xGC'.^2).*(N^2./(1-xGC')-1./(1-xGC').^2 )/N;
 
 EN=((-1).^(N:2*N-1)).*sqrt(1-xGC'.^2).*(N^2./(1+xGC')-1./(1+xGC').^2 )/N;
 
 Eder=[E0;Eder;EN];
 
%  Eder=[Eder,zeros(N+1,1)];