function [MtoL]=interpolMtoL(N,xGLC,xGC)


T = toeplitz((-1).^(0:N)); 
T = T(:,1:N);

 A=repmat(sqrt(1-xGC'.^2) , [N+1,1] );
 B=repmat( xGLC , [1,N] );
 C=repmat( xGC' , [N+1,1] );
 
 MtoL=T.*A./(B-C)/N;
 
%  MtoL=[MtoL,zeros(N+1,1)];