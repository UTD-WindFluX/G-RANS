function [xGC,D1GC]=nodichebGC(N)

xGC = cos(pi*(2*(0:N-1)+1)/(2*N))';

[xGLC,D1GLC]=nodiD1chebGLC(N);
[Eder]=interderivaEtoL(N,xGLC,xGC);

% % [MtoL]=interpolanodale(N,xGC,1);
[MtoL]=interpolMtoL(N,xGLC,xGC);

D1GC=(MtoL\Eder);