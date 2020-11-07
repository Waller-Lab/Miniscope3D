function y = soft(x,T)

y = sign(x).*max(abs(x) - T, 0);
%y = y./(y+T) .* x;


