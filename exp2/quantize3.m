function y = quantize3(x,B)  
% y = quantize3(x,B) 
% The function rounds x into a binary fixed point 
% representation with B bits, including the sign bit. 

% No overflow occurs --- the variables are normalized by the largest value
y = 0;
B1 = B-1;
%s = sign(x);
mx=max(abs(x));
if mx == 0
    y = x;
else
    expon = ceil(log2(max(abs(x))));
    fac = 2^(B1-expon);
    invfac = 2^(-B1+expon);
    xint = round(x*fac);
    y = invfac * xint;
    %y = invfac * xint; % sem overflow
end
%y=x;