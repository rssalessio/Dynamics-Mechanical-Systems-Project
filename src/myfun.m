
function F = myfun(x,xdata)
%F = x(1)+x(2)*xdata+x(3)*xdata.^2+x(4)*xdata.^3+x(5)*xdata.^4;
F= x(1)./xdata +x(2)./xdata.^2;
end