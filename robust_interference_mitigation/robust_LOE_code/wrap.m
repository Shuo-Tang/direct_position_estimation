function x = wrap( x, x_max )

while( sum( abs(x) > x_max ) ~= 0)
    x(abs(x)>x_max)  =   x(abs(x)>x_max) - sign(x(abs(x)>x_max))*2*x_max;
end