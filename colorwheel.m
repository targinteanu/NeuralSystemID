function rgb = colorwheel(theta)

%theta = .75*theta;

r = .5*sin(2*pi*theta) + .5;
g = .5*sin(2*pi*theta - 2*pi/3) + .5;
b = .5*sin(2*pi*theta - 4*pi/3) + .5;

rgb = [r,g,b]*.7;

end