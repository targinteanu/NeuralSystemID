function rgb = colorwheel(theta)

theta = .75*theta;

r = sin(2*pi*theta) + 1;
g = sin(2*pi*theta - 2*pi/3) + 1;
b = sin(2*pi*theta - 4*pi/3) + 1;

r = min(r,1);
g = min(g,1);
b = min(b,1);

rgb = [r,g,b]*.7;

end