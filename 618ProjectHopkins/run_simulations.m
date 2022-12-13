%for testing multiple values of Eeqy for a given flip and pmax
flip = 350;
pmax = 3000;
for g=-0.1E-3:0.3E-3:2.1E-3
    [t, p] = main_Margaret(g,flip, pmax);
    figure()
    plot(t, p)
    title(g)
end
