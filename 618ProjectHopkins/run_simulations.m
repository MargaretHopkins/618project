for g=-0.1E-3:0.3E-3:2.1E-3
    [t, p] = main_Margaret(g);
    figure()
    plot(t, p)
    title(g)
end