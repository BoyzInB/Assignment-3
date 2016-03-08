%{
%%-----Omega-----
w = [1.7 1.75 1.8 1.85 1.9];
t = [3.595 2.983 2.386 1.808 1.25492];

plot(w,t)
title('100x100 grid with different omega')
xlabel('Omega')
ylabel('Time [s]')
%}

%% Speedup

thr = [1 2 3 4 6 8 10 12];
t = [60.42 40.19 38.6 37.12 42.61 50.88 49.18 51.12];
spup = t(1)./t;

plot(thr, spup)
title('Speedup (200x200 grid)')
xlabel('Number of threads')
ylabel('T(1)/T')

