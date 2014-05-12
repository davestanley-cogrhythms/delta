
% Traub GABAB Synapse
t = 1:2000;
y = (1-exp(-(t-10)/38.1)).^4 .* (10.2 * exp(-(t-10)/122) + 1.1 * exp(-(t-10)/587));
figure; plot(t,y)
xlabel('Time (ms)')
ylabel('G')
