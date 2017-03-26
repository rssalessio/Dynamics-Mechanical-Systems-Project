sisma=load('sismaspost.txt'); 
 
t=sisma(:,1);
yo1=sisma(:,2);
yo2=sisma(:,3);
figure
subplot(211)
plot(t,yo1);grid on
xlabel('[s]')
ylabel('Y_O_1 [m]')
subplot(212)
plot(t,yo2);grid on
xlabel('[s]')
ylabel('Y_O_2 [m]')