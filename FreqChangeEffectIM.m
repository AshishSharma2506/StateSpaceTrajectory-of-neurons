%% plotting change in speed-torque char with change in frequency
 
%inputting the motor parameters used in one of the lectures
V_rated = 230/sqrt(3);  % since we plot the per phase char of the IM, rated per phase Voltage 
poles = 4;
f = [20,30,40,50,60,70,80,90]; %fundamental frequency of Supply voltage is 50Hz
VfRatio = V_rated/50; %we have to keep it constant for the Constant Torque period so magnetization I is low

%parameters from manufacturer or from No-load abd Block rotor tes
R1 = 0.095;
R2 = 0.23;
X1 = 0.68;
X2 = 0.672;
Xm = 18.7;

%since impedances depend on frequency,using Fundamental 50Hz value to get Impedances 
L1 = X1/(2*pi*50);
L2 = X2/(2*pi*50);
Lm = Xm/(2*pi*50);

%finding reactances at each freq and creating a vector for each
    x1 = 2*pi*f*L1;
    x2 = 2*pi*f*L2;
    xm = 2*pi*f*Lm;
    

Ns = 120*f/poles; %synchronous speed
omegaS = 4*pi*f/poles; %synchronous speed in rad/sec


%% finding the applied Voltage in each curve according to the frequency

%for f<fundamentalFreq we will decrease V such that V/f is constant
%for f>fundamentalFreq  we will keep the V=constant due to insulation limitation

Vapp = zeros(1,8);

%creating applied voltgae in each case acc to frequency
for i=1:8
   if f(i)<50
       Vapp(i) = VfRatio*f(i);
   
   else
       Vapp(i) = V_rated;
   end
end

%% finding thevenin equivalent circuit values (refer to Induction Motor notes KV5-2)
%thevinin Impedance
Veq = zeros(1,8);
Zth = zeros(1,8);
Zeq = zeros(1,8);
Req = zeros(1,8);
Xeq = zeros(1,8);

for i=1:8
 
    Zth(i) = (1i*xm(i)*(R1+1i*x1(i)))/(R1 + 1i*(x1(i) +xm(i)));
    Zeq(i) = Zth(i);
    Req(i) = real(Zeq(i));
    Xeq(i) = imag(Zeq(i));

    Veq(i) = (Vapp(i)*1i*xm(i))/(R1 + 1i*(x1(i) + xm(i)));
end

%% plotting the Torque Slip char
for j=1:8
    %initialize Slip and Torque vectors
    s = zeros(1,400);
    T = zeros(1,400);
    Speed = zeros(1,400);

    for i=400:-1:1
       k = i-100; 
       s(i)=k/200; %present  value of Slip
       Speed(i) = Ns(j)*(1-s(i)); % rotor speed 

       %calculate Torque
       T(i) = (3/omegaS(j))*(R2/s(i))*(abs(Veq(j))^2/((R2/s(i) + Req(j))^2  + (Xeq(j) + x2(j))^2));

    end
    % plot the speed Torque char
    figure(1)
    hold on
    plot(Speed,T,'linew',2)
    

end

title('Speed-Torque Characteristics'),xlabel('Speed(\omega) - RPM'),ylabel('Torque(\tau) - Nm');
plot(get(gca,'xlim'),[0 0],'k--'),plot([Ns(4) Ns(4)],get(gca,'ylim'),'k--'),plot([0 0],get(gca,'ylim'),'k--');
legend('20Hz','30Hz','40Hz', 'Fundamental Curve(50Hz)','60Hz','70Hz','80Hz','90Hz');
