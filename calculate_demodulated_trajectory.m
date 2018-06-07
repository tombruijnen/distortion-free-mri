function [time,wf] = calculate_demodulated_trajectory(par,varargin)
% Use Michaels code to get the minimal time-waveform where the k-space
% coordinates are demodulated with the first point. This way I set the
% initial conditions similarly.

if nargin <2 
    varargin{1}=1;
end

time={};
wf={};
for ro=1:par.N
    C=[real(par.ktraj(:,ro,1)) imag(par.ktraj(:,ro,1)) zeros(size(par.ktraj(:,ro,1)))];
    C(:,1)=C(:,1)-C(1,1);C(:,2)=C(:,2)-C(1,2);C(:,3)=C(:,3)-C(1,3);
    [~,time_riv,wf{ro}]=minTimeGradient(C,0,0,0,par.Gmax,par.Smax,par.dt,par.ds);      
    time{ro}=linspace(0,time_riv,size(wf{ro},1));
end

for ro=1:par.N
    subplot(223);plot(time{ro},wf{ro},'LineWidth',2);
    axis([0 time{ro}(end) -1.5*par.Gmax 1.5*par.Gmax]);xlabel('Time [ms]');ylabel('G [cT/m]');
    title('Example of k-space at t=tau');grid on;box on
    set(gca,'FontWeight','bold','LineWidth',2,'FontSize',16)
    pause();hold off
end
% END
end