function ktraj = dspi_center_out_radial_trajectory(par,varargin)
%% Design through time radial k-space trajectory and corresponding trajectory
% |o----|  |-----|       |----o|    CENTER OUT RADIAL TRAJECTORY
% |-o---|  |-----|       |---o-|        - N time-points
% |--o--|  |--o--|  ...  |--o--|        - M reconstruction matrix
% |-----|  |-o---|       |-----|        - A radial arms
% |-----|  |o----|       |-----|        - Phi radial angle between arms
%
% Each radial arm has N samples per definition.
% Must have an integer number A 

% Simple numeric method to create our radial/spiral
krt=krt_matrix(par.N);

% Define 2D kt-space trajectory in dimensions: [sample ]
ktraj=[];
kx=linspace(0,par.maxk,par.N);
for a=1:par.A
    for t=1:par.N
        for r=1:par.N
            ktraj(r,t,a)=kx(krt(r,t))*exp(-1j*(par.ga1*(r-1)+par.ga2*(a-1)));
        end
    end
end

if nargin > 1
    figure,subplot(221)
    scatter(real(ktraj(:,varargin{1},1)),imag(ktraj(:,varargin{1},1)),50,par.cco,'filled')
    axis(1.5*[-par.maxk par.maxk -par.maxk par.maxk]);xlabel('KX [1/cm]');ylabel('KY [1/cm]')
    title('Example of single readout');grid on;box on
    set(gca,'FontWeight','bold','LineWidth',2,'FontSize',16)
    subplot(222);
    scatter(real(ktraj(3,:,1)),imag(ktraj(3,:,1)),50,par.cco,'filled')
    axis(1.5*[-par.maxk par.maxk -par.maxk par.maxk]);xlabel('KX [1/cm]');ylabel('KY [1/cm]');
    title('Example of k-space at t=tau');grid on;box on
    set(gca,'FontWeight','bold','LineWidth',2,'FontSize',16)
end
% END
end