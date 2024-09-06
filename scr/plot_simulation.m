clc; clear all; close all;
set(0,'defaulttextinterpreter','latex');
set(0,'defaultLegendInterpreter','latex');
set(0,'defaultAxesTickLabelInterpreter','latex');
set(0,'defaultAxesFontSize',20)

folder = 'Simulations_satAbs-GDD/';
pg = param_gain(0.532, 20e3, 57); plot(pg(:,1),pg(:,2));

sp = 1;
for nn = 14

    sig = [];
    t = load ([folder,'sPPLT_N_',num2str(nn),'/time.dat']);
    for r = 9936:9999    
        sr = load ([folder,'sPPLT_N_',num2str(nn),'/output_signal_',num2str(r),'_r.dat']);
        si = load ([folder,'sPPLT_N_',num2str(nn),'/output_signal_',num2str(r),'_i.dat']);
        sig =  [sig;sr + 1i* si];
    end
    trt = t(end)-t(1);
    SIZE = length(sig);

    T = linspace( -32*trt, 32*trt, SIZE );
    F = linspace( -0.5*SIZE/(64*trt), +0.5*SIZE/(64*trt), SIZE );



    % subplot(2,3,sp )
    hold on
    area(pg(:,1),pg(:,2), 'FaceColor', 'r', FaceAlpha=0.2)
    plot( F, abs(ifftshift(ifft(sig))).^2/max(abs(ifftshift(ifft(sig))).^2), 'Color', 'b' )
    xlabel('Frequency (THz)')
    ylabel('Norm. Spec. Dens. (a.u.)')
    box on; grid on;
    xlim([-15,15])
    title(['$N =~$', num2str(nn)])
    sp = sp + 1;
end