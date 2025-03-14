clear
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%         INPUT SEQUENCE HERE
%
ergocubknee = load('ergocubknee.mat','data');
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

KmOpt = 0.0;
KwOpt = 0.0;
FcOpt_pos = 0.0;
FsOpt_pos = 0.0;
FcOpt_neg = 0.0;
FsOpt_neg = 0.0;
S0Opt = 0.0;
S1Opt = 0.0;
VtOpt = 0.0;

sz = size(ergocubknee.data)
N = sz(2);

time        = ergocubknee.data(1,:);
current     = ergocubknee.data(2,:);
velocity    = ergocubknee.data(3,:);
torque_meas = ergocubknee.data(4,:);

warmup = N/5;

velmax = max(abs(velocity));
velmax_pos = max(velocity);
velmax_neg = min(velocity);
velthr_pos =  5; %velmax/10
velthr_neg = -5; %-velmax/10 %velmax_neg/10;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STAGE 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

errMin = 1E+12;

KmA = 3.0; KmB = 25.0;
KwA = 0.0; KwB = 0.1;
FcA_pos =  2.0;  FcB_pos = 20.0;
FsA_pos =-10.0;  FsB_pos = 10.0;
FsA_neg =-10.0;  FsB_neg = 10.0;
VtA = 0.01*velmax; VtB = velmax;

for k=1:3

    dKm = 0.1*(KmB - KmA);
    dKw = 0.1*(KwB - KwA);
    dFc_pos = 0.1*(FcB_pos - FcA_pos);
    dFs_pos = 0.1*(FsB_pos - FsA_pos);
    dFs_neg = 0.1*(FsB_neg - FsA_neg);
    dVt = 0.1*(VtB - VtA);

    for Km = KmA : dKm : KmB
        for Kw = KwA : dKw : KwB
            for Fc_pos = FcA_pos : dFc_pos : FcB_pos
                for Fs_pos = FsA_pos : dFs_pos : FsB_pos
                    for Fs_neg = FsA_neg : dFs_neg : FsB_neg
                        for Vt = VtA : dVt : VtB

                            err = 0.0;

                            AVold = abs(velocity(1));

                            for i=1:N

                                %if abs(velocity(i)) < AVold

                                    if velocity(i) > 5

                                        torque_calc = Km*current(i) - Kw*velocity(i);
                                        err = err + abs(torque_calc - torque_meas(i) - Fc_pos - exp(-velocity(i)/Vt)*Fs_pos);
                                        if err > errMin
                                            break
                                        end
                                    elseif velocity(i) < -5

                                        torque_calc = Km*current(i) - Kw*velocity(i);
                                        err = err + abs(torque_calc - torque_meas(i) + Fc_pos + exp(velocity(i)/Vt)*Fs_neg);
                                        if err > errMin
                                            break
                                        end
                                    end
                                %end

                                AVold = abs(velocity(i));

                            end

                            if err < errMin
                                errMin = err;
                                KmOpt = Km;
                                KwOpt = Kw;
                                FcOpt_pos = Fc_pos;
                                FsOpt_pos = Fs_pos;
                                FsOpt_neg = Fs_neg;
                                VtOpt = Vt;
                            end
                        end
                    end
                end
            end
        end
    end

    KmA = KmOpt - 2*dKm; KmB = KmOpt + 2*dKm;
    KwA = KwOpt - 2*dKw; KwB = KwOpt + 2*dKw;
    FcA_pos = FcOpt_pos - 2*dFc_pos; FcB_pos = FcOpt_pos + 2*dFc_pos;
    FsA_pos = FsOpt_pos - 2*dFs_pos; FsB_pos = FsOpt_pos + 2*dFs_pos;
    FsA_neg = FsOpt_neg - 2*dFs_neg; FsB_neg = FsOpt_neg + 2*dFs_neg;
    VtA = VtOpt - 2*dVt; VtB = VtOpt + 2*dVt; 

    if KmA < 0
        Kma = 0;
    end
    if KwA < 0
        KwA = 0;
    end
    if FcA_pos < 0
        FcA_pos = 0;
    end
    if VtA < 0
        VtA = 0;
    end

    FcOpt_neg = FcOpt_pos;

    k
    KmOpt
    KwOpt
    FcOpt_pos
    FcOpt_neg
    FsOpt_pos
    FsOpt_neg
    VtOpt
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STAGE 2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

errMin = 1E+12;

S0A = 0.0; S0B = 25.0;
S1A = -25; S1B = 25.0;

for k = 1:3

dS0 = 0.02*(S0B - S0A);
dS1 = 0.02*(S1B - S1A);

for S0 = S0A : dS0 : S0B
for S1 = S1A : dS1 : S1B

err = 0.0;

Z = 0.0;

time_prev = time(1);

for i=1:N

    gV = 0;
    
    if velocity(i) > 0
        gV = FcOpt_pos + FsOpt_pos*exp(-abs(velocity(i)/VtOpt));
    else
        gV = FcOpt_neg + FsOpt_neg*exp(-abs(velocity(i)/VtOpt));
    end

    Zdot = (velocity(i) - abs(velocity(i))*S0*Z/gV);

    Z = Z + (time(i)-time_prev)*Zdot;
    time_prev = time(i);

    friction_calc = KwOpt*velocity(i) + S0*Z + S1*Zdot;

    torque_motor = KmOpt*current(i);

    if i > warmup
        err = err + ((torque_motor - torque_meas(i)) - friction_calc)^2;
    end

    if err > errMin
        break
    end
end

if err < errMin
    errMin = err;
    S0Opt = S0;
    S1Opt = S1;
end

end
end
S0A = S0Opt-dS0; S0B = S0Opt+dS0;
S1A = S1Opt-dS1; S1B = S1Opt+dS1;

end

% Km = KmOpt
% Kw = KwOpt
% S0 = S0Opt
% S1 = S1Opt
% Fc_pos = FcOpt
% Fc_neg = FcOpt
% Fs_pos = FcOpt+FsOpt
% Fs_neg = FcOpt+FsOpt
% Vth = VtOpt

fprintf('Estimated LuGre parameters:\n\n')
fprintf('Km       %f Nm/A\nKw       %f Nm/(rad/s)\n',KmOpt,KwOpt)
fprintf('S0       %f Nm/rad\nS1       %f Nm/(rad/s)\n',S0Opt,S1Opt)
fprintf('Fc_pos   %f Nm\nFc_neg   %f Nm\n',FcOpt_pos,FcOpt_pos)
fprintf('Fs_pos   %f Nm\nFs_neg   %f Nm\nVth      %f rad/s\n\n',FsOpt_pos,FsOpt_neg,Vt)


