% RF_Track reads MADX file
% Objectives
% 1. Check it is matched
% 2. Check IBS rates

RF_Track;
RF_Track.number_of_threads = 8;
IB = IntraBeamScattering(32,32,32); % Sin Ncol

% TAKING SOME PARAMETERS FROM THE MADX
%% Define Twiss parameters
B = load("initial_bunch_25pC.txt");
B0 = Bunch6d([B; B; B; B; B]);
% LATTICE
% OPTION 1
RING = Lattice('../lattice/twiss_withCavity2.tws');

A = RING{'*'};
%RF = RING.get_rf_elements();
%RF{1}.set_phid(90.0);
%for i = 1:numel(A)
%    L = A{i}.get_length();   % call method once
%    if L == 0.0
%        % do nothing
%    elseif L < 1.0
%        A{i}.add_collective_effect(IB);
%        A{i}.set_cfx_nsteps(1);
%
%    else
%        A{i}.add_collective_effect(IB);
%        A{i}.set_cfx_nsteps(2);
%    end
%end
RING.add_collective_effect(IB);
RING.set_cfx_nsteps(10);

emitt_x = [];
emitt_y = [];
emitt_z = [];
sigma_pt = [];
T = [];

num_turns = 5000;
for i=1:num_turns
    emitt_z(end+1,1) = B0.get_info().sigma_t;
    emitt_x(end+1,1) = B0.get_info().emitt_x;
    emitt_y(end+1,1) = B0.get_info().emitt_y;
    sigma_pt(end+1,1) = B0.get_info().sigma_pt;
    B1 = RING.track(B0);
    disp(i)
    disp(B1.get_info().emitt_x)
    %T = [ T ; RING.get_transport_table(['%S %beta_x %beta_y %emitt_x %emitt_y %sigma_x %sigma_y %emitt_4d']) ];
    filename = sprintf('../Results/ps/phase_space_%d.dat', i);
    M = B1.get_phase_space;
    fid = fopen(filename, 'w');
    fprintf(fid, '%.12g ', M.');
    fprintf(fid, '\n');
    fclose(fid);
    B0 = B1;
    save -ascii '../Results/RFT_25pC_nemit_x_IBS_1000.dat' emitt_x
    save -ascii '../Results/RFT_25pC_nemit_y_IBS_1000.dat' emitt_y
end

%save -ascii '../Results/convergence_tests/Ncol/10.dat' emitt_x
%save -ascii '../Results/convergence_tests/mesh/32.dat' emitt_x

