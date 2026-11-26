% RF_Track reads MADX file
% Objectives
% 1. Check it is matched
% 2. Check IBS rates

RF_Track;

% TAKING SOME PARAMETERS FROM THE MADX
%% Define Twiss parameters
B = load("initial_bunch_25pC.txt");
B0 = Bunch6d(B);
% LATTICE
% OPTION 1
RING = Lattice('../lattice/twiss_withCavity2.tws');

A = RING{'*'};

emitt_x = [];
emitt_y = [];
emitt_z = [];
sigma_pt = [];
T = [];

num_turns = 1000;
for i=1:num_turns
    emitt_z(end+1,1) = B0.get_info().sigma_t;
    emitt_x(end+1,1) = B0.get_info().emitt_x;
    emitt_y(end+1,1) = B0.get_info().emitt_y;
    sigma_pt(end+1,1) = B0.get_info().sigma_pt;
    B1 = RING.track(B0);
    disp(i)
    disp(B1.get_info().emitt_x)
    %T = [ T ; RING.get_transport_table(['%S %beta_x %beta_y %emitt_x %emitt_y %sigma_x %sigma_y %emitt_4d']) ];
    filename = sprintf('../Results/ps_noIBS/phase_space_%d.dat', i);
    M = B1.get_phase_space;
    save(filename, 'M', '-ascii');
    B0 = B1;
end

save -ascii '../Results/noIBS_emittx.dat' emitt_x
%save -ascii '../Results/convergence_tests/mesh/32.dat' emitt_x

