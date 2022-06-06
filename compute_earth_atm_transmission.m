%% Read Soruce file


%%

Source_file = '/Users/anbu8374/Documents/LibRadTran/libRadtran-2.0.4/data/solar_flux/kurudz_1.0nm.dat';

fileID = fopen(Source_file);
delimiter = ' ';
format_spec = '%f';
N_columns = 2;

A = textscan(fileID, '%d %f', 'Delimiter',' ');


B = importdata(Source_file);
% wl = A.Var3(3:end);
% flux = A.Var4(3:end);

