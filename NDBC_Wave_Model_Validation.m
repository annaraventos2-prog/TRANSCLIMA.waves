clc; 

%% 1. CARGAR DATOS DE ENTRADA 
if ~exist('Lista_Boyas_2010_2025.csv', 'file') 
error('No se encuentra "Lista_Boyas_2010_2025.csv".'); 
end 
B = readtable('Lista_Boyas_2010_2025.csv'); 
BuoyID_List = B.BuoyID; 
LatB_List = B.Lat; 
LonB_List = B.Lon; 
%% 2. CONFIGURACIÓN Y RUTAS 
folder_model = 'PModel_Huracanes_4'; 
folder_ndbc = 'NDBC_Data'; 
kmdeg = 111.1; 
max_dist_km = 1; 
max_time_diff = hours(1.5); 
d = dir(fullfile(folder_model, '**', '*.csv')); 
allfiles = {d(~[d.isdir] & ~contains({d.name}, 'Lista_Boyas')).name}; 
allpaths = {d(~[d.isdir] & ~contains({d.name}, 'Lista_Boyas')).folder}; 
nf = numel(allfiles); 
if nf == 0 
error('No se encontraron snapshots en %s', folder_model); 
end 
fprintf('Total snapshots detectados: %d\n', nf); 
%% 3. PROCESAMIENTO PARALELO 
q = parallel.pool.DataQueue; 
afterEach(q, @nUpdateProgress); 
count_prog = 0; 
results_all = cell(nf,1); 
tic; 
parfor i = 1:nf 
local_res = {}; 
try 
% --- LECTURA ROBUSTA --- 
opts = detectImportOptions(fullfile(allpaths{i}, allfiles{i})); 
opts.Delimiter = ','; 
% Forzamos ISO_Time a 'char' para que no explote con el año 0002 
opts = setvartype(opts, 'ISO_Time', 'char'); 
M = readtable(fullfile(allpaths{i}, allfiles{i}), opts); 
if isempty(M), send(q, i); continue; end 
% --- CAMBIO CRÍTICO 2: CORRECCIÓN DE FECHA MANUAL --- 
% Convertimos el texto a datetime ignorando el error de año 0002 inicialmente 
t_raw = datetime(M.ISO_Time{1}, 'InputFormat', 'dd/MM/yyyy HH:mm'); 
% Si el año es 0002 (o < 100), lo movemos a la década de 2000 
if year(t_raw) < 100 
t_raw.Year = year(t_raw) + 2000; 
end 
t_mod = t_raw; 
yr = year(t_mod); 
% Coordenadas y Hs del modelo 
latM = M.Latitud; 
lonM = M.Longitud; 
hsM = M.Hs_m; 
% Centro del snapshot 
lat0 = mean(latM); lon0 = mean(lonM); 
% Filtrar boyas cercanas (coarse) 
dg = hypot(LonB_List - lon0, LatB_List - lat0); 
idx_boyas = find(dg <= 6); 
if isempty(idx_boyas), send(q, i); continue; end 
for k = 1:numel(idx_boyas) 
id_actual = BuoyID_List(idx_boyas(k)); 
fname = fullfile(folder_ndbc, num2str(id_actual), sprintf('%d_%d.txt', id_actual, yr)); 
if exist(fname, 'file') 
D = read_ndbc_80s(fname); 
if isempty(D), continue; end 
% CRITERIO TEMPORAL 
dt = abs(D.time - t_mod); 
[dtmin, ii] = min(dt); 
if dtmin <= max_time_diff 
% CRITERIO ESPACIAL 
xB = (LonB_List(idx_boyas(k))-lon0)*cosd(lat0)*kmdeg; 
yB = (LatB_List(idx_boyas(k))-lat0)*kmdeg; 
xM = (lonM - lon0).*cosd(lat0)*kmdeg; 
yM = (latM - lat0)*kmdeg; 
[idx_nn, d_real_km] = dsearchn([xM yM], [xB yB]); 
if d_real_km <= max_dist_km 
hs_obs = D.raw(ii, D.col); 
% Evitar valores erróneos de boya (99.0 o NaN) 
if hs_obs < 30 && ~isnan(hs_obs) 
hs_mod = hsM(idx_nn); 
% Solo guardar si el modelo no es NaN en ese punto 
if ~isnan(hs_mod) 
local_res(end+1,:) = {t_mod, allfiles{i}, id_actual, d_real_km, hs_mod, hs_obs}; 
end 
end 
end 
end 
end 
end 
catch ME 
% fprintf('Error en archivo %s: %s\n', allfiles{i}, ME.message); 
end 
results_all{i} = local_res; 
send(q, i); 
end 
%% 4. RESULTADOS FINALES 
results = vertcat(results_all{:}); 
if isempty(results) 
fprintf('\n!!! 0 COINCIDENCIAS !!!\n'); 
else 
T = cell2table(results, 'VariableNames', {'Tiempo','Snapshot','Boya','Dist_km','Modelo','Real'}); 
T.Error = T.Modelo - T.Real; 
writetable(T, 'Validacion_Final_Criterio_1.5h_4.csv'); 
fprintf('\nEXITO: %d coincidencias encontradas.\n', height(T)); 
% Gráfico 
figure('Color','w'); 
scatter(T.Real, T.Modelo, 'filled', 'MarkerFaceAlpha', 0.5); hold on; 
max_val = max([max(T.Real), max(T.Modelo)]); 
line([0 max_val], [0 max_val], 'Color', 'r', 'LineWidth', 2); 
xlabel('Hs Boya (m)'); ylabel('Hs Modelo (m)'); 
title('Validación Modelo (2010-2025)'); grid on; axis equal; 
end 
function nUpdateProgress(~) 
count_prog = count_prog + 1; 
if mod(count_prog, 200) == 0 || count_prog == nf 
fprintf(' -> Avance: %.1f%%\n', (count_prog/nf)*100); 
end 
end 
end 
function D = read_ndbc_80s(fname) 
D = []; 
try 
% Leer matriz de datos 
data = readmatrix(fname, 'NumHeaderLines', 1); 
if isempty(data), return; end 
% Leer cabecera para encontrar WVHT 
fid = fopen(fname,'r'); 
h = fgetl(fid); fclose(fid); 
cols = strsplit(strtrim(h)); 
c_idx = find(strcmpi(cols, 'WVHT') | strcmpi(cols, 'Hs'), 1); 
if isempty(c_idx), return; end 
% --- CAMBIO CRÍTICO 3: CORRECCIÓN AÑOS BOYA --- 
yy_raw = data(:,1); 
yy = yy_raw; 
% Lógica robusta: años < 70 son 2000s, años >= 70 y < 100 son 1900s 
yy(yy_raw < 70) = yy_raw(yy_raw < 70) + 2000; 
yy(yy_raw >= 70 & yy_raw < 100) = yy_raw(yy_raw >= 70 & yy_raw < 100) + 1900; 
mm = data(:,2); dd = data(:,3); hh = data(:,4); 
D.time = datetime(yy, mm, dd, hh, 0, 0); 
D.raw = data; 
D.col = c_idx; 
catch 
return; 
end
