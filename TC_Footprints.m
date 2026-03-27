clear; clc;

% 1. RUTAS
parentFolder = 'PModel_TC_1980_2025';
outputFolder = 'Footprint_TC_1980_2025';

if ~exist(outputFolder, 'dir'), mkdir(outputFolder); end

% 2. OBTENER LISTA DE HURACANES
stormDirs = dir(parentFolder);
stormDirs = stormDirs([stormDirs.isdir] & ~startsWith({stormDirs.name}, '.'));

%% -------------------------------------------------
% 3️. MALLA FIJA
%% -------------------------------------------------
res = 0.1; 
lon_vec = -100:res:-10;
lat_vec = 0:res:50;
[LonGrid, LatGrid] = meshgrid(lon_vec, lat_vec);

fprintf('Iniciando/Reanudando procesamiento de %d huracanes...\n', length(stormDirs));

%% -------------------------------------------------
% 4️. BUCLE PRINCIPAL
%% -------------------------------------------------
for i = 1:length(stormDirs)
    sid = stormDirs(i).name;
    outputFileName = fullfile(outputFolder, [sid '_Hs_max.csv']);
    
    if exist(outputFileName, 'file')
        fprintf('-> Saltando %s (ya existe)\n', sid);
        continue; 
    end
    
    % Extraer año real del nombre de la carpeta (primeros 4 dígitos después del nombre)
    % Ej: ALBERTO_1994... -> extrae 1994
    yearStr = regexp(sid, '\d{4}', 'match', 'once');
    
    currentPath = fullfile(parentFolder, sid);
    csvFiles = dir(fullfile(currentPath, '*.csv'));
    if isempty(csvFiles), continue; end
    
    % Ordenar archivos (P1, P2...)
    [~, idx] = sort(cellfun(@(x) str2double(regexp(x, '\d+', 'match', 'once')), {csvFiles.name}));
    csvFiles = csvFiles(idx);
    
    Hs_max_grid = zeros(size(LonGrid));
    Time_max_grid = strings(size(LonGrid));
    
    try
        for k = 1:length(csvFiles)
            filePath = fullfile(currentPath, csvFiles(k).name);
            
            % --- SOLUCIÓN A LOS WARNINGS ---
            % Detectamos opciones y forzamos que las fechas se lean como TEXTO (string)
            opts = detectImportOptions(filePath);
            posiblesTiempos = {'ISO_Time', 'Time', 'Date', 'tiempo'};
            colTime = posiblesTiempos{ismember(posiblesTiempos, opts.VariableNames)};
            opts = setvartype(opts, colTime, 'string'); 
            
            T = readtable(filePath, opts);
            
            % Captura y limpieza de tiempo
            t_str = T.(colTime)(1);
            
            % CORRECCIÓN DE AÑO (Si sale 0094 -> 1994)
            if contains(t_str, '/00')
                t_str = strrep(t_str, '/00', ['/' yearStr(1:2)]);
            end
            
            lat = T.Latitud; lon = T.Longitud; hs = T.Hs_m;
            hs(isnan(hs)) = 0;
            lon(lon > 180) = lon(lon > 180) - 360;
            
            F = scatteredInterpolant(double(lon), double(lat), double(hs), 'linear', 'none');
            Hs_interp = F(LonGrid, LatGrid);
            Hs_interp(isnan(Hs_interp)) = 0;
            
            mask = Hs_interp > Hs_max_grid;
            Hs_max_grid(mask) = Hs_interp(mask);
            Time_max_grid(mask) = t_str;
        end
        
        % EXPORTACIÓN ORDENADA
        mask_final = Hs_max_grid > 0.1;
        if any(mask_final(:))
            T_out = table(Time_max_grid(mask_final), ...
                          LatGrid(mask_final), ...
                          LonGrid(mask_final), ...
                          Hs_max_grid(mask_final), ...
                          'VariableNames', {'Tiempo', 'Latitud', 'Longitud', 'Hs_max'});
            
            % Ordenar por fecha real antes de guardar
            try
                % Intentamos convertir a fecha para ordenar correctamente
                fechas = datetime(T_out.Tiempo, 'InputFormat', 'dd/MM/yyyy HH:mm');
                [~, sIdx] = sort(fechas);
                T_out = T_out(sIdx, :);
            catch
                % Si falla, se queda con el orden de los archivos
            end
            
            writetable(T_out, outputFileName);
            fprintf('✓ Finalizado: %s\n', sid);
        end
        
    catch ME
        fprintf('❌ Error en %s: %s\n', sid, ME.message);
    end
end
disp('==================================================');
disp('PROCESO COMPLETADO SIN WARNINGS');
