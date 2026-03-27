clear; clc;

% ===============================================================
% CONFIGURACION
% ===============================================================

archivo_csv = 'IBTrACS_NA_1980_2025.csv';

opts = detectImportOptions(archivo_csv);
opts.VariableNamingRule = 'preserve';
datos_totales = readtable(archivo_csv,opts);
datos_totales.SID = categorical(datos_totales.SID);

carpeta_maestra = 'PModel_Era5';
if ~exist(carpeta_maestra,'dir'), mkdir(carpeta_maestra); end

[lista_sid,idx_unicos] = unique(datos_totales.SID,'stable');
nombres_asociados = datos_totales.NAME(idx_unicos);

fprintf('--- INICIANDO %d HURACANES ---\n',numel(lista_sid));

p = gcp('nocreate');
if isempty(p)
    parpool('local',3);
end


% ===============================================================
% LOOP HURACAN
% ===============================================================

for h = 1:numel(lista_sid)

    sid_actual    = lista_sid(h);
    nombre_actual = nombres_asociados{h};

    id_huracan = regexprep( ...
        sprintf('%s_%s',string(nombre_actual),string(sid_actual)), ...
        '[^\w-]','_');

    ruta_subcarpeta = fullfile(carpeta_maestra,id_huracan);
    
    % SALTAR SI YA EXISTE:
    if exist(ruta_subcarpeta,'dir')
        fprintf('Saltando %s (ya procesado)\n', id_huracan);
        continue; 
    end
    
    mkdir(ruta_subcarpeta);


    data = datos_totales(datos_totales.SID==sid_actual,:);
    n = height(data);

    fprintf('Procesando %s (%d pasos)\n',id_huracan,n-1);

    % ---- extracción a arrays (MUY importante)
    LAT  = data.LAT;
    LON  = data.LON;
    WIND = data.USA_WIND;
    RMW  = data.USA_RMW;
    ISO  = data.ISO_TIME;

    Tout = cell(n-1,1);

    parfor i = 1:n-1

        Lat0 = LAT(i);
        Lon0 = LON(i);
        Lat1 = LAT(i+1);
        Lon1 = LON(i+1);

        Vmax = WIND(i)*0.514444;

        iso_time_val = ISO(i);
        if iscell(iso_time_val), iso_time_val = iso_time_val{1}; end

        dist_km = haversine(Lat0,Lon0,Lat1,Lon1);
        Vfm = (dist_km*1000)/(3*3600);

        rmw = RMW(i);

if isnumeric(rmw) && isscalar(rmw) && ~isnan(rmw) && rmw > 0
    Rmax = rmw * 1852;
else
    Rmax = 30000;
end

        R34 = ((2.5*(Rmax/1000))+124)*1000;

        try
            Tout{i} = PModel_fast_identical( ...
                Vmax,Vfm,Rmax,R34,Lat0,Lon0,i,ruta_subcarpeta,iso_time_val);
        catch ME
            fprintf('Error %s paso %d: %s\n',id_huracan,i,ME.message);
        end
    end

    % ---- escritura fuera del parfor
    for i=1:n-1
        if ~isempty(Tout{i})
            writetable(Tout{i}, ...
                fullfile(ruta_subcarpeta,Tout{i}.Properties.Description));
        end
    end
end

disp('--- FINALIZADO ---');

% ===============================================================
% MODELO (RAPIDO + IDENTICO)
% ===============================================================

function T = PModel_fast_identical( ...
    Vmax,Vfm,Rmax,R34,Lat0,Lon0,paso,ruta,iso_time)

persistent MATCACHE_ZZ MATCACHE_DATA km_per_deg Req ...
           a b c d e f g h i0 C ...
           Vmax_runs Dp_runs Vfm_runs R34_runs R_runs

if isempty(MATCACHE_ZZ)

    MATCACHE_ZZ   = containers.Map;
    MATCACHE_DATA = containers.Map;

    km_per_deg = 111.1;
    Req = 30;

    a=0.54; b=-169; c=-1442; d=0.3; e=14.3; f=-43;
    g=9600; h=4470; i0=1e5; C=0.1;

    Vmax_runs = [17 30 40 50 65 78];
    Dp_runs   = [10 30 50 70 110 150];
    Vfm_runs  = [0 2.5 5 7.5 10 12.5 15];
    R34_runs  = [200 300 400];
    R_runs    = [15 30 60];
end

% ================= LIMITES =================
Vmax = min(max(Vmax,17),78);
Vfm  = min(max(Vfm,0),15);
Rmax = min(max(Rmax,15e3),60e3);
R34  = min(max(R34,200e3),400e3);

% ================= FISICA ==================
F_P32 = (a*Vmax^3 + b*Vmax^2 + c*Vfm^2 + ...
         d*Vmax^2*Vfm + e*Vmax*Vfm^2 + ...
         f*Vmax*Vfm + g*Vmax + h*Vfm + i0)*exp(C*Vfm);

lambda = 0.85*log10(Rmax/30e3)+1;
gamma  = 0.65*log10(R34/300e3)+1;

F = (F_P32*lambda*gamma)/1000;
Hs_max = 0.89*((0.0016*((9.81*F*1000)^0.5)*Vmax)/9.81);

% ================= RUNS ==================
iV  = find_indices(Vmax_runs,Vmax);
iF  = find_indices(Vfm_runs,Vfm);
i34 = find_indices(R34_runs,R34/1000);
iR  = find_indices(R_runs,Rmax/1000);

% ================= CACHE ZZ ==================
keyZZ = sprintf('%d_%d_%d_%d',iV(1),iF(1),i34(1),iR(1));

if isKey(MATCACHE_ZZ,keyZZ)

    Cc = MATCACHE_ZZ(keyZZ);
    ZZ  = Cc.ZZ;
    xrm = Cc.xrm;
    yrm = Cc.yrm;

else
    % --- carga identica a la original
    keyD = sprintf('%d_%d_%d_%d', ...
        Dp_runs(iV(1)),round(Vfm_runs(iF(1))*10), ...
        R34_runs(i34(1)),R_runs(iR(1)));

    if isKey(MATCACHE_DATA,keyD)
        D = MATCACHE_DATA(keyD);
    else
        D.Zc = cell(2,2,2,2);
        for a1=1:2
            for a2=1:2
                for a3=1:2
                    for a4=1:2
                        fname = sprintf('wave_%d_%d_%d_%d_Hs_Hsmax.mat', ...
                            Dp_runs(iV(a1)), ...
                            round(Vfm_runs(iF(a2))*10), ...
                            R34_runs(i34(a3)), ...
                            R_runs(iR(a4)));
                        tmp = load(fullfile('Hsmax_diagrams',fname));
                        D.Zc{a1,a2,a3,a4} = tmp.Z;
                        D.xrm = tmp.xrm;
                        D.yrm = tmp.yrm;
                    end
                end
            end
        end
        MATCACHE_DATA(keyD) = D;
    end

    % --- interpolación punto a punto
    nx = length(D.xrm); ny = length(D.yrm);
    ZZ = zeros(nx,ny);
    for xx=1:nx
        for yy=1:ny
            z = zeros(2,2,2,2);
            for a1=1:2
                for a2=1:2
                    for a3=1:2
                        for a4=1:2
                            z(a1,a2,a3,a4) = D.Zc{a1,a2,a3,a4}(xx,yy);
                        end
                    end
                end
            end
            ZZ(xx,yy) = interpn( ...
                Vmax_runs(iV),Vfm_runs(iF), ...
                R34_runs(i34),R_runs(iR), ...
                z,Vmax,Vfm,R34/1000,Rmax/1000);
        end
    end

    MATCACHE_ZZ(keyZZ) = struct('ZZ',ZZ,'xrm',D.xrm,'yrm',D.yrm);
    xrm = D.xrm; yrm = D.yrm;
end

% ================= GEOGRAFIA ==================
if Lat0<0
    Zfinal = fliplr(ZZ');
else
    Zfinal = ZZ';
end

Hs = Zfinal * Hs_max;

[GX,GY] = meshgrid(xrm*Req,yrm*Req);
LatGrid = Lat0 + GY/km_per_deg;
LonGrid = Lon0 + GX./(km_per_deg*cosd(Lat0));

% ================= TABLA ==================
[~,id_h] = fileparts(ruta);
fname = sprintf('%s_P%03d.csv',id_h,paso);

T = table( ...
    repmat(string(iso_time),numel(Hs),1), ...
    LatGrid(:),LonGrid(:),Hs(:), ...
    'VariableNames',{'ISO_Time','Latitud','Longitud','Hs_m'});

T.Properties.Description = fname;
end

% ===============================================================
% UTILIDADES
% ===============================================================
function idx = find_indices(a,v)
    [~,s] = sort(abs(a-v));
    idx = sort(s(1:2));
end

function d = haversine(lat1,lon1,lat2,lon2)
    R = 6371;
    dLat = deg2rad(lat2-lat1);
    dLon = deg2rad(lon2-lon1);
    a = sin(dLat/2)^2 + ...
        cosd(lat1)*cosd(lat2)*sin(dLon/2)^2;
    d = R*2*atan2(sqrt(a),sqrt(1-a));
end