close all;
clear all;

runset = dlmread('E:\LSM_GW\GIS\HCDN_Sites_Plasticity_Estimate_matlab.csv');

%Read in precip spatial data
inpath = strcat('E:\LSM_GW\LSM_JoFlo\Precip_Livneh\');
infile = strcat(inpath,strcat('prec.1999.nc'));
precip_lat = double(ncread(infile,'lat'));
precip_lon = double(ncread(infile,'lon'));

%Read in air temp spatial data
inpath = strcat('E:\LSM_GW\LSM_JoFlo\Temp_Livneh\');
infile = strcat(inpath,strcat('tmax.1999.nc'));
test = ncinfo(infile);
temp_lat = double(ncread(infile,'lat'));
temp_lon = double(ncread(infile,'lon'));

plotcounter = 1;
for site = 1:length(runset(:,1))
    
    disp(num2str(runset(site,1)));

    %Extract precip TS
    sitelat = runset(site,3);
    sitelon = runset(site,4);
    [minValue,lat_index] = min(abs(transpose(precip_lat) - sitelat));
	[minValue,lon_index] = min(abs(transpose(precip_lon) - sitelon));
    if sitelon < 0
        [minValue,lon_index] = min(abs(transpose(precip_lon) - (360+sitelon)));
    end
    Precip_All = zeros(1,1);
    for year = 1999:2010
        inpath = strcat('E:\LSM_GW\LSM_JoFlo\Precip_Livneh\');
        infile = strcat(inpath,strcat('prec.',num2str(year),'.nc'));
        Precip = ncread(infile,'prec');
        Precip_pixel = squeeze(Precip(lon_index,lat_index,:));
        %Precip(Precip < 0) = 0;
        Precip_All = [Precip_All; Precip_pixel];
    end      
    Precip_All(1) = [];

	%Extract temperature TS
    [minValue,lat_index] = min(abs(transpose(temp_lat) - sitelat));
	[minValue,lon_index] = min(abs(transpose(temp_lon) - sitelon));
    if sitelon < 0
        [minValue,lon_index] = min(abs(transpose(temp_lon) - (360+sitelon)));
    end
    Tmax_All = zeros(1,1);
    Tmin_All = zeros(1,1);
    for year = 1999:2010
        inpath = strcat('E:\LSM_GW\LSM_JoFlo\Temp_Livneh\');
        infile = strcat(inpath,strcat('tmax.',num2str(year),'.nc'));
        Tmax = ncread(infile,'tmax');
        Tmax_pixel = squeeze(Tmax(lon_index,lat_index,:));
        Tmax_All = [Tmax_All; Tmax_pixel];
        infile = strcat(inpath,strcat('tmin.',num2str(year),'.nc'));
        Tmin = ncread(infile,'tmin');
        Tmin_pixel = squeeze(Tmin(lon_index,lat_index,:));
        Tmin_All = [Tmin_All; Tmin_pixel];
    end      
    Tmax_All(1) = [];
    Tmin_All(1) = [];
    
    ForcingData(:,1) = Precip_All;
    ForcingData(:,2) = Tmax_All;
    ForcingData(:,3) = Tmin_All;
    
    outname = strcat('E:\LSM_GW\LSM_JoFlo\MetForcing_Livneh\',num2str(runset(site,1)),'_MetData.csv');
    dlmwrite(outname,ForcingData)

end