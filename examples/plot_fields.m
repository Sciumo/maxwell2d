% This matlab script plots the dielectric constant field and the Poynting
% vector field for electromagnetic-wave simulations by the maxwell2d
% program. To run this script you need the matlab netcdf toolbox to be
% installed.

% Do we print figures immediately?
is_print = 1;
visible = 'on';

% Get a list of the netcdf files in the current directory
files = dir('*.nc');

% ...or specify the file directly:
%clear files; files(1).name = 'gradient.nc';

% Loop through them
for ifile = 1:length(files)

file = files(ifile).name;
root_name = file(1:end-3);

% Load data
[data,att] = load_nc_struct(file);

% Determine if the scattered component of the field has been calculated
% explicitly 
if isfield(data,'Sx_scat')
  is_scat = 1;
else
  is_scat = 0;
end

is_point_osc = ~isempty(findstr(att.global.config,'point_oscillator'));

% Make the Poynting vector lengths proportional to the amplitude rather
% than the intensity of the wave, so that the dynamic range is more
% discernable 
Smag_plot = (data.Sx.^2 + data.Sy.^2).^(0.25);
Sx_plot = data.Sx./Smag_plot;
Sy_plot = data.Sy./Smag_plot;

if is_scat
  Smag_scat_plot = (data.Sx_scat.^2 + data.Sy_scat.^2).^(0.25);
  Sx_scat_plot = data.Sx_scat./Smag_scat_plot;
  Sy_scat_plot = data.Sy_scat./Smag_scat_plot;
end

% The nominal x and y axes
[nx,ny] = size(Smag_plot);
x = [0.5:nx-0.5];
y = [0.5:ny-0.5];

% Determine the resolution at which the vectors will be plotted
ns = 10;
ix = ns+2:ns:nx-ns;
iy = ns+4:ns:ny-ns;

% The domain of the fields to show; note that the lower end of the domain
% is often messy due to the presence of the oscillator here
nns = 1;
iix = 3:nns:nx-2;
iiy = 8:nns:ny-3;

% Get a list of the complex dielectric constants in the domain
epsilon = data.epsilon_r - i.*data.epsilon_i;
epsilon = round(epsilon.*1e6)./1e6;
eps_list = unique(epsilon);


% Plot the dielectric constant
figure(1)
set(gcf,'Units','inches','position',[0.5 0.5 3.8 2.5],'Resize','on',...
        'papertype','a4','paperorientation','portrait',...
        'paperposition',[0.5 0.5 3.8 2.5],'visible',visible);
clf
set(gcf,'defaultaxesfontsize',12,'defaulttextfontsize',12);
if length(eps_list) > 8
  % If more than 8 different dielectric constant values, treat them as a
  % continuous distribution and use a colour bar of the real values
  ha = axes('position',[0 0 0.65 1]);
  pcolor(x(iix), y(iiy), data.epsilon_r(iiy, iix));
  shading flat
  daspect([1 1 1]);
  set(gca,'layer','top','xtick',[],'ytick',[]);
  set(gca,'visible','off');
  axis([x([1 end]) y([1 end])]);
  c = caxis;
  if c > 1
    c(1) = 1;
  end
  if c(2) < 3;
    c(2) = 3;
  end
  caxis(c);
  h = axes('position',[0.67 0.15 0.04 0.7]);
  axes(ha)
  colorbar(h);
else 
  % If there are 8 or fewer values of dielectric constant, treat them as
  % a discrete list, and report both the real and imaginary parts
  ha = axes('position',[0 0 0.8 1]);
  epstrans = zeros(size(data.epsilon_r));
  for ii = 1:length(eps_list)
    index = find(epsilon == eps_list(ii));
    epstrans(index) = ii;
  end
  n = ii;
  pcolor(x(iix), y(iiy), epstrans(iiy, iix));
  shading flat
  caxis([1 max([length(eps_list) 6])]);
  set(gca,'layer','top','xtick',[],'ytick',[]);
  set(gca,'visible','off');

  daspect([1 1 1]);
  hold on
  pcolor(x(end).*[1.01 1.05], y(end).*linspace(0.25, 0.75, n+1), ...
         [1:n+1;1:n+1]');
  axis([x([1 end]).*[1 1.2] y([1 end])]);
  for ii = 1:n
    text(x(end).*1.05, y(end).*0.25 + ((ii-0.5).* y(end).*0.5)./n, ...
         [' ' num2str(eps_list(ii))]);
  end
end
chiljet
drawnow
if is_print
  print('-dpng','-r80','-painter',[root_name '_epsilon.png']);
end


% Plot the Poynting vector, with the colours indicating (the square root
% of) its magnitude
figure(2)
set(gcf,'Units','inches','position',[0.5 0.5 4 4],'Resize','on',...
        'papertype','a4','paperorientation','portrait',...
        'paperposition',[0.5 0.5 4 4],'visible',visible);
clf
axes('position',[0 0 1 1]);
pcolor(x(iix), y(iiy), Smag_plot(iiy,iix));
shading flat
hold on
quiver(x(ix), y(iy), Sx_plot(iy,ix), Sy_plot(iy,ix),0.9,'k');
chiljet
set(gca,'layer','top','xtick',[],'ytick',[]);
axis([x([1 end]) y([1 end])]);
daspect([1 1 1]);
set(gca,'visible','off');
if is_point_osc
  caxis(caxis./2);
end
thecaxis = caxis;
drawnow
if is_print
  print('-dpng','-r100',[root_name '_poynting.png']);
end

% Plot the Poynting vector of the scattered field
if is_scat
  figure(3)
  set(gcf,'Units','inches','position',[0.5 0.5 4 4],'Resize','on',...
          'papertype','a4','paperorientation','portrait',...
          'paperposition',[0.5 0.5 4 4],'visible',visible);
  clf
  axes('position',[0 0 1 1]);

  pcolor(x(iix), y(iiy), Smag_scat_plot(iiy,iix));
  shading flat
  hold on
  quiver(x(ix), y(iy), Sx_scat_plot(iy,ix), Sy_scat_plot(iy,ix),0.9,'k');
  chiljet
  set(gca,'layer','top','xtick',[],'ytick',[]);
  axis([x([1 end]) y([1 end])]);
  daspect([1 1 1]);
  set(gca,'visible','off');
  caxis(thecaxis);
  drawnow
  if is_print
    print('-dpng','-r100',[root_name '_poynting_scat.png']);
  end
end

end
