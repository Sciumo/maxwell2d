function [data, attribute, dimension] = load_nc_struct(nc_file, names);
% load_nc_struct -- Load NetCDF variables and attributes.
%
% [data, attribute] = load_nc_struct('nc_file') loads all variables of
%   'nc_file' into structure 'data' and all attributes into structure
%   'attribute', so variable 'X' could then be accessed using 'data.X'
%   and attribute 'long_name' of 'X' could be accessed with
%   'attribute.X.long_name'.  Global attributes are accessed via
%   'attribute.global.*'.
%
%   To use this function you need the NetCDF Toolbox for
%   Matlab-5. This can be downloaded free from:
%     http://crusty.er.usgs.gov/~cdenham/MexCDF/nc4ml5.html
 
if nargin < 1, help(mfilename), return, end

result = [];
if nargout > 0, data = []; attribute = []; dimension = []; end

if nargout > 2, require_dimensions = 1; else require_dimensions = 0; end

f = netcdf(nc_file, 'nowrite');
if isempty(f), return, end
disp(['Loading:']);
disp(['   ' name(f)]);

if require_dimensions
  disp(['Dimensions:']);
  dims = dim(f);
  dimension = [];
  for ii = 1:length(dims)
    disp(['     ' name(dims{ii}) ' = ' num2str(length(dims{ii}))]);
    eval(['dimension.' name(dims{ii}) ' = length(dims{ii});']);
  end
end

if nargin < 2
  names = ncnames(var(f));
end

allnames = ncnames(var(f));

disp(['Variables:']);
for ii = 1:length(names)
  if any(strcmp(names{ii},allnames))
    newname = names{ii};
    newname(find(newname == '-')) = '_';
    missing_value = f{names{ii}}.missing_value(:);
    FillValue_ = f{names{ii}}.FillValue_(:);
    if isempty(missing_value) 
      missing_value = FillValue_;
    end
    if ~isempty(missing_value)
      if ischar(missing_value)
	missing_value = str2num(missing_value);
      end
      ncvar = f{names{ii}}(:);
      ncvar(find(ncvar == missing_value)) = NaN;
%      eval(['data.' names{ii} ' = ncvar;']);
      eval(['data.' newname ' = ncvar;']);
      clear ncvar      
    else
%      eval(['data.' names{ii} ' = f{names{ii}}(:);']);
      eval(['data.' newname ' = f{names{ii}}(:);']);
    end
    
    if nargout > 1
      % Add the dimensions
      if require_dimensions
	eval(['attribute.' names{ii} '.dimensions = {};']);
	dims = dim(f{names{ii}});
	for jj = 1:length(dims)
	  eval(['attribute.' names{ii} '.dimensions{jj} = name(dims{jj});']);
	end
      end
      % Do attributes
      atts = att(f{names{ii}});
      for jj = 1:length(atts)
	% Check for underscore starting attribute name
	attname = name(atts{jj});
	if strcmp(attname,'_FillValue')
	  attname = 'FillValue_';
	elseif attname(1) == '_'
	  warning([names{ii} ':' attname ' changed to ' names{ii} ':X' attname]);
	  attname = ['X' attname];
	end
	if ischar(atts{jj}(:))
          % eval(['attribute.' names{ii} '.' attname ' = clean_up_string(atts{jj}(:));']);
	  eval(['attribute.' newname '.' attname ' = clean_up_string(atts{jj}(:));']);
	else
          % eval(['attribute.' names{ii} '.' attname ' = atts{jj}(:);']);
	  eval(['attribute.' newname '.' attname ' = atts{jj}(:);']);
	end
      end
      clear atts
    end
  
    the_size = size(f{names{ii}}(:));
    long_name = clean_up_string(f{names{ii}}.long_name(:));
    if ~isempty(long_name)
      long_name = clean_up_string(long_name);
      long_name = [' (' long_name ')'];
    end
    units = clean_up_string(f{names{ii}}.units(:));
    if strcmp(units,'1')
      % dimensionless
      units = '';
    end
    
    names{ii} = newname;
    
    namefill = blanks(max(0,30-length(newname)));
    
    if the_size(1) == 1 & the_size(2) == 1
      disp([namefill newname ':' 9 num2str(f{names{ii}}(1)) ' ' units 9 9 ...
	    long_name]);
    elseif length(the_size) > 2
      disp([namefill newname ':' 9 '[' num2str(the_size(1)) 'x' ...
	    num2str(the_size(2)) 'x' num2str(the_size(3)) '] ' units ...
	    9 long_name]); 
    else
      disp([namefill newname ':' 9 '[' num2str(the_size(1)) 'x' ...
	    num2str(the_size(2)) '] ' units 9 long_name]); 
    end
  end
end

% Do global attributes
disp('Global attributes:');
attnames = ncnames(att(f));
for ii = 1:length(attnames)
 if isempty(find(attnames{ii} == '/' | attnames{ii} == '-' | attnames{ii} == '+' | attnames{ii} == ' ') > 0)
   eval(['attr = f.' attnames{ii} '(:);']);
   if ischar(attr)
     attr = clean_up_string(attr);
   end
   eval(['attribute.global.' attnames{ii} ' = attr;']);
   namefill = blanks(max(0,14-length(attnames{ii})));
   disp([namefill attnames{ii} ': ' num2str(attr)]);
 else
   disp(['Attribute name ' attnames{ii} ' incompatible with Matlab - not loaded.'])
 end
end
close(f)

function newstr = clean_up_string(oldstr)
newstr = num2str(oldstr);
if length(newstr) > 1
  if newstr(end-1) == '\' & newstr(end) == '0'
    newstr = deblank(newstr(1:end-2));
  end
end
