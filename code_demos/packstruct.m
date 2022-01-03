function structout = packstruct(varargin)
% packs variables listed in varargin into structout
% e.g.: struct = packstruct(x,y,z)

% Copyright (c) Andrew Y. Chen 2012-2013


Nvar = nargin;

for vari = 1:Nvar
    str = ['structout.',inputname(vari),' = varargin{vari};'];
    eval(str);
end


