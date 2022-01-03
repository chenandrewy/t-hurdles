function structout = add2struct(struct,varargin)
% adds variables listed in varargin to sstruct
% e.g.: struct = packstruct(structin,x,y,z)

% Andrew Chen 2016 Sept
structout = struct;

Nvar = length(varargin);

for vari = 1:Nvar
    % you need the +1 to account for struct!
    str = ['structout.',inputname(vari+1),' = varargin{vari};'];
    eval(str);
end


