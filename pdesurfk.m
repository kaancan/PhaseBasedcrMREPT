function  h=pdesurfk(p,t,u)
%PDESURF Surface plot of PDE solution.
%       PDESURF has been removed. Use PDEPLOT instead.
%
%       PDESURF(P,T,U) plots a 3-D surface of PDE node or triangle data.
%       If U is a column vector, node data is assumed, and continuous
%       style and interpolated shading are used.  If U is a row vector,
%       triangle data is assumed, and discontinuous style and flat shading
%       are used.
%
%       H=PDESURF(P,T,U) additionally returns a handle to the patch objects.
%
%       See also PDEPLOT

%       A. Nordmark 4-26-94, LL 1-30-95.
%       Copyright 1994-2024 The MathWorks, Inc.

% me = MException(message('pde:pdesurf:FunctionToBeRemoved'));
% me = me.addCorrection(matlab.lang.correction.ReplaceIdentifierCorrection('PDESURF', 'PDEPLOT'));
% throw(me);


if size(t,2)==size(u,2)
  cont='discontinuous';
  shad='flat';
elseif size(p,2)==size(u,1)
  cont='continuous';
  shad='interp';
else
  error(message('pde:pdesurf:InvalidSolFormat'));
end

hh=pdeplot(p,[],t,'xydata',u,'xystyle',shad,...
           'zdata',u,'zstyle',cont,'colorbar','off');

if nargout==1
  h=hh;
end

