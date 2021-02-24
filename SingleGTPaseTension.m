function out = SingleGTPaseTension
out{1} = @init;
out{2} = @fun_eval;
out{3} = [];
out{4} = [];
out{5} = [];
out{6} = [];
out{7} = [];
out{8} = [];
out{9} = [];

% --------------------------------------------------------------------------
function dydt = fun_eval(t,kmrgd,bta,b,gam,GT,ell0,phi,Gh,epsilon,alfa,hilln,hillp)
dydt=[(b+bta/(1+exp(-alfa*(kmrgd(2)-(ell0-phi*(kmrgd(1)^hillp)/(Gh^hillp+kmrgd(1)^hillp)))))+gam*(kmrgd(1)^hilln)/(1+kmrgd(1)^hilln))*(GT-kmrgd(1))-kmrgd(1);-epsilon*(kmrgd(2)-(ell0-phi*(kmrgd(1)^hillp)/(Gh^hillp+kmrgd(1)^hillp)));];

% --------------------------------------------------------------------------
function [tspan,y0,options] = init
handles = feval(SingleGTPaseTension);
y0=[0;1];
options = odeset('Jacobian',[],'JacobianP',[],'Hessians',[],'HessiansP',[]);
tspan = [0 10];

% --------------------------------------------------------------------------
function jac = jacobian(t,kmrgd,k0,gam,k,del,hill,T)
% --------------------------------------------------------------------------
function jacp = jacobianp(t,kmrgd,k0,gam,k,del,hill,T)
% --------------------------------------------------------------------------
function hess = hessians(t,kmrgd,k0,gam,k,del,hill,T)
% --------------------------------------------------------------------------
function hessp = hessiansp(t,kmrgd,k0,gam,k,del,hill,T)
%---------------------------------------------------------------------------
function tens3  = der3(t,kmrgd,k0,gam,k,del,hill,T)
%---------------------------------------------------------------------------
function tens4  = der4(t,kmrgd,k0,gam,k,del,hill,T)
%---------------------------------------------------------------------------
function tens5  = der5(t,kmrgd,k0,gam,k,del,hill,T)
