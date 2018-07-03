function [sys,x0,str,ts,simStateCompliance] = missile(t,x,u,flag)
%SFUNTMPL General MATLAB S-Function Template
%   With MATLAB S-functions, you can define you own ordinary differential
%   equations (ODEs), discrete system equations, and/or just about
%   any type of algorithm to be used within a Simulink block diagram.
%
%   The general form of an MATLAB S-function syntax is:
%       [SYS,X0,STR,TS,SIMSTATECOMPLIANCE] = SFUNC(T,X,U,FLAG,P1,...,Pn)
%
%   What is returned by SFUNC at a given point in time, T, depends on the
%   value of the FLAG, the current state vector, X, and the current
%   input vector, U.
%
%   FLAG   RESULT             DESCRIPTION
%   -----  ------             --------------------------------------------
%   0      [SIZES,X0,STR,TS]  Initialization, return system sizes in SYS,
%                             initial state in X0, state ordering strings
%                             in STR, and sample times in TS.
%   1      DX                 Return continuous state derivatives in SYS.
%   2      DS                 Update discrete states SYS = X(n+1)
%   3      Y                  Return outputs in SYS.
%   4      TNEXT              Return next time hit for variable step sample
%                             time in SYS.
%   5                         Reserved for future (root finding).
%   9      []                 Termination, perform any cleanup SYS=[].
%
%
%   The state vectors, X and X0 consists of continuous states followed
%   by discrete states.
%
%   Optional parameters, P1,...,Pn can be provided to the S-function and
%   used during any FLAG operation.
%
%   When SFUNC is called with FLAG = 0, the following information
%   should be returned:
%
%      SYS(1) = Number of continuous states.
%      SYS(2) = Number of discrete states.
%      SYS(3) = Number of outputs.
%      SYS(4) = Number of inputs.
%               Any of the first four elements in SYS can be specified
%               as -1 indicating that they are dynamically sized. The
%               actual length for all other flags will be equal to the
%               length of the input, U.
%      SYS(5) = Reserved for root finding. Must be zero.
%      SYS(6) = Direct feedthrough flag (1=yes, 0=no). The s-function
%               has direct feedthrough if U is used during the FLAG=3
%               call. Setting this to 0 is akin to making a promise that
%               U will not be used during FLAG=3. If you break the promise
%               then unpredictable results will occur.
%      SYS(7) = Number of sample times. This is the number of rows in TS.
%
%
%      X0     = Initial state conditions or [] if no states.
%
%      STR    = State ordering strings which is generally specified as [].
%
%      TS     = An m-by-2 matrix containing the sample time
%               (period, offset) information. Where m = number of sample
%               times. The ordering of the sample times must be:
%
%               TS = [0      0,      : Continuous sample time.
%                     0      1,      : Continuous, but fixed in minor step
%                                      sample time.
%                     PERIOD OFFSET, : Discrete sample time where
%                                      PERIOD > 0 & OFFSET < PERIOD.
%                     -2     0];     : Variable step discrete sample time
%                                      where FLAG=4 is used to get time of
%                                      next hit.
%
%               There can be more than one sample time providing
%               they are ordered such that they are monotonically
%               increasing. Only the needed sample times should be
%               specified in TS. When specifying more than one
%               sample time, you must check for sample hits explicitly by
%               seeing if
%                  abs(round((T-OFFSET)/PERIOD) - (T-OFFSET)/PERIOD)
%               is within a specified tolerance, generally 1e-8. This
%               tolerance is dependent upon your model's sampling times
%               and simulation time.
%
%               You can also specify that the sample time of the S-function
%               is inherited from the driving block. For functions which
%               change during minor steps, this is done by
%               specifying SYS(7) = 1 and TS = [-1 0]. For functions which
%               are held during minor steps, this is done by specifying
%               SYS(7) = 1 and TS = [-1 1].
%
%      SIMSTATECOMPLIANCE = Specifices how to handle this block when saving and
%                           restoring the complete simulation state of the
%                           model. The allowed values are: 'DefaultSimState',
%                           'HasNoSimState' or 'DisallowSimState'. If this value
%                           is not speficified, then the block's compliance with
%                           simState feature is set to 'UknownSimState'.


%   Copyright 1990-2010 The MathWorks, Inc.

%
% The following outlines the general structure of an S-function.
%
switch flag,

  %%%%%%%%%%%%%%%%%%
  % Initialization %
  %%%%%%%%%%%%%%%%%%
  case 0,
    [sys,x0,str,ts,simStateCompliance]=mdlInitializeSizes;

  %%%%%%%%%%%%%%%
  % Derivatives %
  %%%%%%%%%%%%%%%
  case 1,
    sys=mdlDerivatives(t,x,u);

  %%%%%%%%%%
  % Update %
  %%%%%%%%%%
  case 2,
    sys=mdlUpdate(t,x,u);

  %%%%%%%%%%%
  % Outputs %
  %%%%%%%%%%%
  case 3,
    sys=mdlOutputs(t,x,u);

  %%%%%%%%%%%%%%%%%%%%%%%
  % GetTimeOfNextVarHit %
  %%%%%%%%%%%%%%%%%%%%%%%
  case 4,
    sys=mdlGetTimeOfNextVarHit(t,x,u);

  %%%%%%%%%%%%%
  % Terminate %
  %%%%%%%%%%%%%
  case 9,
    sys=mdlTerminate(t,x,u);

  %%%%%%%%%%%%%%%%%%%%
  % Unexpected flags %
  %%%%%%%%%%%%%%%%%%%%
  otherwise
    DAStudio.error('Simulink:blocks:unhandledFlag', num2str(flag));

end

% end sfuntmpl

%
%=============================================================================
% mdlInitializeSizes
% Return the sizes, initial conditions, and sample times for the S-function.
%=============================================================================
%
function [sys,x0,str,ts,simStateCompliance]=mdlInitializeSizes

%
% call simsizes for a sizes structure, fill it in and convert it to a
% sizes array.
%
% Note that in this example, the values are hard coded.  This is not a
% recommended practice as the characteristics of the block are typically
% defined by the S-function parameters.
%
sizes = simsizes;

sizes.NumContStates  = 12;
sizes.NumDiscStates  = 0;
sizes.NumOutputs     = 6;
sizes.NumInputs      = 3;
sizes.DirFeedthrough = 1;
sizes.NumSampleTimes = 1;   % at least one sample time is needed

sys = simsizes(sizes);

%
% initialize the initial conditions
%
x0  = [100, zeros(1,5),pi/3,zeros(1,5)];

%
% str is always an empty matrix
%
str = [];

%
% initialize the array of sample times
%
ts  = [0 0];

% Specify the block simStateCompliance. The allowed values are:
%    'UnknownSimState', < The default setting; warn and assume DefaultSimState
%    'DefaultSimState', < Same sim state as a built-in block
%    'HasNoSimState',   < No sim state
%    'DisallowSimState' < Error out when saving or restoring the model sim state
simStateCompliance = 'UnknownSimState';

% end mdlInitializeSizes

%
%=============================================================================
% mdlDerivatives
% Return the derivatives for the continuous states.
%=============================================================================
%
function sys=mdlDerivatives(t,x,u)
global M_;
m = 40000;
Cx = 0;
Cya = 2.0;
Cywz = 10.2;
Czwy = -Cywz;
mxb = 0;mxwx = 0;mxwy = 0;
myb = 11.4;mywx = 0;mywy=-58.68;
mza = myb;
mzwz = mywy;
Czb = -2.0;
vx=x(1);vy=x(2);vz=x(3);wx=x(4);wy=x(5);wz=x(6);
theta = x(7);pusy = x(8); phi = x(9); 
xc = -0.9;  yc = 0; zc = 0;
Jx = 289.3; Jy = 6771.8; Jz = 6771.8;
fvx = m*(vz*wy - vy*wz + yc*wx*wy + zc*wx*wz - xc*(wy^2+wz^2));
fvy = m*(vx*wz - vz*wx + xc*wx*wy + zc*wy*wz - yc*(wz^2+wx^2));
fvz = m*(vy*wx - vx*wy + xc*wx*wz + yc*wy*wz - zc*(wy^2+wx^2));
fwx = (Jz-Jy)*wy*wz;
fwy = m*xc*(vx*wy-vy*wx)+(Jx-Jz)*wx*wz;
fwz = m*xc*(vx*wz-vz*wx)+(Jy-Jx)*wx*wy;
dtheta = wy*sin(phi) + wz*cos(phi);
dpusy = 1/cos(theta)*(wy*cos(phi)-wz*sin(phi));
dphi = wx - wy*tan(theta)*cos(phi)+wz*tan(theta)*sin(phi);
L = 12;
S = pi*1^2;  V = S*10.2+S*1.8/3;
q = 1025*V^2/2;
B = 1025*9.8*V;
G = m * 9.8;
a = -atan(vy/vx);   b = atan(vz/sqrt(vx^2+vy^2));
wx_ = wx*L/V; wz_ = wz*L/V;  wy_ = wy*L/V; qS = q*S;
xx = -fvx-Cx*qS -(G-B)*sin(theta);
yy = -fvy + qS*(Cya*a + Cywz*wz_)+(G - B)*cos(theta)*cos(phi);
zz = -fvz + qS*(Czb*b + Czwy*wy_) + (G - B)*cos(theta)*sin(phi);
wxx = -fwx + qS*L*(mxb*b+mxwx*wx_ + mxwy*wy_)+B*cos(theta)*0;%b
wyy = -fwy + qS*L*(myb*b + mywx*wx_ + mywy*wy_)-B*xc*cos(theta)*sin(phi);%b
wzz = -fwz + qS*L*(mza*a + mzwz*wz_)+B*(0-xc*cos(theta)*cos(phi));%a
Cye = 1.467; Czr= -1.467; mxr = 0; mxd = -0.022; myr = -0.615; mze= -0.615;
Mb = [0,0,0;Cye,0,0;0,Czr,0;0,mxr,mxd;0,myr,0;mze,0,0];
ue = u(1);ur = u(2); ud = u(3);    %Ë®Æ½¶æ½Ç ´¹Ö±¶æ½Ç ²î¶¯¶æ½Ç
dv = M_\([xx;yy;zz;wxx;wyy;wzz] + Mb*100*[ue;ur;ud]);
dp = [vx*cos(theta)*cos(pusy) + vy*(sin(pusy)*sin(phi)-sin(theta)*cos(pusy)*cos(phi))+vz*(sin(pusy)*cos(phi)+sin(theta)*cos(pusy)*sin(phi));
    vx*sin(theta)+vy*cos(theta)*cos(phi)-vz*cos(theta)*sin(phi);
    -vx*cos(theta)*sin(pusy)+vy*(cos(pusy)*sin(phi)+sin(theta)*sin(pusy)*cos(phi))+vz*(cos(pusy)*cos(phi)-sin(theta)*sin(pusy)*sin(phi))];
sys = [dv;dtheta;dpusy;dphi;dp]';

% end mdlDerivatives

%
%=============================================================================
% mdlUpdate
% Handle discrete state updates, sample time hits, and major time step
% requirements.
%=============================================================================
%
function sys=mdlUpdate(t,x,u)

sys = [];

% end mdlUpdate

%
%=============================================================================
% mdlOutputs
% Return the block outputs.
%=============================================================================
%
function sys=mdlOutputs(t,x,u)

sys = [x(7),x(8),x(9),x(4),x(5),x(6)];

% end mdlOutputs

%
%=============================================================================
% mdlGetTimeOfNextVarHit
% Return the time of the next hit for this block.  Note that the result is
% absolute time.  Note that this function is only used when you specify a
% variable discrete-time sample time [-2 0] in the sample time array in
% mdlInitializeSizes.
%=============================================================================
%
function sys=mdlGetTimeOfNextVarHit(t,x,u)

sampleTime = 1;    %  Example, set the next hit to be one second later.
sys = t + sampleTime;

% end mdlGetTimeOfNextVarHit

%
%=============================================================================
% mdlTerminate
% Perform any end of simulation tasks.
%=============================================================================
%
function sys=mdlTerminate(t,x,u)

sys = [];

% end mdlTerminate
