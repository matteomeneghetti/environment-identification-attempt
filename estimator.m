function environment(block)

setup(block);

%endfunction

function setup(block)

block.NumInputPorts  = 4;
block.NumOutputPorts = 3;

block.SetPreCompInpPortInfoToDynamic;
block.SetPreCompOutPortInfoToDynamic;

% motor command (u)
block.InputPort(1).Dimensions        = 1;
block.InputPort(1).DatatypeID  = 0;  % double
block.InputPort(1).Complexity  = 'Real';
block.InputPort(1).DirectFeedthrough = true;

% motor position (qm or A)
block.InputPort(2).Dimensions        = 1;
block.InputPort(2).DatatypeID  = 0;  % double
block.InputPort(2).Complexity  = 'Real';
block.InputPort(2).DirectFeedthrough = true;

% Environment torque (tau_e or b)
block.InputPort(3).Dimensions        = 1;
block.InputPort(3).DatatypeID  = 0;  % double
block.InputPort(3).Complexity  = 'Real';
block.InputPort(3).DirectFeedthrough = true;

% Time
block.InputPort(4).Dimensions        = 1;
block.InputPort(4).DatatypeID  = 0;  % double
block.InputPort(4).Complexity  = 'Real';
block.InputPort(4).DirectFeedthrough = true;

% Ke_est
block.OutputPort(1).Dimensions       = 1;
block.OutputPort(1).DatatypeID  = 0; % double
block.OutputPort(1).Complexity  = 'Real';
block.OutputPort(1).SamplingMode = 'Sample';

% DOB output (DOB_y)
block.OutputPort(2).Dimensions       = 1;
block.OutputPort(2).DatatypeID  = 0; % double
block.OutputPort(2).Complexity  = 'Real';
block.OutputPort(2).SamplingMode = 'Sample';

% x0
block.OutputPort(3).Dimensions       = 1;
block.OutputPort(3).DatatypeID  = 0; % double
block.OutputPort(3).Complexity  = 'Real';
block.OutputPort(3).SamplingMode = 'Sample';

% Register parameters
block.NumDialogPrms     = 0;

% Register sample times
%  [0 offset]            : Continuous sample time
%  [positive_num offset] : Discrete sample time
%
%  [-1, 0]               : Inherited sample time
%  [-2, 0]               : Variable sample time
block.SampleTimes = [0 0];

% Specify the block simStateCompliance. The allowed values are:
%    'UnknownSimState', < The default setting; warn and assume DefaultSimState
%    'DefaultSimState', < Same sim state as a built-in block
%    'HasNoSimState',   < No sim state
%    'CustomSimState',  < Has GetSimState and SetSimState methods
%    'DisallowSimState' < Error out when saving or restoring the model sim state
block.SimStateCompliance = 'DefaultSimState';

block.RegBlockMethod('PostPropagationSetup',    @DoPostPropSetup);
block.RegBlockMethod('InitializeConditions', @InitializeConditions);
block.RegBlockMethod('Start', @Start);
block.RegBlockMethod('Outputs', @Outputs);     % Required
block.RegBlockMethod('Update', @Update);
block.RegBlockMethod('Derivatives', @Derivatives);
block.RegBlockMethod('Terminate', @Terminate); % Required
block.RegBlockMethod('SetInputPortSamplingMode', @SetInpPortFrameData);

%end setup

function DoPostPropSetup(block)
block.NumDworks = 5;
  
  block.Dwork(1).Name            = 'running';
  block.Dwork(1).Dimensions      = 1;
  block.Dwork(1).DatatypeID      = 8; % boolean
  block.Dwork(1).Complexity      = 'Real';
  block.Dwork(1).UsedAsDiscState = true;

  block.Dwork(2).Name            = 'Ke_ext';
  block.Dwork(2).Dimensions      = 1;
  block.Dwork(2).DatatypeID      = 0; % double
  block.Dwork(2).Complexity      = 'Real';
  block.Dwork(2).UsedAsDiscState = true;

  block.Dwork(3).Name            = 'Pk';
  block.Dwork(3).Dimensions      = 1;
  block.Dwork(3).DatatypeID      = 0; % double
  block.Dwork(3).Complexity      = 'Real';
  block.Dwork(3).UsedAsDiscState = true;

  block.Dwork(4).Name            = 'DOB_y';
  block.Dwork(4).Dimensions      = 1;
  block.Dwork(4).DatatypeID      = 0; % double
  block.Dwork(4).Complexity      = 'Real';
  block.Dwork(4).UsedAsDiscState = true;

  block.Dwork(5).Name            = 'x0';
  block.Dwork(5).Dimensions      = 1;
  block.Dwork(5).DatatypeID      = 0; % double
  block.Dwork(5).Complexity      = 'Real';
  block.Dwork(5).UsedAsDiscState = true;

function InitializeConditions(block)
%end InitializeConditions


function Start(block)
block.Dwork(1).Data = false;
block.Dwork(2).Data = 0;
block.Dwork(3).Data = 0;
block.Dwork(4).Data = 0;
block.Dwork(5).Data = 0;

%end Start

function Outputs(block)

block.OutputPort(1).Data = block.Dwork(2).Data;
block.OutputPort(2).Data = block.Dwork(4).Data;
block.OutputPort(3).Data = block.Dwork(5).Data;
%end Outputs

function Update(block)

running = block.Dwork(1).Data;
Ak = block.InputPort(2).Data; % theta_m
theta_m = Ak;
M = size(Ak, 1);
b = block.InputPort(3).Data; % tau_e
time = block.InputPort(4).Data;

for i = 1:50
    if ~running
        % Initial step (k = 0)
        Pk = inv(Ak.' * Ak);
        if isfinite(Pk)
            block.Dwork(1).Data = true; % running
            block.Dwork(3).Data = Pk; % Pk
            xk = (Ak.' * Ak) \ (Ak.' * b);
            block.Dwork(2).Data = xk; % xk
        end
    else
    xk = block.Dwork(2).Data;
    Pk = block.Dwork(3).Data/0.9;
    % Use SMW to get Pk
    Pk = Pk - Pk * Ak.' * inv(eye(M) + Ak*Pk*Ak.') * Ak * Pk;
    xk = xk + Pk*Ak.' * (b - Ak*xk);
    block.Dwork(2).Data = xk;
    block.Dwork(3).Data = Pk;
    
    wc = 10; w2 = 100;
    Jm = 0.0104; dm = 0.0;
    tau_e = b;
    Ke_ext = block.Dwork(2).Data;
    u = block.InputPort(1).Data;
    DOB_input = tau_e / Ke_ext;
    
    Q = FilterNthOrder(2, [1 2*wc, w2], [0 0 w2], 0.001);
    invQ = FilterNthOrder(2, [1 2*wc w2], [Jm*w2 dm*w2 0], 0.001);
    DOB_y = invQ.process(DOB_input, time) - Q.process(u, time);
    motor = FilterNthOrder(2, [Jm dm 0], [0 0 -1], 0.001);
    
    x0 = motor.process(DOB_y, time);
    % Ak = theta_m - x0; % theta_m - x0
    block.Dwork(4).Data = DOB_y;
    block.Dwork(5).Data = Pk;
    end
end



%end Update

function Derivatives(block)

%end Derivatives

function Terminate(block)

%end Terminate

function SetInpPortFrameData(block, idx, fd)
  
  block.InputPort(idx).SamplingMode = fd;
  block.OutputPort(1).SamplingMode  = fd;
  
%endfunction