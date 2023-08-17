function [R] = SetRheoParam(R,M,N)
% ======================================================================= %
% Function to read in the rheological parameters from an excel file       %
% Rheo.xlsx. The file is stored in the DATA directory.                    %
% All parameters in the data file are (and should be!) scaled to SI units!%
%                                                                         %
% ======================================================================= %
disp('    -> Read in rheological parameters ...')
% ----------------------------------------------------------------------- %
% ----------------------------------------------------------------------- %
R.Param     = 'dry';        % Only used for gse so far!
% File name of rheology data -------------------------------------------- %
fname           = '../../FDCSGm/Rheology/DATA/Rheo.xlsx';
range           = 'A2:M3';
[~,~,raw1]= xlsread(fname,1,range);
[~,~,raw2]= xlsread(fname,2,range);
[~,~,raw3]= xlsread(fname,3,range);
[~,~,raw4]= xlsread(fname,4,range);
% ----------------------------------------------------------------------- %
% Allocate data --------------------------------------------------------- %
switch lower(R.type)
    case 'br99'
        line    = strcmp(raw1(:,2),'Br99');
    case 'hk03'
        line    = strcmp(raw1(:,2),'HK03');
        R.ElM   = cell2mat(raw1(line,12));  % [ J/mol ]
        R.EfM   = cell2mat(raw2(line,12));  % [ J/mol ]
end
% Dislocation creep data ------------------------------------------------ %
R.AlM   = cell2mat(raw1(line,9));           % [ s^-1 Pa^-n ]
R.nM    = cell2mat(raw1(line,10));
R.QlM   = cell2mat(raw1(line,11));          % [ J/mol ]
R.VlM   = cell2mat(raw1(line,13));          % [ J/mol ]
R.AlM   = 0.5*3^((R.nM+1)/2)*R.AlM;
if isnan(R.QlM)
    R.QlM   = R.ElM;
end
% Diffusion Creep data -------------------------------------------------- %
R.AfM   = cell2mat(raw2(line,9));           % [ s^-1 Pa^-1 m^m ]
R.mM    = cell2mat(raw2(line,10));
R.QfM   = cell2mat(raw2(line,11));          % [ J/mol ]
R.VfM   = cell2mat(raw2(line,13));          % [ J/mol ]
R.AfM   = 0.5*3*R.AfM;
if isnan(R.QfM)
    R.QfM   = R.EfM;
end
% ----------------------------------------------------------------------- %

switch lower(M.Type)
    case 'continental'
        line    =   strcmp(raw3(:,2),'Sch09');
        % Rheological parameter - dry Granite - Schmalholz et al. 2009
        R.AUC   =   cell2mat(raw3(line,9));     % [ s^-1 Pa^-1 m^m ]
        R.nUC   =   cell2mat(raw3(line,10));
        R.QUC   =   cell2mat(raw3(line,11));    % [ J/mol ]
        R.EUC   =   cell2mat(raw3(line,12));    % [ J/mol ]
        R.VUC   =   cell2mat(raw3(line,13));    % [ J/mol ]
        R.AUC   =   0.5*3^((R.nUC+1)/2)*R.AUC;  % Correction for deformation experiment
        if isnan(R.QUC)
            R.QUC   =   R.EUC;
        end
        
        if M.zLC ~= 'NaN'
            line    =   strcmp(raw4(:,2),'Sch09');
            % Rheological parameter - Diabase - Schmalholz et al. 2009
            R.ALC   =   cell2mat(raw4(line,9));     % [ s^-1 Pa^-1 m^m ]
            R.nLC   =   cell2mat(raw4(line,10));
            R.QLC   =   cell2mat(raw4(line,11));    % [ J/mol ]
            R.ELC   =   cell2mat(raw4(line,12));    % [ J/mol ]
            R.VLC   =   cell2mat(raw4(line,13));    % [ J/mol ]
            R.ALC   =   0.5*3^((R.nLC+1)/2)*R.ALC;  % Correction for deformation experiment
            if isnan(R.QLC)
                R.QLC   =   R.ELC;
            end
        end
end
clear fname range raw1 raw2

switch lower(M.Type)
    case 'oceanic'
        % Dislocation creep --------------------------------------------- %
        R.Al    =   R.AlM;
        R.n     =   R.nM;
        R.Vl    =   R.VlM;
        R.El    =   R.ElM;
        R.Ql    =   R.QlM;
        % Diffusion Creep ----------------------------------------------- %
        R.Af    =   R.AfM; 
        R.m     =   R.mM;
        R.Ef    =   R.EfM; 
        R.Vf    =   R.VfM;
        R.Qf    =   R.QfM; 
    case 'continental'
        R.Al    =   zeros(N.nz,1);  R.n     =   zeros(N.nz,1);
        R.Vl    =   zeros(N.nz,1);  R.El    =   zeros(N.nz,1);
        R.Ql    =   zeros(N.nz,1);  
        R.Af    =   zeros(N.nz,1);  R.m     =   zeros(N.nz,1);
        R.Ef    =   zeros(N.nz,1);  R.Vf    =   zeros(N.nz,1);
        R.Qf    =   zeros(N.nz,1); 
        % Dislocation creep --------------------------------------------- %
        % Upper crust ---
        R.Al(M.UCind)   =   R.AUC;
        R.n(M.UCind)    =   R.nUC;
        R.Vl(M.UCind)   =   R.VUC;
        R.El(M.UCind)   =   R.EUC;
        R.Ql(M.UCind)   =   R.QUC;
        % Lower Crust ---
        R.Al(M.LCind)   =   R.ALC;
        R.n(M.LCind)    =   R.nLC;
        R.Vl(M.LCind)   =   R.VLC;
        R.El(M.LCind)   =   R.ELC;
        R.Ql(M.LCind)   =   R.QLC;
        % Mantle ---
        R.Al(M.Mind)    =   R.AlM;
        R.n(M.Mind)     =   R.nM;
        R.Vl(M.Mind)    =   R.VlM;
        R.El(M.Mind)    =   R.ElM;
        R.Ql(M.Mind)    =   R.QlM;
%         % Diffusion Creep ----------------------------------------------- %
%         % Upper Crust ---
%         R.Af(M.Mind)    =   R.AfM; 
%         R.m(M.Mind)     =   R.mM;
%         R.Ef(M.Mind)    =   R.EfM; 
%         R.Vf(M.Mind)    =   R.VfM;
%         R.Qf(M.Mind)    =   R.QfM; 
%         % Lower Crust ---
%         R.Af(M.Mind)    =   R.AfM; 
%         R.m(M.Mind)     =   R.mM;
%         R.Ef(M.Mind)    =   R.EfM; 
%         R.Vf(M.Mind)    =   R.VfM;
%         R.Qf(M.Mind)    =   R.QfM; 
        % Mantle ---
        R.Af(:) =   R.AfM; 
        R.m(:)  =   R.mM;
        R.Ef(:) =   R.EfM; 
        R.Vf(:) =   R.VfM;
        R.Qf(:) =   R.QfM; 
end

% ======================================================================= %
% ============================ END ====================================== %
% ======================================================================= %
end
