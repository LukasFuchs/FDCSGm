function [G,R] = SetGSEParam(G,R)
% ======================================================================= %
% Define constants for different grain size evolution models.             %
% ======================================================================= %
disp('    -> Read in gse parameters ...')
G.Kam97     = []; G.Br99    = []; G.AE07    = []; G.Be09    = [];
G.Ro11      = []; G.Dan17   = []; G.MB17    = []; G.Ge20    = [];
if max(strcmpi(G.GSE,'kam97'))
            G.GSE   = 'Kam97';
    % Grain growth ------------------------------------------------------ %
    G.Kam97{1}  = 1e-8;                     % k0 - [ m^2/s ]
    G.Kam97{2}  = 2e5;                      % Q  - [ J/mol ]
    % Grain reduction --------------------------------------------------- %
    G.Kam97{3}  = 2;                        % lambda
    G.Kam97{4}  = 1e5;                      % B  - [ Pa^p m ]
    G.Kam97{5}  = 1.18;                     % p
end
if max(strcmpi(G.GSE,'br99a'))    
            G.GSE   = 'Br99a';    
    G.growth    = 'no';
    % Grain reduction --------------------------------------------------- %
    G.Br99{3}   = 0.4;                      % epsilon_T
    G.Br99{4}   = 15;                       % B  - [ mm MPa^p ]
    G.Br99{5}   = 4/3;                      % p
    G.Br99{4}   = G.Br99{4}*10^(3*(2*G.Br99{5}-1));     % [ m Pa^p ]
end
if max(strcmpi(G.GSE,{'br99b'}))
    switch lower(G.GSE)
        case 'br99b'
            G.GSE       = 'Br99b';
            G.growth    = 'yes';
        case 'all'
            G.GSE       = 'all';
    end
    G.growth    = 'yes';
    % Grain growth ------------------------------------------------------ %
    G.Br99{1}   = 2e-2;                     % [ mm^2/s ]
    G.Br99{2}   = 200e3;                    % [ J/mol ]
    G.Br99{1}   = G.Br99{1}*1e-6;           % [ m^2/s ]
    % Grain reduction --------------------------------------------------- %
    G.Br99{3}   = 0.4;                      % epsilon_T
    G.Br99{4}   = 15;                       % B  - [ mm MPa^p ]
    G.Br99{5}   = 4/3;                      % p
    G.Br99{4}   = G.Br99{4}*10^(3*(2*G.Br99{5}-1));     % [ m Pa^p ]
end
if max(strcmpi(G.GSE,{'ae07','be09','all'}))
    if max(strcmpi(G.GSE,{'ae07','all'}))
        switch lower(G.GSE)
            case 'ae07'
                G.GSE   = 'AE07';
            case 'all'
                G.GSE   = 'all';
        end
        % Grain growth -------------------------------------------------- %
        G.AE07{1}   = 7e4;                  % [ m^p s^-1 ]
        G.AE07{2}   = 520e3;                % [ J mol^-1 ]
        G.AE07{3}   = 2;                    % Grain growth exponent
        % Grain reduction ----------------------------------------------- %
        G.AE07{4}   = 1;                    % gamma - [ J m^-2 ]
        G.AE07{5}   = 0.1;                  % lambda
        G.AE07{6}   = 3;                    % c
    end
    if max(strcmpi(G.GSE,{'be09','all'}))
        switch lower(G.GSE)
            case 'be09'
                G.GSE   = 'Be09';
            case 'all'
                G.GSE   = 'all';
        end
        % Grain growth -------------------------------------------------- %
        G.Be09{1}   = 1.5e-5;               % [ m^p s^-1 ]
        G.Be09{2}   = 350e3;                % [ J mol^-1 ]
        G.Be09{3}   = 3;                    % Grain growth exponent
        % Grain reduction ----------------------------------------------- %
        G.Be09{4}   = 1;                    % gamma - [ J m^-2 ]
        G.Be09{5}   = 0.1;                  % lambda
        G.Be09{6}   = 3;                    % c
    end
end
if max(strcmpi(G.GSE,{'ro11','all'}))
    switch lower(G.GSE)
        case 'ro11'
            G.GSE   = 'Ro11';
        case 'all'
            G.GSE   = 'all';
    end
    
    % Grain growth ------------------------------------------------------ %
    switch R.Param
        case 'wet'
            G.Ro11{1}   = 1.64e4;           % [ �m^p s^-1 ]
            G.Ro11{2}   = 160e3;            % [ J/mol ]
        case 'dry'
            G.Ro11{1}   = 2e4;              % [ �m^p s^-1 ]
            G.Ro11{2}   = 200e3;            % [ J/mol ]
    end
    G.Ro11{3}       = 2;                    % Grain growth exponent
    % Grain reduction --------------------------------------------------- %
    G.Ro11{4}       = 1;                    % gamma - [ J m^-2 ]
    G.Ro11{5}       = 2.054;                % lambda2
    G.Ro11{6}       = 5.053;                % lamdba3
    G.Ro11{7}       = 2;                    % af
    G.Ro11{8}       = 2.9;                  % b
    G.Ro11{1}       = G.Ro11{1}*10^(-6*G.Ro11{3});  % [ m^p s^-1 ]
end
if max(strcmpi(G.GSE,{'br14','all'}))
    switch lower(G.GSE)
        case 'br14'
            G.GSE   = 'BR14';
        case 'all'
            G.GSE   = 'all';
    end
end
if max(strcmpi(G.GSE,{'dan17','all'}))
    switch lower(G.GSE)
        case 'dan17'
            G.GSE   = 'Dan17';
        case 'all'
            G.GSE   = 'all';
    end
    % Grain growth ------------------------------------------------------ %
    G.Dan17{1}      = 1.92e-10;             % [ m^p s^-1 ]
    G.Dan17{2}      = 400e3;                % [ J mol^-1 ]
    G.Dan17{3}      = 3;                    % Grain growth exponent
    % Grain reduction --------------------------------------------------- %
    G.Dan17{4}      = 1;                    % gamma - [ J m^-2 ]
    G.Dan17{5}      = 0.1;                  % lambda
    G.Dan17{6}      = 3;                    % c
end
if max(strcmpi(G.GSE,{'MB17','all'}))
    switch lower(G.GSE)
        case 'mb17'
            G.GSE   = 'MB17';
        case 'all'
            G.GSE   = 'all';
    end
    % Grain growth ------------------------------------------------------ %
    G.MB17{1}       = 2e4;                  % [ microm^p s^-1 ]
    G.MB17{2}       = 300e3;                % [ J/mol ]
    G.MB17{3}       = 2;                    % Grain growth exponent
    G.MB17{1}       = G.MB17{1}*10^(-6*G.MB17{3});  % [ m^p s^-1 ]
    % Interface coarsening ---------------------------------------------- %
    G.MB17{4}       = 4;                    % Exponent
    G.MB17{5}       = G.MB17{4}/G.MB17{3}/250;     % [ microm^(q-p) ]
    G.MB17{5}       = G.MB17{5}*...
        (1e-6)^(G.MB17{4}-G.MB17{3}); % [ m(q-p) ]
    phi1            = 0.4;                  % Phase 1 vol. fract.
    phi2            = 0.6;                  % Phase 2 vol. fract.
    G.MB17{6}       = 3*phi1*phi2;          % Phase distr. funct.
    % Reduction --------------------------------------------------------- %
    G.MB17{7}       = 1;                    % gamma - [ J m^-2 ]
    G.MB17{8}       = 1e-1;                 % Rm*
    G.MB17{9}       = 2.054;                % lambda2
    G.MB17{10}      = 5.053;                % lamdba3
    G.MB17{11}      = 1e-6;                 % F0
    G.MB17{12}      = 4;                    % Mixing transition exp.
end
if max(strcmpi(G.GSE,{'Ge20','all'}))
    switch lower(G.GSE)
        case 'ge20'
            G.GSE   = 'Ge20';
        case 'all'
            G.GSE   = 'all';
    end
    % Grain growth ------------------------------------------------------ %
    G.Ge20{1}       = 2e4;                  % [ microm^p s^-1 ]
    G.Ge20{2}       = 300e3;                % [ J/mol ]
    G.Ge20{3}       = 2;                    % Grain growth exponent
    G.Ge20{1}       = G.Ge20{1}*10^(-6*G.Ge20{3});  % [ m^p s^-1 ]
    % Interface coarsening ---------------------------------------------- %
    G.Ge20{4}       = 4;                    % Exponent
    G.Ge20{5}       = G.Ge20{1}/100*...
        G.Ge20{4}/G.Ge20{3};        % [ GI -  m^p s^-1 ]
    phiOL           = 0.6;                  % Olivine vol. fract.
    phiPX           = 0.4;                  % pyroxene vol. fract.
    G.Ge20{6}       = 3*phiOL*phiPX;        % eta - interface area density
    % Reduction --------------------------------------------------------- %
    G.Ge20{7}       = 1;                    % Interface surface tension
    %       gamma - [ J m^-2 ]
end
% ======================================================================= %
% ============================ END ====================================== %
% ======================================================================= %
end
