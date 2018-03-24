function constants= get_constants(config)
% load all the constants required for simulation

% if config.useSPICE=True, then use the up-to-date earth properties
if config.useSPICE
    constants.radii=cspice_bodvrd( 'EARTH', 'RADII', 3);
else
    constants.radii=[6378.137,6378.137,6378.137];
end
constants.f = (constants.radii(1)-constants.radii(3))/constants.radii(1);

constants.mu      = 3.986004418e5;     % Gravitational Const
constants.Re      = constants.radii(1);          % Earth radius (km)

% Canonical Units
constants.muCan   = 1;
constants.RU      = constants.Re;
constants.TU      = sqrt(constants.RU^3 / constants.mu);
constants.VU      = constants.RU/constants.TU;

constants.a_tol=1e-5;
constants.a_maxiter=150;

constants.trueA2normA=(constants.TU^2/constants.RU);
constants.normA2trueA=(constants.RU/constants.TU^2);

constants.trueV2normV=(constants.TU/constants.RU);
constants.normV2trueV=(constants.RU/constants.TU);

constants.trueX2normX=(1/constants.RU);
constants.normX2trueX=(constants.RU);

constants.trueT2normT=(1/constants.TU);
constants.normT2trueT=(constants.TU);
