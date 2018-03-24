function normalized_constants= normalize_constants(constants)

normalized_constants.radii=constants.radii/constants.radii(1);
normalized_constants.f=constants.f ;

normalized_constants.mu      = 1;     % Gravitational Const
normalized_constants.Re      = 1;          % Earth radius (km)


normalized_constants.a_tol=1e-5;
normalized_constants.a_maxiter=150;

normalized_constants.trueA2normA=(constants.TU^2/constants.RU);
normalized_constants.normA2trueA=(constants.RU/constants.TU^2);

normalized_constants.trueV2normV=(constants.TU/constants.RU);
normalized_constants.normV2trueV=(constants.RU/constants.TU);

normalized_constants.trueX2normX=(1/constants.RU);
normalized_constants.normX2trueX=(constants.RU);

normalized_constants.trueT2normT=(1/constants.TU);
normalized_constants.normT2trueT=(constants.TU);