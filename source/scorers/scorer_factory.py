def create_scorer(params):
    if params.scorer_type == 'LogPTestCase':
        from scorers.LogPOctanol.logp_octanol_scorer import LogPOctanolWaterPartitionCoef
        return LogPOctanolWaterPartitionCoef(params)
    elif params.scorer_type == 'SeaLikeTanimoto':
        from scorers.Sea_like_Tanimoto.sea_like_tanimoto_scorer import TanimotoSeaLikeCoef
        return TanimotoSeaLikeCoef(params)
    elif params.scorer_type == 'ampl_model':
        from scorers.AMPL_pred_model.ampl_pred_model import ampl_pred_model
        return ampl_pred_model(params)
    else:
        raise ValueError("Unknown scorer_type: %s" % params.scorer_type)
