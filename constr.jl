using JuMP, Ipopt

#eq for PAM optimizer
# set up model
model = Model(Ipopt.Optimizer)
set_optimizer_attribute(model, "max_iter", 25000)

# Set the environment:
#@NLparameter(model, T   == 20.0 + 273.0)   # K,  temperature
@NLparameter(model, ePlam_env == 1) # M,  external laminarin concentration
@NLparameter(model, ePalg_env == 1) # M,  external alginate concentration
println("eP_env = ", value(ePalg_env))

#variables to optimize
@variables(model, begin
    logmu, (start = -15) # 1/s, growth rate

    ePlam_surf_rel >= 0, (start = 0.5)  # M, laminarin polysaccharide at cell surface
    eOlam_rel >= 0, (start = 0.1)       # number, extracellular laminarin oligosaccharide 
    pOlam_rel >= 0, (start = 0.1)       # number, periplasmic laminarin oligosaccharide 
    pdbOlam_rel >= 0, (start = 0.1)     # number, periplasmic debranched laminarin oligosaccharide
    pGlc_rel >= 0, (start = 0.1)        # number, periplasmic glucose

    ePalg_surf_rel >= 0, (start = 0.5)  # M, alginate polysaccharide at cell surface
    eLMO_rel >= 0, (start = 0.1)        # number, extracellular large M-rich oligosaccharide 
    eLGO_rel >= 0, (start = 0.1)        # number, extracellular large G-rich oligosaccharide
    eSMO_rel >= 0, (start = 0.1)        # number, extracellular small M-rich oligosaccharide
    eSGO_rel >= 0, (start = 0.1)        # number, extracellular small G-rich oligosaccharide 
    pSMO_rel >= 0, (start = 0.1)        # number, periplasmic M-rich oligosaccharide
    pSGO_rel >= 0, (start = 0.1)        # number, periplasmic G-rich oligosaccharide
    pU_rel >= 0, (start = 0.1)          # number, periplasmic unsaturated uronate
    
    cGlc_rel >= 0, (start = 0.1)    # number, intracellular glucose
    cU_rel >= 0, (start = 0.1)      # number, intracellular unsaturated uronate
    cPyr_rel >= 0, (start = 0.1)    # number, intracellular pyruvate
    cC_rel >= 0, (start = 0.1)      # number, cellular carbon 
    cP_rel >= 0, (start = 0.1)      # number, cellular protein 
    lm_rel >= 0, (start = 0.1)      # number, cellular membrane
    
    phi_GH16 >=0.0, (start = 0.05)      # relative proteome allocation to extracellular endolaminarase GH16
    phi_TBDTlam >=0.0, (start = 0.1)    # relative proteome allocation to laminarin TBDT in outer membrane
    phi_GH30 >=0.0, (start = 0.01)      # relative proteome allocation to periplasmic debranching laminarase GH30
    phi_GH17 >=0.0, (start = 0.01)      # relative proteome allocation to periplasmic exolaminarase GH17

    phi_PL7 >=0.0, (start = 0.05)       # relative proteome allocation to extracellular sentinel endo-alginate lyase PL7
    phi_PL5 >=0.0, (start = 0.05)       # relative proteome allocation to extracellular M-specific endo-alginate lyase PL5
    phi_PL38 >=0.0, (start = 0.05)      # relative proteome allocation to extracellular G-specific endo-alginate lyase PL38
    phi_TBDTalg >=0.0, (start = 0.1)    # relative proteome allocation to alginate TBDT in outer membrane
    phi_PL17 >=0.0, (start = 0.05)      # relative proteome allocation to periplasmic M-specific exo-alginate lyase PL17
    phi_PL6 >=0.0, (start = 0.05)       # relative proteome allocation to periplasmic G-specific exo-alginate lyase PL6

    phi_UGlc >= 0.0, (start = 0.03)     # relative proteome allocation to glucose uptake
    phi_UUro >= 0.0, (start = 0.03)     # relative proteome allocation to unsaturated uronate uptake
    phi_ACB >= 0.0, (start = 0.06)      # relative proteome allocation to alginate-catabolism bridge enzymes
    phi_C >= 0.0, (start = 0.06)        # relative proteome allocation to catabolism
    phi_P >= 0.0, (start = 0.06)        # relative proteome allocation to protein synthesis
    phi_M >= 0.0, (start = 0.06)        # relative proteome allocation to membrane synthesis
    
    lm_inner >= 0.1, (start = 0.5) # fraction of lipids in inner membrane
    C_to_ATP >= 0, (start = 0.5) # mol/mol, carbon to ATP investment ratio
end)
    
# Set the objective: maximize the log growth rate.
@objective(model, Max, logmu)

# Define intermediate calculations.
@NLexpressions(model, begin
    # metabolite numbers
    ePlam_surf, ePlam_surf_rel * ePlam_env  # number, laminarin polysaccharide at cell surface
    eOlam, eOlam_rel * eOlam_max            # number, extracellular laminarin oligosaccharide
    pOlam, pOlam_rel * pOlam_max            # number, periplasmic laminarin oligosaccharide
    pdbOlam, pdbOlam_rel * pdbOlam_max      # number, periplasmic debranched oligosaccharide
    pGlc, pGlc_rel * pGlc_max               # number, periplasmic glucose

    ePalg_surf, ePalg_surf_rel * ePalg_env  # number, alginate polysaccharide at cell surface
    eLMO, eLMO_rel * (eLO_max/2)            # number, extracellular large M-rich oligosaccharide
    eLGO, eLGO_rel * (eLO_max/2) # number, extracellular large G-rich oligosaccharide
    eSMO, eSMO_rel * (eSO_max/2) # number, extracellular small M-rich oligosaccharide
    eSGO, eSGO_rel * (eSO_max/2) # number, extracellular small G-rich oligosaccharide
    pSMO, pSMO_rel * (pSO_max/2) # number, periplasmic small M-rich oligosaccharide
    pSGO, pSGO_rel * (pSO_max/2) # number, periplasmic oligosaccharide
    pU, pU_rel * pU_max          # number, periplasmic unsaturated uronate

    cGlc, cGlc_rel * cGlc_max    # number, intracellular glucose
    cU, cU_rel * cU_max          # number, intracellular unsaturated uronate
    cPyr, cPyr_rel * cPyr_max    # number, intracellular pyruvate
    cC, cC_rel * cC_max          # number, cellular carbon
    cP, cP_rel * cP_max          # number, cellular protein
    lm, lm_rel * lm_max          # number, cellular membrane

    # protein numbers
    GH16, phi_GH16 * cP * cC_mw/GH16_mw/Qcp       # number, extracellular endolaminarase GH16
    TBDTlam, phi_TBDTlam * cP * cC_mw/TBDT_mw/Qcp # number, laminarin TBDT in outer membrane
    GH30, phi_GH30 * cP * cC_mw/GH30_mw/Qcp       # number, periplasmic debranching laminarase GH30
    GH17, phi_GH17 * cP * cC_mw/GH17_mw/Qcp       # number, periplasmic exolaminarase GH17

    PL7, phi_PL7 * cP * cC_mw/PL7_mw/Qcp          # number, extracellular sentinel endo-alginate lyase PL7
    PL5, phi_PL5 * cP * cC_mw/PL5_mw/Qcp          # number, extracellular M-specific endo-alginate lyase PL5
    PL38, phi_PL38 * cP * cC_mw/PL38_mw/Qcp       # number, extracellular G-specific endo-alginate lyase PL38
    TBDTalg, phi_TBDTalg * cP * cC_mw/TBDT_mw/Qcp # number, alginate TBDT in outer membrane
    PL17, phi_PL17 * cP * cC_mw/PL17_mw/Qcp       # number, periplasmic M-specific exo-alginate lyase PL17
    PL6, phi_PL6 * cP * cC_mw/PL6_mw/Qcp          # number, periplasmic G-specific exo-alginate lyase PL6

    UGlc, phi_UGlc * cP * cC_mw/UT_mw/Qcp # number, glucose uptake transporters
    UUro, phi_UUro * cP * cC_mw/UT_mw/Qcp # number, uronate uptake transporters
    ACB, phi_ACB * cP * cC_mw/ACB_mw/Qcp  # number, alginate-catabolism bridge enzymes
    CE, phi_C * cP * cC_mw/C_mw/Qcp       # number, catabolism enzymes
    PE, phi_P * cP * cC_mw/P_mw/Qcp       # number, protein synthesis enzymes
    ME, phi_M * cP * cC_mw/M_mw/Qcp       # number, membrane synthesis enzymes
    
    # cell parameters
    area_cy, sa_lm*lm/2*lm_inner + sa_UT*(UGlc+UUro)                    # μm2, cytoplasm surface area
    r_cy, sqrt(area_cy / (4*pi))                                        # μm, cytoplasm radius
    volume_cy, 4/3*pi*r_cy^3                                            # μm3, cytoplasm volume
    area_pp, sa_lm*lm/2*(1-lm_inner) + sa_TBDT*(TBDTlam+TBDTalg)        # μm2, periplasm surface area
    r_pp, sqrt(area_pp/(4*pi)) -r_cy                                    # μm, periplasm thickness
    volume_pp, 4/3*pi*(r_pp+r_cy)^3 - volume_cy                         # μm3, periplasm volume
    volume_excelrp, 4/3*pi*(r_pp+r_cy+delta)^3 - volume_cy - volume_pp  # μm3, extracellular reaction space volume
    
    # metabolite concentrations
    eOlam_M, eOlam/(volume_excelrp*1e-15*Av)    # M, extracellular laminarin oligosaccharide concentration
    pOlam_M, pOlam/(volume_pp*1e-15*Av)         # M, periplasmic laminarin oligosaccharide concentration
    pdbOlam_M, pdbOlam/(volume_pp*1e-15*Av)     # M, periplasmic debranched laminarin oligosaccharide concentration
    pGlc_M, pGlc/(volume_pp*1e-15*Av)           # M, periplasmic glucose concentration

    eLMO_M, eLMO/(volume_excelrp*1e-15*Av)      # M, extracellular large M-rich oligosaccharide concentration
    eLGO_M, eLGO/(volume_excelrp*1e-15*Av)      # M, extracellular large G-rich oligosaccharide concentration
    eSMO_M, eSMO/(volume_excelrp*1e-15*Av)      # M, extracellular small M-rich oligosaccharide concentration
    eSGO_M, eSGO/(volume_excelrp*1e-15*Av)      # M, extracellular small G-rich oligosaccharide concentration
    pSMO_M, pSMO/(volume_pp*1e-15*Av)           # M, periplasmic small M-rich oligosaccharide concentration
    pSGO_M, pSGO/(volume_pp*1e-15*Av)           # M, periplasmic small G-rich oligosaccharide concentration
    pU_M, pU/(volume_pp*1e-15*Av)               # M, periplasmic unsaturated uronate concentration

    cGlc_M, cGlc/(volume_cy*1e-15*Av)           # M, intracellular glucose concentration
    cU_M, cU/(volume_cy*1e-15*Av)               # M, intracellular unsaturated uronate concentration
    cPyr_M, cPyr/(volume_cy*1e-15*Av)           # M, intracellular pyruvate concentration
    cC_M, cC/(volume_cy*1e-15*Av)               # M, cellular carbon concentration
    cP_M, cP/(volume_cy*1e-15*Av)               # M, cellular protein concentration
    lm_M, lm/(volume_cy*1e-15*Av)               # M, cellular membrane concentration

    # Protein concentrations and distributions
    GH16_M, GH16/(volume_excelrp*1e-15*Av)          # M, endolaminarase GH16 concentration in extracellular reaction space
    TBDTlam_M, TBDTlam/(volume_excelrp*1e-15*Av)    # M, laminarin TBDT concentration in outer membrane
    GH30_M, GH30/(volume_pp*1e-15*Av)               # M, debranching laminarase GH30 concentration in periplasm
    GH17_M, GH17/(volume_pp*1e-15*Av)               # M, exolaminarase GH17 concentration in periplasm

    PL7_M, PL7/(volume_excelrp*1e-15*Av)            # M, sentinel endo-alginate lyase PL7 concentration in extracellular reaction space
    PL5_M, PL5/(volume_excelrp*1e-15*Av)            # M, M-specific endo-alginate lyase PL5 concentration in extracellular reaction space
    PL38_M, PL38/(volume_excelrp*1e-15*Av)          # M, G-specific endo-alginate lyase PL38 concentration in extracellular reaction space
    TBDTalg_M, TBDTalg/(volume_excelrp*1e-15*Av)    # M, alginate TBDT concentration in outer membrane
    PL17_M, PL17/(volume_pp*1e-15*Av)               # M, M-specific exo-alginate lyase PL17 concentration in periplasm
    PL6_M, PL6/(volume_pp*1e-15*Av)                 # M, G-specific exo-alginate lyase PL6 concentration in periplasm
    
    UGlc_M, UGlc/(volume_pp*1e-15*Av)               # M, glucose uptake transporter concentration in periplasm
    UUro_M, UUro/(volume_pp*1e-15*Av)               # M, uronate uptake transporter concentration in periplasm
    ACB_M, ACB/(volume_cy*1e-15*Av)                 # M, alginate-catabolism bridge enzyme concentration in cytoplasm
    CE_M, CE/(volume_cy*1e-15*Av)                   # M, catabolism enzyme concentration in cytoplasm
    PE_M, PE/(volume_cy*1e-15*Av)                   # M, protein synthesis enzyme concentration in cytoplasm
    ME_M, ME/(volume_cy*1e-15*Av)                   # M, membrane synthesis enzyme concentration in cytoplasm
    
    # reaction rates
    v_ENDOlam, GH16_kcat * GH16_M * ePlam_surf / (GH16_km + ePlam_surf)  # M/s, extracellular laminarin hydrolysis
    v_IMPRTlam, TBDT_kcat * TBDTlam_M * eOlam_M / (TBDT_km + eOlam_M)    # M/s, laminarin uptake
    v_DEBRlam, GH30_kcat * GH30_M * pOlam_M / (GH30_km + pOlam_M)        # M/s, periplasmic laminarin debranching
    v_EXOlam, GH17_kcat * GH17_M * pdbOlam_M / (GH17_km + pdbOlam_M)     # M/s, periplasmic laminarin hydrolysis

    v_SENTalg, PL7_kcat * PL7_M * ePalg_surf / (PL7_km + ePalg_surf)        # M/s, extracellular alginate hydrolysis
    v_ENDOalg_M, PL5_kcat * PL5_M * eLMO_M / (PL5_km + eLMO_M)           # M/s, extracellular large M-rich oligosaccharide hydrolysis
    v_ENDOalg_G, PL38_kcat * PL38_M * eLGO_M / (PL38_km + eLGO_M)        # M/s, extracellular large G-rich oligosaccharide hydrolysis
    v_IMPRTalg_M, TBDT_kcat * TBDTalg_M * eSMO_M / (TBDT_km + eSMO_M)    # M/s, extracellular small M-rich oligosaccharide uptake THIS MAY NEED TO BE COMPETITIVE
    v_IMPRTalg_G, TBDT_kcat * TBDTalg_M * eSGO_M / (TBDT_km + eSGO_M)    # M/s, extracellular small G-rich oligosaccharide uptake THIS MAY NEED TO BE COMPETITIVE
    v_EXOalg_M, PL17_kcat * PL17_M * pSMO_M / (PL17_km + pSMO_M)         # M/s, periplasmic small M-rich oligosaccharide hydrolysis
    v_EXOalg_G, PL6_kcat * PL6_M * pSGO_M / (PL6_km + pSGO_M)            # M/s, periplasmic small G-rich oligosaccharide hydrolysis

    v_UGlc, UT_kcat * UGlc_M * pGlc_M / (UT_km + pGlc_M)                 # M/s, glucose uptake
    v_UUro, UT_kcat * UUro_M * pU_M / (UT_km + pU_M)                     # M/s, uronate uptake
    v_ACB, ACB_kcat * ACB_M * cU_M / (ACB_km + cU_M)                     # M/s, alginate-catabolism bridge
    v_CGlc, C_to_ATP * C_kcat * CE_M * cGlc_M / (C_km + cGlc_M)          # M/s, catabolism from glucose
    v_GGlc, (1 - C_to_ATP) * C_kcat * CE_M * cGlc_M / (C_km + cGlc_M)    # M/s, ATP production from glucose
    v_CPyr, C_to_ATP * C_kcat * CE_M * cPyr_M / (C_km + cPyr_M)          # M/s, catabolism from pyruvate
    v_GPyr, (1 - C_to_ATP) * C_kcat * CE_M * cPyr_M / (C_km + cPyr_M)    # M/s, ATP production from pyruvate
    v_P, P_kcat * PE_M * cC_M / (P_km + cC_M)                            # M/s, protein synthesis
    v_M, M_kcat * ME_M * cC_M / (M_km + cC_M)                            # M/s, membrane synthesis

    # cell density
    D_cy, (UGlc_M+UUro_M)*UT_mw + CE_M*C_mw + PE_M*P_mw + ME_M*M_mw + ACB_M*ACB_mw + cC_M*cC_mw*3 + cGlc_M*Glc_mw + cU_M*Uro_mw + cPyr_M*Pyr_mw + lm_M*lm_inner*lm_mw # g/L, cytoplasm density, assuming c makes up 1/3 of cellular carbon metabolite
    D_pp, (TBDTlam_M+TBDTalg_M)*TBDT_mw + GH30_M*GH30_mw + GH17_M*GH17_mw + PL17_M*PL17_mw + PL6_M*PL6_mw + pOlam_M*Olam_mw + pdbOlam_M*dbOlam_mw + pSMO_M*SOalg_mw + pSGO_M*SOalg_mw + pU_M*Uro_mw + lm_M*(1-lm_inner)*lm_mw # g/L, cell density, assuming c makes up 30% of cell mass
    
    # diffusion processes
    diff_rate_Plam, DC_Plam * (ePlam_env-ePlam_surf)*1e-15*Av * (4*pi*(r_cy+r_pp+delta)^2) / 50     # number/s, diffusion rate of substrate polysaccharide (eP) to cell surface based on Karp-Boss et al. 1996 
    oligoloss_rate_lam, DC_Olam * (eOlam_M)*volume_excelrp*1e-15*Av * 4*pi*(r_cy+r_pp+delta)^2 / 50 # number/s, diffusion rate of oligosaccharide (eO) from cell surface based on Karp-Boss et al. 1996
    diff_rate_Palg, DC_Palg * (ePalg_env-ePalg_surf)*1e-15*Av * (4*pi*(r_cy+r_pp+delta)^2) / 50     # number/s, diffusion rate of substrate polysaccharide (eP) to cell surface based on Karp-Boss et al. 1996 
    oligoloss_rate_LM, DC_LOalg * (eLMO_M)*volume_excelrp*1e-15*Av * 4*pi*(r_cy+r_pp+delta)^2 / 50   # number/s, diffusion rate of oligosaccharide (eO) from cell surface based on Karp-Boss et al. 1996
    oligoloss_rate_LG, DC_LOalg * (eLGO_M)*volume_excelrp*1e-15*Av * 4*pi*(r_cy+r_pp+delta)^2 / 50   # number/s, diffusion rate of oligosaccharide (eO) from cell surface based on Karp-Boss et al. 1996
    oligoloss_rate_SM, DC_SOalg * (eSMO_M)*volume_excelrp*1e-15*Av * 4*pi*(r_cy+r_pp+delta)^2 / 50   # number/s, diffusion rate of oligosaccharide (eO) from cell surface based on Karp-Boss et al. 1996
    oligoloss_rate_SG, DC_SOalg * (eSGO_M)*volume_excelrp*1e-15*Av * 4*pi*(r_cy+r_pp+delta)^2 / 50   # number/s, diffusion rate of oligosaccharide (eO) from cell surface based on Karp-Boss et al. 1996
   
    # production rates
    ATP_rate, (stoich_matrix[18, 4]*v_IMPRTlam*volume_excelrp + stoich_matrix[18, 7]*(v_IMPRTalg_M+v_IMPRTalg_G)*volume_excelrp + stoich_matrix[18, 9]*v_UGlc*volume_pp + stoich_matrix[18, 12]*v_UUro*volume_pp + stoich_matrix[18, 10]*v_CGlc*volume_cy + stoich_matrix[18, 11]*v_GGlc*volume_cy + stoich_matrix[18, 14]*v_CPyr*volume_cy + stoich_matrix[18, 15]*v_GPyr*volume_cy + stoich_matrix[18, 16]*v_P*volume_cy + stoich_matrix[18, 17]*v_M*volume_cy)*1e-15*Av # number/s, ATP production rate
    eOlam_rate, (stoich_matrix[2, 1]*v_ENDOlam + stoich_matrix[2, 4]*v_IMPRTlam)*volume_excelrp*1e-15*Av # number/s, extracellular laminarin oligosaccharide production rate
    eLMO_rate, (stoich_matrix[4, 2]*v_SENTalg + stoich_matrix[4, 3]*v_ENDOalg_M)*volume_excelrp*1e-15*Av # number/s, extracellular large M-rich oligosaccharide production rate
    eLGO_rate, (stoich_matrix[5, 2]*v_SENTalg + stoich_matrix[5, 3]*v_ENDOalg_G)*volume_excelrp*1e-15*Av # number/s, extracellular large G-rich oligosaccharide production rate
    eSMO_rate, (stoich_matrix[6, 3]*v_ENDOalg_M + stoich_matrix[6, 7]*v_IMPRTalg_M)*volume_excelrp*1e-15*Av # number/s, extracellular small M-rich oligosaccharide production rate
    eSGO_rate, (stoich_matrix[7, 3]*v_ENDOalg_G + stoich_matrix[7, 7]*v_IMPRTalg_G)*volume_excelrp*1e-15*Av # number/s, extracellular small G-rich oligosaccharide production rate
    pOlam_rate, (stoich_matrix[8, 4]*v_IMPRTlam + stoich_matrix[8, 5]*v_DEBRlam)*volume_pp*1e-15*Av # number/s, periplasmic laminarin oligosaccharide production rate
    pdbOlam_rate, (stoich_matrix[9, 5]*v_DEBRlam + stoich_matrix[9, 6]*v_EXOlam)*volume_pp*1e-15*Av # number/s, periplasmic debranched laminarin oligosaccharide production rate
    pGlc_rate, (stoich_matrix[10, 5]*v_DEBRlam + stoich_matrix[10, 6]*v_EXOlam + stoich_matrix[10, 9]*v_UGlc)*volume_pp*1e-15*Av # number/s, periplasmic glucose production rate
    pSMO_rate, (stoich_matrix[11, 7]*v_IMPRTalg_M + stoich_matrix[11, 8]*v_EXOalg_M)*volume_pp*1e-15*Av # number/s, periplasmic small M-rich oligosaccharide production rate
    pSGO_rate, (stoich_matrix[12, 7]*v_IMPRTalg_G + stoich_matrix[12, 8]*v_EXOalg_G)*volume_pp*1e-15*Av # number/s, periplasmic small G-rich oligosaccharide production rate
    pU_rate, (stoich_matrix[13, 8]*(v_EXOalg_M+v_EXOalg_G) + stoich_matrix[13, 12]*v_UUro)*volume_pp*1e-15*Av # number/s, periplasmic unsaturated uronate production rate
    cGlc_rate, (stoich_matrix[14, 9]*v_UGlc + stoich_matrix[14, 10]*v_CGlc + stoich_matrix[14, 11]*v_GGlc)*volume_cy*1e-15*Av # number/s, intracellular glucose production rate
    cU_rate, (stoich_matrix[15, 12]*v_UUro + stoich_matrix[15, 13]*v_ACB)*volume_cy*1e-15*Av # number/s, intracellular unsaturated uronate production rate
    cPyr_rate, (stoich_matrix[16, 13]*v_ACB + stoich_matrix[16, 14]*v_CPyr + stoich_matrix[16, 15]*v_GPyr)*volume_cy*1e-15*Av # number/s, intracellular pyruvate production rate
    cC_rate, (stoich_matrix[17, 10]*v_CGlc + stoich_matrix[17, 14]*v_CPyr + stoich_matrix[17, 16]*v_P + stoich_matrix[17,17]*v_M)*volume_cy*1e-15*Av # number/s, cellular carbon production rate
    cP_rate, (stoich_matrix[19, 16]*v_P)*volume_cy*1e-15*Av # number/s, cellular protein production rate
    lm_rate, (stoich_matrix[20, 17]*v_M)*volume_cy*1e-15*Av # number/s, cellular membrane production rate
end)
    
# Define the constraints.
# Linear constraints
@constraints(model, begin
    # proteome allocation constraint
    phi_GH16 + phi_TBDTlam + phi_GH30 + phi_GH17 + phi_PL7 + phi_PL5 + phi_PL38 + phi_TBDTalg + phi_PL17 + phi_PL6 + phi_UGlc + phi_UUro + phi_ACB + phi_C + phi_P + phi_M == 0.5
end)

# Nonlinear constraints
@NLconstraints(model, begin
    # constraints on the cell membrane architechture and composition
    (UGlc + UUro) / (lm/2*lm_inner) <= Mmax          # upper boundary of transporter to lipid ratio in IM
    (TBDTlam + TBDTalg)/ (lm/2*(1-lm_inner)) <= Mmax   # upper boundary of transporter to lipid ratio in OM
    lm_inner + (1-lm_inner) == 1.0      # lipid partitioning between IM and OM
    
    # cell stoichiometry constraints
    eOlam_rate - exp(logmu) * eOlam - oligoloss_rate_lam == 0   # number/s, steady-state balanced, growth-rate constraint for extracellular laminarin oligosaccharide
    eLMO_rate - exp(logmu) * eLMO - oligoloss_rate_LM == 0      # number/s, steady-state balanced, growth-rate constraint for extracellular large M-rich oligosaccharide
    eLGO_rate - exp(logmu) * eLGO - oligoloss_rate_LG == 0      # number/s, steady-state balanced, growth-rate constraint for extracellular large G-rich oligosaccharide
    eSMO_rate - exp(logmu) * eSMO - oligoloss_rate_SM == 0      # number/s, steady-state balanced, growth-rate constraint for extracellular small M-rich oligosaccharide
    eSGO_rate - exp(logmu) * eSGO - oligoloss_rate_SG == 0      # number/s, steady-state balanced, growth-rate constraint for extracellular small G-rich oligosaccharide
    pOlam_rate - exp(logmu) * pOlam == 0            # number/s, steady-state balanced, growth-rate constraint for periplasmic laminarin oligosaccharide
    pdbOlam_rate - exp(logmu) * pdbOlam == 0        # number/s, steady-state balanced, growth-rate constraint for periplasmic debranched laminarin oligosaccharide
    pGlc_rate - exp(logmu) * pGlc == 0              # number/s, steady-state balanced, growth-rate constraint for periplasmic glucose
    pSMO_rate - exp(logmu) * pSMO == 0              # number/s, steady-state balanced, growth-rate constraint for periplasmic small M-rich oligosaccharide
    pSGO_rate - exp(logmu) * pSGO == 0              # number/s, steady-state balanced, growth-rate constraint for periplasmic small G-rich oligosaccharide
    pU_rate - exp(logmu) * pU == 0                  # number/s, steady-state balanced, growth-rate constraint for periplasmic unsaturated uronate
    cGlc_rate - exp(logmu) * cGlc == 0              # number/s, steady-state balanced, growth-rate constraint for intracellular glucose
    cU_rate - exp(logmu) * cU == 0                  # number/s, steady-state balanced, growth-rate constraint for intracellular unsaturated uronate
    cPyr_rate - exp(logmu) * cPyr == 0              # number/s, steady-state balanced, growth-rate constraint for intracellular pyruvate
    cC_rate - exp(logmu) * cC == 0                  # number/s, steady-state balanced, growth-rate constraint for cellular carbon
    cP_rate - exp(logmu) * cP == 0                  # number/s, steady-state balanced, growth-rate constraint for cellular protein
    lm_rate - exp(logmu) * lm == 0                  # number/s, steady-state balanced, growth-rate constraint for cellular membrane
    ATP_rate >= 0     # constraint for energy in ATP

    # proteome allocation constraints
    phi_GH16 * v_P * cC_mw/GH16_mw/Qcp * volume_cy*1e-15*Av - exp(logmu) * GH16 == 0     # proteome allocation constraint for extracellular endolaminarase GH16
    phi_TBDTlam * v_P * cC_mw/TBDT_mw/Qcp * volume_cy*1e-15*Av - exp(logmu) * TBDTlam == 0     # proteome allocation constraint for TBDT in outer membrane
    phi_GH30 * v_P * cC_mw/GH30_mw/Qcp * volume_cy*1e-15*Av - exp(logmu) * GH30 == 0     # proteome allocation constraint for periplasmic debranching laminarase GH30
    phi_GH17 * v_P * cC_mw/GH17_mw/Qcp * volume_cy*1e-15*Av - exp(logmu) * GH17 == 0     # proteome allocation constraint for periplasmic exolaminarase GH17
    phi_PL7 * v_P * cC_mw/PL7_mw/Qcp * volume_cy*1e-15*Av - exp(logmu) * PL7 == 0     # proteome allocation constraint for extracellular sentinel endo-alginate lyase PL7
    phi_PL5 * v_P * cC_mw/PL5_mw/Qcp * volume_cy*1e-15*Av - exp(logmu) * PL5 == 0     # proteome allocation constraint for extracellular M-specific endo-alginate lyase PL5
    phi_PL38 * v_P * cC_mw/PL38_mw/Qcp * volume_cy*1e-15*Av - exp(logmu) * PL38 == 0     # proteome allocation constraint for extracellular G-specific endo-alginate lyase PL38
    phi_TBDTalg * v_P * cC_mw/TBDT_mw/Qcp * volume_cy*1e-15*Av - exp(logmu) * TBDTalg == 0     # proteome allocation constraint for alginate TBDT in outer membrane
    phi_PL17 * v_P * cC_mw/PL17_mw/Qcp * volume_cy*1e-15*Av - exp(logmu) * PL17 == 0     # proteome allocation constraint for periplasmic M-specific exo-alginate lyase PL17
    phi_PL6 * v_P * cC_mw/PL6_mw/Qcp * volume_cy*1e-15*Av - exp(logmu) * PL6 == 0     # proteome allocation constraint for periplasmic G-specific exo-alginate lyase PL6
    phi_UGlc * v_P * cC_mw/UT_mw/Qcp * volume_cy*1e-15*Av - exp(logmu) * UGlc == 0     # proteome allocation constraint for glucose uptake transporters
    phi_UUro * v_P * cC_mw/UT_mw/Qcp * volume_cy*1e-15*Av - exp(logmu) * UUro == 0     # proteome allocation constraint for uronate uptake transporters
    phi_ACB * v_P * cC_mw/ACB_mw/Qcp * volume_cy*1e-15*Av - exp(logmu) * ACB == 0     # proteome allocation constraint for alginate-catabolism bridge enzymes
    phi_C * v_P * cC_mw/C_mw/Qcp * volume_cy*1e-15*Av - exp(logmu) * CE == 0     # proteome allocation constraint for catabolism
    phi_P * v_P * cC_mw/P_mw/Qcp * volume_cy*1e-15*Av - exp(logmu) * PE == 0     # proteome allocation constraint for protein synthesis
    phi_M * v_P * cC_mw/M_mw/Qcp * volume_cy*1e-15*Av - exp(logmu) * ME == 0     # proteome allocation constraint for membrane synthesis

    # diffusion constraint
    stoich_matrix[1,1]*v_ENDOlam*volume_excelrp*1e-15*Av + diff_rate_Plam == 0 # number/s, laminarin polysaccharide at cell surface
    stoich_matrix[3,2]*v_SENTalg*volume_excelrp*1e-15*Av + diff_rate_Palg == 0 # number/s, alginate polysaccharide at cell surface

    # upper and lower boundaries of physical parameters
    D_cy <= Dmax                # maximum density constraint; g/L
    D_cy >= Dmin                # minimum density constraint; g/L
    D_pp <= Dmax                # maximum density constraint; g/L
    D_pp >= Dmin                # minimum density constraint; g/L
    r_cy >= 0.02                # minimum cell cytoplasm radius set to 0.05 μm
 #   r_cy <= r_cy_max            # maximum cell cytoplasm radius set to 0.6 μm
    r_pp >= 0.01                # minimum cell periplasm radius set to 0.01 μm
    # boundaries of metabolite pools
    C_to_ATP >= 0           # minimum carbon to ATP investment ratio
    C_to_ATP <= 1           # maximum carbon to ATP investment ratio
    ePlam_surf_rel >= 0     # minimum laminarin polysaccharide concentration at cell surface
    ePlam_surf_rel <= 1.0   # maximum laminarin polysaccharide concentration at cell surface
    ePalg_surf_rel >= 0     # minimum alginate polysaccharide concentration at cell surface
    ePalg_surf_rel <= 1.0   # maximum alginate polysaccharide concentration at cell surface
    eOlam_rel <= 1.0        # maximum extracellular laminarin oligosaccharide concentration
    pOlam_rel <= 1.0        # maximum periplasmic laminarin oligosaccharide concentration
    pdbOlam_rel <= 1.0      # maximum periplasmic debranched laminarin oligosaccharide concentration
    pGlc_rel <= 1.0         # maximum periplasmic glucose concentration
    eLMO_rel <= 1.0         # maximum extracellular large M-rich oligosaccharide concentration
    eLGO_rel <= 1.0         # maximum extracellular large G-rich oligosaccharide concentration
    eSMO_rel <= 1.0         # maximum extracellular small M-rich oligosaccharide concentration
    eSGO_rel <= 1.0         # maximum extracellular small G-rich oligosaccharide concentration
    pSMO_rel <= 1.0         # maximum periplasmic small M-rich oligosaccharide concentration
    pSGO_rel <= 1.0         # maximum periplasmic small G-rich oligosaccharide concentration
    pU_rel <= 1.0           # maximum periplasmic unsaturated uronate concentration
    cGlc_rel <= 1.0         # maximum intracellular glucose concentration
    cU_rel <= 1.0           # maximum intracellular unsaturated uronate concentration
    cPyr_rel <= 1.0         # maximum intracellular pyruvate concentration
    cC_rel <= 1.0           # maximum cellular carbon concentration
    cP_rel <= 1.0           # maximum cellular protein concentration
    lm_rel <= 1.0           # maximum cellular membrane concentration
end)
# Solve the optimization problem.
status = JuMP.optimize!(model)