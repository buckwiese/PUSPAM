using JuMP, Ipopt

#Random.seed!(20)  # Set a random seed for reproducibility

#eq for PAM optimizer
# set up model
model = Model(Ipopt.Optimizer)
set_optimizer_attribute(model, "max_iter", 35000)

# Set the environment:
#@NLparameter(model, T   == 20.0 + 273.0)   # K,  temperature
@NLparameter(model, ePlam_env == 1) # M,  external laminarin concentration
@NLparameter(model, ePalg_env == 1) # M,  external alginate concentration
@NLparameter(model, cell_dist == 50) # micrometers, distance between cells

#variables to optimize
@variables(model, begin
    logmu, (start = -15) # 1/s, growth rate

    ePlam_surf_rel >= 0, (start = 0.25)  # M, laminarin polysaccharide at cell surface
    eOlam_rel >= 0, (start = 0.1)       # number, extracellular laminarin oligosaccharide 
    pOlam_rel >= 0, (start = 0.1)       # number, periplasmic laminarin oligosaccharide 
    pdbOlam_rel >= 0, (start = 0.1)     # number, periplasmic debranched laminarin oligosaccharide
    pGlc_rel >= 0, (start = 0.1)        # number, periplasmic glucose

    ePalg_surf_rel >= 0, (start = 0.005)  # M, alginate polysaccharide at cell surface
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
    
    phi_GH16A >=0.0, (start = 0.11)      # relative proteome allocation to extracellular endolaminarase GH16A
    phi_GH16B >=0.0, (start = 0.11)      # relative proteome allocation to extracellular endolaminarase GH16B
    phi_GH16C >=0.0, (start = 0.11)      # relative proteome allocation to extracellular endolaminarase GH16C
    phi_TBDTlam >=0.0, (start = 0.08)    # relative proteome allocation to laminarin TBDT in outer membrane
    phi_GH30A >=0.0, (start = 0.001)      # relative proteome allocation to periplasmic debranching laminarase GH30A
    phi_GH30B >=0.0, (start = 0.001)      # relative proteome allocation to periplasmic debranching laminarase GH30B
    phi_GH17A >=0.0, (start = 0.001)      # relative proteome allocation to periplasmic exolaminarase GH17A
    phi_GH17B >=0.0, (start = 0.001)      # relative proteome allocation to periplasmic exolaminarase GH17B

    phi_PL7A >=0.0, (start = 0.001)       # relative proteome allocation to extracellular sentinel endo-alginate lyase PL7A
    phi_PL7B >=0.0, (start = 0.001)       # relative proteome allocation to extracellular sentinel endo-alginate lyase PL7B
    phi_PL7C >=0.0, (start = 0.001)       # relative proteome allocation to extracellular sentinel endo-alginate lyase PL7C
    phi_PL5A >=0.0, (start = 0.002)       # relative proteome allocation to extracellular M-specific endo-alginate lyase PL5A
    phi_PL5B >=0.0, (start = 0.002)       # relative proteome allocation to extracellular M-specific endo-alginate lyase PL5B
    phi_PL5C >=0.0, (start = 0.002)       # relative proteome allocation to extracellular M-specific endo-alginate lyase PL5C
    phi_PL38A >=0.0, (start = 0.001)      # relative proteome allocation to extracellular G-specific endo-alginate lyase PL38A
    phi_PL38B >=0.0, (start = 0.001)      # relative proteome allocation to extracellular G-specific endo-alginate lyase PL38B
    phi_PL38C >=0.0, (start = 0.001)      # relative proteome allocation to extracellular G-specific endo-alginate lyase PL38C
    phi_TBDTalg >=0.0, (start = 0.025)    # relative proteome allocation to alginate TBDT in outer membrane
    phi_PL17A >=0.0, (start = 0.001)      # relative proteome allocation to periplasmic M-specific exo-alginate lyase PL17A
    phi_PL17B >=0.0, (start = 0.001)      # relative proteome allocation to periplasmic M-specific exo-alginate lyase PL17B
    phi_PL6A >=0.0, (start = 0.001)       # relative proteome allocation to periplasmic G-specific exo-alginate lyase PL6A
    phi_PL6B >=0.0, (start = 0.001)       # relative proteome allocation to periplasmic G-specific exo-alginate lyase PL6B

    phi_UGlc >= 0.0, (start = 0.002)     # relative proteome allocation to glucose uptake
    phi_UUro >= 0.0, (start = 0.001)     # relative proteome allocation to unsaturated uronate uptake
    phi_ACB >= 0.0, (start = 0.001)      # relative proteome allocation to alginate-catabolism bridge enzymes
    phi_C >= 0.0, (start = 0.01)        # relative proteome allocation to catabolism
    phi_P >= 0.0, (start = 0.05)        # relative proteome allocation to protein synthesis
    phi_M >= 0.0, (start = 0.2)        # relative proteome allocation to membrane synthesis
    
    lm_inner >= 0.1, (start = 0.55) # fraction of lipids in inner membrane
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
    GH16A, phi_GH16A * cP * cC_mw/GH16A_mw/Qcp       # number, extracellular endolaminarase GH16A
    GH16B, phi_GH16B * cP * cC_mw/GH16B_mw/Qcp       # number, extracellular endolaminarase GH16B
    GH16C, phi_GH16C * cP * cC_mw/GH16C_mw/Qcp       # number, extracellular endolaminarase GH16C
    TBDTlam, phi_TBDTlam * cP * cC_mw/TBDT_mw/Qcp # number, laminarin TBDT in outer membrane
    GH30A, phi_GH30A * cP * cC_mw/GH30A_mw/Qcp       # number, periplasmic debranching laminarase GH30A
    GH30B, phi_GH30B * cP * cC_mw/GH30B_mw/Qcp       # number, periplasmic debranching laminarase GH30B
    GH17A, phi_GH17A * cP * cC_mw/GH17A_mw/Qcp       # number, periplasmic exolaminarase GH17A
    GH17B, phi_GH17B * cP * cC_mw/GH17B_mw/Qcp       # number, periplasmic exolaminarase GH17B

    PL7A, phi_PL7A * cP * cC_mw/PL7A_mw/Qcp          # number, extracellular sentinel endo-alginate lyase PL7A
    PL7B, phi_PL7B * cP * cC_mw/PL7B_mw/Qcp          # number, extracellular sentinel endo-alginate lyase PL7B
    PL7C, phi_PL7C * cP * cC_mw/PL7C_mw/Qcp          # number, extracellular sentinel endo-alginate lyase PL7C
    PL5A, phi_PL5A * cP * cC_mw/PL5A_mw/Qcp          # number, extracellular M-specific endo-alginate lyase PL5A
    PL5B, phi_PL5B * cP * cC_mw/PL5B_mw/Qcp          # number, extracellular M-specific endo-alginate lyase PL5B
    PL5C, phi_PL5C * cP * cC_mw/PL5C_mw/Qcp          # number, extracellular M-specific endo-alginate lyase PL5C
    PL38A, phi_PL38A * cP * cC_mw/PL38A_mw/Qcp       # number, extracellular G-specific endo-alginate lyase PL38A
    PL38B, phi_PL38B * cP * cC_mw/PL38B_mw/Qcp       # number, extracellular G-specific endo-alginate lyase PL38B
    PL38C, phi_PL38C * cP * cC_mw/PL38C_mw/Qcp       # number, extracellular G-specific endo-alginate lyase PL38C
    TBDTalg, phi_TBDTalg * cP * cC_mw/TBDT_mw/Qcp # number, alginate TBDT in outer membrane
    PL17A, phi_PL17A * cP * cC_mw/PL17A_mw/Qcp       # number, periplasmic M-specific exo-alginate lyase PL17A
    PL17B, phi_PL17B * cP * cC_mw/PL17B_mw/Qcp       # number, periplasmic M-specific exo-alginate lyase PL17B
    PL6A, phi_PL6A * cP * cC_mw/PL6A_mw/Qcp          # number, periplasmic G-specific exo-alginate lyase PL6A
    PL6B, phi_PL6B * cP * cC_mw/PL6B_mw/Qcp          # number, periplasmic G-specific exo-alginate lyase PL6B

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
    GH16A_M, GH16A/(volume_excelrp*1e-15*Av)          # M, endolaminarase GH16A concentration in extracellular reaction space
    GH16B_M, GH16B/(volume_excelrp*1e-15*Av)          # M, endolaminarase GH16B concentration in extracellular reaction space
    GH16C_M, GH16C/(volume_excelrp*1e-15*Av)          # M, endolaminarase GH16C concentration in extracellular reaction space
    TBDTlam_M, TBDTlam/(volume_excelrp*1e-15*Av)      # M, laminarin TBDT concentration in outer membrane
    GH30A_M, GH30A/(volume_pp*1e-15*Av)               # M, debranching laminarase GH30A concentration in periplasm
    GH30B_M, GH30B/(volume_pp*1e-15*Av)               # M, debranching laminarase GH30B concentration in periplasm
    GH17A_M, GH17A/(volume_pp*1e-15*Av)               # M, exolaminarase GH17A concentration in periplasm
    GH17B_M, GH17B/(volume_pp*1e-15*Av)               # M, exolaminarase GH17B concentration in periplasm

    PL7A_M, PL7A/(volume_excelrp*1e-15*Av)            # M, sentinel endo-alginate lyase PL7A concentration in extracellular reaction space
    PL7B_M, PL7B/(volume_excelrp*1e-15*Av)            # M, sentinel endo-alginate lyase PL7B concentration in extracellular reaction space
    PL7C_M, PL7C/(volume_excelrp*1e-15*Av)            # M, sentinel endo-alginate lyase PL7C concentration in extracellular reaction space
    PL5A_M, PL5A/(volume_excelrp*1e-15*Av)            # M, M-specific endo-alginate lyase PL5A concentration in extracellular reaction space
    PL5B_M, PL5B/(volume_excelrp*1e-15*Av)            # M, M-specific endo-alginate lyase PL5B concentration in extracellular reaction space
    PL5C_M, PL5C/(volume_excelrp*1e-15*Av)            # M, M-specific endo-alginate lyase PL5C concentration in extracellular reaction space
    PL38A_M, PL38A/(volume_excelrp*1e-15*Av)          # M, G-specific endo-alginate lyase PL38A concentration in extracellular reaction space
    PL38B_M, PL38B/(volume_excelrp*1e-15*Av)          # M, G-specific endo-alginate lyase PL38B concentration in extracellular reaction space
    PL38C_M, PL38C/(volume_excelrp*1e-15*Av)          # M, G-specific endo-alginate lyase PL38C concentration in extracellular reaction space
    TBDTalg_M, TBDTalg/(volume_excelrp*1e-15*Av)    # M, alginate TBDT concentration in outer membrane
    PL17A_M, PL17A/(volume_pp*1e-15*Av)               # M, M-specific exo-alginate lyase PL17A concentration in periplasm
    PL17B_M, PL17B/(volume_pp*1e-15*Av)               # M, M-specific exo-alginate lyase PL17B concentration in periplasm
    PL6A_M, PL6A/(volume_pp*1e-15*Av)                 # M, G-specific exo-alginate lyase PL6A concentration in periplasm
    PL6B_M, PL6B/(volume_pp*1e-15*Av)                 # M, G-specific exo-alginate lyase PL6B concentration in periplasm
    
    UGlc_M, UGlc/(volume_pp*1e-15*Av)               # M, glucose uptake transporter concentration in periplasm
    UUro_M, UUro/(volume_pp*1e-15*Av)               # M, uronate uptake transporter concentration in periplasm
    ACB_M, ACB/(volume_cy*1e-15*Av)                 # M, alginate-catabolism bridge enzyme concentration in cytoplasm
    CE_M, CE/(volume_cy*1e-15*Av)                   # M, catabolism enzyme concentration in cytoplasm
    PE_M, PE/(volume_cy*1e-15*Av)                   # M, protein synthesis enzyme concentration in cytoplasm
    ME_M, ME/(volume_cy*1e-15*Av)                   # M, membrane synthesis enzyme concentration in cytoplasm
    
    # reaction rates
    v_ENDOlamA, GH16A_kcat * GH16A_M * ePlam_surf / (GH16A_km + ePlam_surf)  # M/s, extracellular laminarin hydrolysis by GH16A
    v_ENDOlamB, GH16B_kcat * GH16B_M * ePlam_surf / (GH16B_km + ePlam_surf)  # M/s, extracellular laminarin hydrolysis by GH16B
    v_ENDOlamC, GH16C_kcat * GH16C_M * ePlam_surf / (GH16C_km + ePlam_surf)  # M/s, extracellular laminarin hydrolysis by GH16C
    v_ENDOlam, v_ENDOlamA + v_ENDOlamB + v_ENDOlamC                          # M/s, extracellular laminarin hydrolysis
    v_IMPRTlam, TBDT_kcat * TBDTlam_M * eOlam_M / (TBDT_km + eOlam_M)        # M/s, laminarin uptake
    v_DEBRlamA, GH30A_kcat * GH30A_M * pOlam_M / (GH30A_km + pOlam_M)        # M/s, periplasmic laminarin debranching by GH30A
    v_DEBRlamB, GH30B_kcat * GH30B_M * pOlam_M / (GH30B_km + pOlam_M)        # M/s, periplasmic laminarin debranching by GH30B
    v_DEBRlam, v_DEBRlamA + v_DEBRlamB                                       # M/s, periplasmic laminarin debranching
    v_EXOlamA, GH17A_kcat * GH17A_M * pdbOlam_M / (GH17A_km + pdbOlam_M)     # M/s, periplasmic laminarin hydrolysis by GH17A
    v_EXOlamB, GH17B_kcat * GH17B_M * pdbOlam_M / (GH17B_km + pdbOlam_M)     # M/s, periplasmic laminarin hydrolysis by GH17B
    v_EXOlam, v_EXOlamA + v_EXOlamB                                          # M/s, periplasmic laminarin hydrolysis

    v_SENTalgA, PL7A_kcat * PL7A_M * ePalg_surf / (PL7A_km + ePalg_surf)     # M/s, extracellular alginate hydrolysis by PL7A
    v_SENTalgB, PL7B_kcat * PL7B_M * ePalg_surf / (PL7B_km + ePalg_surf)     # M/s, extracellular alginate hydrolysis by PL7B
    v_SENTalgC, PL7C_kcat * PL7C_M * ePalg_surf / (PL7C_km + ePalg_surf)     # M/s, extracellular alginate hydrolysis by PL7C
    v_SENTalg, v_SENTalgA + v_SENTalgB + v_SENTalgC                          # M/s, extracellular alginate hydrolysis
    v_ENDOalg_MA, PL5A_kcat * PL5A_M * eLMO_M / (PL5A_km + eLMO_M)           # M/s, extracellular large M-rich oligosaccharide hydrolysis by PL5A
    v_ENDOalg_MB, PL5B_kcat * PL5B_M * eLMO_M / (PL5B_km + eLMO_M)           # M/s, extracellular large M-rich oligosaccharide hydrolysis by PL5B
    v_ENDOalg_MC, PL5C_kcat * PL5C_M * eLMO_M / (PL5C_km + eLMO_M)           # M/s, extracellular large M-rich oligosaccharide hydrolysis by PL5C
    v_ENDOalg_M, v_ENDOalg_MA + v_ENDOalg_MB + v_ENDOalg_MC                  # M/s, extracellular large M-rich oligosaccharide hydrolysis
    v_ENDOalg_GA, PL38A_kcat * PL38A_M * eLGO_M / (PL38A_km + eLGO_M)        # M/s, extracellular large G-rich oligosaccharide hydrolysis by PL38A
    v_ENDOalg_GB, PL38B_kcat * PL38B_M * eLGO_M / (PL38B_km + eLGO_M)        # M/s, extracellular large G-rich oligosaccharide hydrolysis by PL38B
    v_ENDOalg_GC, PL38C_kcat * PL38C_M * eLGO_M / (PL38C_km + eLGO_M)        # M/s, extracellular large G-rich oligosaccharide hydrolysis by PL38C
    v_ENDOalg_G, v_ENDOalg_GA + v_ENDOalg_GB + v_ENDOalg_GC                  # M/s, extracellular large G-rich oligosaccharide hydrolysis
    fraction_IMPRTalg_M, eSMO_M / (eSMO_M + eSGO_M)                          # fraction of extracellular small M-rich oligosaccharide that is M-rich
    v_IMPRTalg_M, fraction_IMPRTalg_M * TBDT_kcat * TBDTalg_M * (eSMO_M+eSGO_M) / (TBDT_km + eSMO_M + eSGO_M)        # M/s, extracellular small M-rich oligosaccharide uptake
    v_IMPRTalg_G, (1-fraction_IMPRTalg_M) * TBDT_kcat * TBDTalg_M * (eSMO_M+eSGO_M) / (TBDT_km + eSMO_M + eSGO_M)    # M/s, extracellular small G-rich oligosaccharide uptake
    v_EXOalg_MA, PL17A_kcat * PL17A_M * pSMO_M / (PL17A_km + pSMO_M)         # M/s, periplasmic small M-rich oligosaccharide hydrolysis by PL17A
    v_EXOalg_MB, PL17B_kcat * PL17B_M * pSMO_M / (PL17B_km + pSMO_M)         # M/s, periplasmic small M-rich oligosaccharide hydrolysis by PL17B
    v_EXOalg_M, v_EXOalg_MA + v_EXOalg_MB                                    # M/s, periplasmic small M-rich oligosaccharide hydrolysis
    v_EXOalg_GA, PL6A_kcat * PL6A_M * pSGO_M / (PL6A_km + pSGO_M)            # M/s, periplasmic small G-rich oligosaccharide hydrolysis by PL6A
    v_EXOalg_GB, PL6B_kcat * PL6B_M * pSGO_M / (PL6B_km + pSGO_M)            # M/s, periplasmic small G-rich oligosaccharide hydrolysis by PL6B
    v_EXOalg_G, v_EXOalg_GA + v_EXOalg_GB                                    # M/s, periplasmic small G-rich oligosaccharide hydrolysis

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
    D_pp, (TBDTlam_M+TBDTalg_M)*TBDT_mw + GH30A_M*GH30A_mw + GH30B_M*GH30B_mw + GH17A_M*GH17A_mw + GH17B_M*GH17B_mw + PL17A_M*PL17A_mw + PL17B_M*PL17B_mw + PL6A_M*PL6A_mw + PL6B_M*PL6B_mw + pOlam_M*Olam_mw + pdbOlam_M*dbOlam_mw + pSMO_M*SOalg_mw + pSGO_M*SOalg_mw + pU_M*Uro_mw + lm_M*(1-lm_inner)*lm_mw # g/L, cell density, assuming c makes up 30% of cell mass
    
    # diffusion processes
    diff_rate_Plam, DC_Plam * (ePlam_env-ePlam_surf)*1e-15*Av * (4*pi*(r_cy+r_pp+delta)^2) / cell_dist     # number/s, diffusion rate of substrate polysaccharide (eP) to cell surface based on Karp-Boss et al. 1996 
    oligoloss_rate_lam, DC_Olam * (eOlam_M)*volume_excelrp*1e-15*Av * 4*pi*(r_cy+r_pp+delta)^2 / 50 # number/s, diffusion rate of oligosaccharide (eO) from cell surface based on Karp-Boss et al. 1996
    diff_rate_Palg, DC_Palg * (ePalg_env-ePalg_surf)*1e-15*Av * (4*pi*(r_cy+r_pp+delta)^2) / cell_dist     # number/s, diffusion rate of substrate polysaccharide (eP) to cell surface based on Karp-Boss et al. 1996 
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
    phi_GH16A + phi_GH16B + phi_GH16C + phi_TBDTlam + phi_GH30A + phi_GH30B + phi_GH17A + phi_GH17B + phi_PL7A + phi_PL7B + phi_PL7C + phi_PL5A + phi_PL5B + phi_PL5C + phi_PL38A + phi_PL38B + phi_PL38C + phi_TBDTalg + phi_PL17A + phi_PL17B + phi_PL6A + phi_PL6B + phi_UGlc + phi_UUro + phi_ACB + phi_C + phi_P + phi_M == 0.5
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
    phi_GH16A * v_P * cC_mw/GH16A_mw/Qcp * volume_cy*1e-15*Av - exp(logmu) * GH16A == 0     # proteome allocation constraint for extracellular endolaminarase GH16A
    phi_GH16B * v_P * cC_mw/GH16B_mw/Qcp * volume_cy*1e-15*Av - exp(logmu) * GH16B == 0     # proteome allocation constraint for extracellular endolaminarase GH16B
    phi_GH16C * v_P * cC_mw/GH16C_mw/Qcp * volume_cy*1e-15*Av - exp(logmu) * GH16C == 0     # proteome allocation constraint for extracellular endolaminarase GH16C
    phi_TBDTlam * v_P * cC_mw/TBDT_mw/Qcp * volume_cy*1e-15*Av - exp(logmu) * TBDTlam == 0     # proteome allocation constraint for TBDT in outer membrane
    phi_GH30A * v_P * cC_mw/GH30A_mw/Qcp * volume_cy*1e-15*Av - exp(logmu) * GH30A == 0     # proteome allocation constraint for periplasmic debranching laminarase GH30A
    phi_GH30B * v_P * cC_mw/GH30B_mw/Qcp * volume_cy*1e-15*Av - exp(logmu) * GH30B == 0     # proteome allocation constraint for periplasmic debranching laminarase GH30B
    phi_GH17A * v_P * cC_mw/GH17A_mw/Qcp * volume_cy*1e-15*Av - exp(logmu) * GH17A == 0     # proteome allocation constraint for periplasmic exolaminarase GH17A
    phi_GH17B * v_P * cC_mw/GH17B_mw/Qcp * volume_cy*1e-15*Av - exp(logmu) * GH17B == 0     # proteome allocation constraint for periplasmic exolaminarase GH17B
    phi_PL7A * v_P * cC_mw/PL7A_mw/Qcp * volume_cy*1e-15*Av - exp(logmu) * PL7A == 0     # proteome allocation constraint for extracellular sentinel endo-alginate lyase PL7A
    phi_PL7B * v_P * cC_mw/PL7B_mw/Qcp * volume_cy*1e-15*Av - exp(logmu) * PL7B == 0     # proteome allocation constraint for extracellular sentinel endo-alginate lyase PL7B
    phi_PL7C * v_P * cC_mw/PL7C_mw/Qcp * volume_cy*1e-15*Av - exp(logmu) * PL7C == 0     # proteome allocation constraint for extracellular sentinel endo-alginate lyase PL7C
    phi_PL5A * v_P * cC_mw/PL5A_mw/Qcp * volume_cy*1e-15*Av - exp(logmu) * PL5A == 0     # proteome allocation constraint for extracellular M-specific endo-alginate lyase PL5A
    phi_PL5B * v_P * cC_mw/PL5B_mw/Qcp * volume_cy*1e-15*Av - exp(logmu) * PL5B == 0     # proteome allocation constraint for extracellular M-specific endo-alginate lyase PL5B
    phi_PL5C * v_P * cC_mw/PL5C_mw/Qcp * volume_cy*1e-15*Av - exp(logmu) * PL5C == 0     # proteome allocation constraint for extracellular M-specific endo-alginate lyase PL5C
    phi_PL38A * v_P * cC_mw/PL38A_mw/Qcp * volume_cy*1e-15*Av - exp(logmu) * PL38A == 0     # proteome allocation constraint for extracellular G-specific endo-alginate lyase PL38A
    phi_PL38B * v_P * cC_mw/PL38B_mw/Qcp * volume_cy*1e-15*Av - exp(logmu) * PL38B == 0     # proteome allocation constraint for extracellular G-specific endo-alginate lyase PL38B
    phi_PL38C * v_P * cC_mw/PL38C_mw/Qcp * volume_cy*1e-15*Av - exp(logmu) * PL38C == 0     # proteome allocation constraint for extracellular G-specific endo-alginate lyase PL38C
    phi_TBDTalg * v_P * cC_mw/TBDT_mw/Qcp * volume_cy*1e-15*Av - exp(logmu) * TBDTalg == 0     # proteome allocation constraint for alginate TBDT in outer membrane
    phi_PL17A * v_P * cC_mw/PL17A_mw/Qcp * volume_cy*1e-15*Av - exp(logmu) * PL17A == 0     # proteome allocation constraint for periplasmic M-specific exo-alginate lyase PL17A
    phi_PL17B * v_P * cC_mw/PL17B_mw/Qcp * volume_cy*1e-15*Av - exp(logmu) * PL17B == 0     # proteome allocation constraint for periplasmic M-specific exo-alginate lyase PL17B
    phi_PL6A * v_P * cC_mw/PL6A_mw/Qcp * volume_cy*1e-15*Av - exp(logmu) * PL6A == 0     # proteome allocation constraint for periplasmic G-specific exo-alginate lyase PL6A
    phi_PL6B * v_P * cC_mw/PL6B_mw/Qcp * volume_cy*1e-15*Av - exp(logmu) * PL6B == 0     # proteome allocation constraint for periplasmic G-specific exo-alginate lyase PL6B
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