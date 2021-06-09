bad_samples = set([
    # These two need checking:
    #'cortex_neuron_rp10', "colon_rp4", # hyper outliers
    'HeadlessEmbryoE11', # Changed my mind about including this
    'ventrical_zone_tissue_rp1', # Mbd-seq, moron!
    'Neuron_NCAMp_rp1', 'Neuron_NCAMp_rp2', # ESC-derived
    'striatial neurons', # Fuck you MEFs! Cause me so many fucking problems and contributed to a rejection.
    'cardiomyocytes', # ESC-derived
    'cardiac precursoes', # turned out to be ESC-derived
    # <1e6 mapped reads:
    'megakaryoblast_rp2',
    # Weird distribution:
    'cortex_neuron_rp10', 'cortex_neuron_rp9',
    # This is identical to 'adipocytes white':
    'adipocytes_rp1',
    # outliers, I take the majority view:
    'proB_rp3', 'EpiSC_rp5', 'HSC_rp17', 'HSC_rp18', 'HSC_rp19', 'HSC_rp20',
    # Not clear why these failed:
    'cortical_plate_rp1', 'cortex_neuron_rp3', 'hippocampus_neuron_rp12', 'hippocampus_neuron_rp13', 'hippocampus_neuron_rp14',
    'colon_rp4', 'colon_rp5', 'granule_cells_rp1', 'E105_medial_nasal_prominence_rp1', 'granulocyte_rp7',
    'sperm_rp1', # clusters with epididymis, probably they took too much tissue when dissecting.
    'quadriceps_rp1', 'quadriceps_rp2', 'quadriceps_rp3',
    # These are probably fibroblasts:
    'microglia_rp6', 'microglia_rp7',
    'microglia_rp4', 'microglia_rp5', # These appear heavily contaminated with fibroblasts
    # I changed my mind on these:
    'pgc_rp1', 'pgc_rp2', 'pgc_rp3', 'pgc_rp4', 'pgc_rp5', 'pgc_rp6', 'pgc_rp7', 'pgc_rp8', 'pgc_rp9',
    'leptotene_spermatocytes_rp1', 'leptotene_spermatocytes_rp2',
    'pachytene_spermatocytes_rp1', 'pachytene_spermatocytes_rp2',
    'typeA_spermatogonia_rp1', 'typeA_spermatogonia_rp2',
    'typeB_spermatogonia_rp1', 'typeB_spermatogonia_rp2',
    'prim_typeA_spermatogonia_rp1', 'prim_typeA_spermatogonia_rp2',
    'genital_fatpad_rp1',
    'genital_fatpad_rp2',
    'genital_fatpad_rp3',
    'genital_fatpad_rp4',
    'genital_fatpad_rp5',
    'genital_fatpad_rp6',
    'cranial_neural_crest_mandible_rp1', 'cranial_neural_crest_mandible_rp2', 'cranial_neural_crest_mandible_rp3', # These look like some in vitro ESC-derived thing
    'cerebral_cortex_endothelium_rp1', 'cerebral_cortex_endothelium_rp2', # Also seem ESC-like
    'skin_rp1', # I think they took a lot of adipose when they took this sample, as it is very close to fatpad and mesoderm tissues. It is also an outlier in the PCA
    # They never submitted to GEO:
    'dental_pulp_tissue_rp1',
    #'prostate_luminal_cells_rp1', 'prostate_luminal_cells_rp2', 'prostate_luminal_cells_rp4', 'prostate_luminal_cells_rp5',
    #'ventrical_zone_tissue_rp1', 'microglia_rp4', 'microglia_rp5',
    # Have skewed distributions or weird looking profiles or are otherwise weird:
    "SS_Late_blastocyst_E1_C5", "SS_Late_blastocyst_E1_C6", "SS_Late_blastocyst_E2_C1",
    "SS_Late_blastocyst_E2_C16", "SS_Late_blastocyst_E2_C17", "SS_Late_blastocyst_E2_C5",
    "SS_Late_blastocyst_E2_C7", "SS_Late_blastocyst_E1_C13",
    "SS_Late_blastocyst_E1_C2", "SS_Late_blastocyst_E1_C21", "SS_Late_blastocyst_E1_C20", "SS_Late_blastocyst_E1_C23",
    "SS_Late_blastocyst_E1_C4", "SS_Late_blastocyst_E2_C12",
    "SS_Early_blastocyst_E2_C10", "SS_Early_blastocyst_E2_C12",
    "SS_Early_blastocyst_E2_C15", "SS_Early_blastocyst_E2_C6",
    "SS_Early_blastocyst_E2_C8", "SS_Early_blastocyst_E3_C10",
    "SS_Early_blastocyst_E3_C12", "SS_Early_blastocyst_E3_C14",
    "SS_Early_blastocyst_E3_C17", "SS_Early_blastocyst_E3_C7",
    "SS_Early_blastocyst_E3_C8", "SS_Early_blastocyst_E4_C1",
    "SS_Early_blastocyst_E4_C2", "SS_Early_blastocyst_E4_C3",
    "SS_Early_blastocyst_E4_C5", "SS_Early_blastocyst_E4_C6",
    "SS_Embryo2C_early_E0_C1", "SS_Embryo2C_early_E0_C2",
    "SS_Embryo2C_early_E1_C1", "SS_Embryo2C_early_E2_C2",
    "SS_Embryo2C_early_E3_C1", "SS_Embryo2C_late_E5_C1",
    "SS_Embryo2C_late_E5_C2", "SS_Embryo2C_late_E6_C1", "SS_Embryo2C_late_E8_C1",
    "SS_Embryo2C_mid_E0_C1", "SS_Embryo2C_mid_E0_C2",
    "SS_Embryo8C_E1_C4", "SS_Embryo8C_E1_C5",
    "SS_Embryo8C_E2_C6", "SS_Embryo8C_E2_C7", "SS_Embryo8C_E5_C4", "SS_Embryo8C_E8_C3",
    "SS_Mid_blastocyst_E1_C1", "SS_Mid_blastocyst_E1_C16",
    "SS_Mid_blastocyst_E1_C18", "SS_Mid_blastocyst_E2_C8",
    "SS_Embryo2C_mid_E3_C1", "SS_Embryo8C_E5_C7", "SS_Late_blastocyst_E1_C14",
    "SS_Mid_blastocyst_E1_C2", "SS_Mid_blastocyst_E1_C22", "SS_Mid_blastocyst_E1_C3",
    "SS_Mid_blastocyst_E1_C4", "SS_Mid_blastocyst_E1_C6",
    "SS_Mid_blastocyst_E2_C11", "SS_Mid_blastocyst_E2_C7",
    "SS_Mid_blastocyst_E3_C1", "SS_Mid_blastocyst_E3_C10",
    "SS_Mid_blastocyst_E3_C12", "SS_Mid_blastocyst_E3_C14",
    "SS_Mid_blastocyst_E3_C17", "SS_Mid_blastocyst_E3_C3",
    "SS_Mid_blastocyst_E3_C9", "SS_Morula_E1_C12",
    "SS_Morula_E1_C13", "SS_Morula_E4_C10",
    "SS_Morula_E4_C11", "SS_Morula_E4_C2",
    "SS_Morula_E4_C8", "SS_Morula_E4_C9",
    "SS_Morula_E5_C3", "SS_Morula_E5_C4",
    "SS_Morula_E6_C10", "SS_Morula_E6_C11",
    "SS_Morula_E6_C5", "SS_Morula_E6_C9", "SS_Early_blastocyst_E3_C15",
    "SS_Early_blastocyst_E3_C16", "SS_Embryo2C_mid_E3_C2", "SS_Embryo2C_mid_E4_C1",
    "SS_Embryo2C_mid_E4_C2", "SS_Embryo2C_mid_E5_C1", "SS_Embryo4C_E1_C1", "SS_Embryo8C_E8_C4",
    "SS_Mid_blastocyst_E3_C15"])

if __name__ == '__main__':
    print(len(bad_samples))
