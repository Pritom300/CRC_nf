#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.input_csv      = "${projectDir}/data/input/GSE41258_series_matrix.txt.gz"
params.valid_matrix   = "${projectDir}/data/input/GSE14333_series_matrix.txt.gz"
params.valid_clinical = "${projectDir}/data/input/GSE14333_clinical_data.csv"
params.matrix_110223  = "${projectDir}/data/input/GSE110223_series_matrix.txt.gz"
params.matrix_164191  = "${projectDir}/data/input/GSE164191_series_matrix.txt.gz"
params.matrix_142987  = "${projectDir}/data/input/GSE142987_sample_count_matrix.txt.gz"
params.gpl_platform   = "${projectDir}/data/input/GPL570-55999(2).txt"
params.rds_dir        = "${projectDir}/results/rds"
params.tables_dir     = "${projectDir}/results/tables"
params.figures_dir    = "${projectDir}/results/figures"
params.logs_dir       = "${projectDir}/results/logs"

process R_LoadData {
    container 'progression_crc_r:latest'
    publishDir "${params.tables_dir}", mode: 'copy', pattern: "*.csv"
    publishDir "${params.logs_dir}",   mode: 'copy', pattern: "sessionInfo.txt"

    input:
    path gz_file

    output:
    path "GSE41258_Classified_Metadata.csv", emit: metadata
    path "GSE41258_Sample_Table.csv"
    path "Unclassified_Samples.csv", optional: true
    path "sessionInfo.txt"

    script:
    """
    Rscript /project/scripts/r/1.Data_Classification.R \$PWD/${gz_file} \$PWD/
    """
}

process R_LimmaAnalysis {
    container 'progression_crc_r:latest'
    publishDir "${params.rds_dir}",    mode: 'copy', pattern: "*.rds"
    publishDir "${params.tables_dir}", mode: 'copy', pattern: "LIMMA_Results/DEG_*.csv"
    publishDir "${params.figures_dir}", mode: 'copy', pattern: "Figure2_*.png"

    input:
    path gz_file
    path classified_metadata

    output:
    path "expr_train.rds",            emit: expr_train_rds
    path "pheno_train.rds",           emit: pheno_train_rds
    path "expr_test.rds",             emit: expr_test_rds
    path "pheno_test.rds",            emit: pheno_test_rds
    path "expr_clean.rds",            emit: expr_clean_rds
    path "pheno_clean.rds",           emit: pheno_clean_rds
    path "LIMMA_Results/DEG_*.csv",   emit: deg_csvs
    path "LIMMA_Results/",            emit: limma_results_dir  
    path "Figure2_*.png"

    script:
    """
    Rscript "/project/scripts/r/2.Differential Expression Analysis With LIMMA.R" \
        \$PWD/${gz_file} \$PWD/${classified_metadata} \$PWD/
    """
}

process R_Annotation {
    container 'progression_crc_r:latest'
    publishDir "${params.tables_dir}", mode: 'copy', pattern: "*_Q1Mapped.csv"

    input:
    path single_deg_file

    output:
    path "*_Q1Mapped.csv", optional: true

    script:
    """
    Rscript /project/scripts/r/3.Probe_Annotation.R \$PWD/${single_deg_file} \$PWD/
    """
}

process R_StepMiner {
    container 'progression_crc_r:latest'
    publishDir "${params.rds_dir}",    mode: 'copy', pattern: "*.rds"
    publishDir "${params.tables_dir}", mode: 'copy', pattern: "*.csv"

    input:
    path expr_clean_rds
    path pheno_clean_rds
    path deg_csvs 

    output:
    path "StepMiner_Q1_Final_Results.csv",  emit: stepminer_csv
    path "F_Threshold_Sensitivity.csv"
    path "ordered_samples.rds",             emit: ordered_samples_rds
    path "expr_final.rds",                  emit: expr_final_rds
    path "step_results.rds",                emit: step_results_rds
    path "expr_clean.rds",                  emit: expr_clean_rds_out
    path "pheno_clean.rds",                 emit: pheno_clean_rds_out

    script:
    """
    Rscript /project/scripts/r/4.StepMiner_Analysis.R \$PWD/
    """
}

process R_TransitionMapping {
    container 'progression_crc_r:latest'
    publishDir "${params.rds_dir}",    mode: 'copy', pattern: "*.rds"
    publishDir "${params.tables_dir}", mode: 'copy', pattern: "StepMiner_*.csv"
    publishDir "${params.figures_dir}", mode: 'copy', pattern: "*.png"

    input:
    path ordered_samples_rds
    path pheno_clean_rds
    path stepminer_csv

    output:
    path "StepMiner_Q8_Final_Results.csv",  emit: q8_csv
    path "stage_boundaries.rds",            emit: stage_boundaries_rds
    path "ordered_samples.rds",             emit: ordered_samples_rds_out
    path "pheno_clean.rds",                 emit: pheno_clean_rds_out
    path "Figure3_CRC_Progression_Transition.png"

    script:
    """
    Rscript /project/scripts/r/5.Transition_Mapping.R
    """
}
process R_StageClassification {
    container 'progression_crc_r:latest'
    publishDir "${params.rds_dir}",    mode: 'copy', pattern: "*.rds"
    publishDir "${params.tables_dir}", mode: 'copy', pattern: "Stage_Classification_*.csv"
    publishDir "${params.figures_dir}", mode: 'copy', pattern: "*.png"

    input:
    path expr_train_rds
    path expr_test_rds
    path pheno_test_rds
    path expr_clean_rds
    path pheno_clean_rds
    path q8_csv

    output:
    path "Stage_Classification_Final.rds",      emit: trained_rf_model
    path "Stage_Classification_Importance.csv", emit: importance_csv
    path "Figure4A_Confusion_Matrix_Heatmap.png"
    path "Figure4B_Top15_Feature_Importance.png"
    path "selected_genes.rds",                  emit: selected_genes_rds
    path "X_train.rds",                         emit: x_train_rds
    path "y_train.rds",                         emit: y_train_rds
    path "X_test.rds",                          emit: x_test_rds
    path "y_test.rds",                          emit: y_test_rds
    path "common_genes.rds",                    emit: common_genes_rds
    path "predictions.rds"
    path "accuracy.rds",                        emit: accuracy_rds
    path "expr_train.rds"
    path "expr_test.rds",                       emit: expr_test_rds_final
    path "pheno_test.rds",                      emit: pheno_test_rds_final
    path "expr_clean.rds",                      emit: expr_clean_rds_final
    path "pheno_clean.rds",                     emit: pheno_clean_rds_final

    script:
    """
    Rscript /project/scripts/r/6.Stage_Classification.R
    """
}

process R_PermutationStabilityCheck {
    container 'progression_crc_r:latest'
    publishDir "${params.tables_dir}",  mode: 'copy', pattern: "*.csv"
    publishDir "${params.figures_dir}", mode: 'copy', pattern: "*.png"

    input:
    path pheno_clean_rds
    path ordered_samples_rds
    path expr_final_rds
    path stage_boundaries_rds
    path step_results_rds
    path selected_genes_rds

    output:
    path "StepMiner_Permutation_Stability.csv"
    path "FigureS1_StepMiner_Permutation_Stability.png"

    script:
    """
    Rscript /project/scripts/r/7.StepMiner_Permutation_Stability_Check.R
    """
}

process R_ExternalValidation {
    container 'progression_crc_r:latest'
    publishDir "${params.rds_dir}",    mode: 'copy', pattern: "*.rds"
    publishDir "${params.tables_dir}", mode: 'copy', pattern: "*.csv"
    publishDir "${params.figures_dir}", mode: 'copy', pattern: "*.png"

    input:
    path rf_model
    path x_train_rds
    path valid_matrix
    path valid_clinical

    output:
    path "External_Validation_GSE14333_Results.csv"
    path "Figure6_GSE14333_Prediction_Distribution.png"
    path "accuracy_primary.rds",       emit: accuracy_primary_rds
    path "predictions_GSE41333.rds",   emit: predictions_gse_rds

    script:
    """
    Rscript /project/scripts/r/8.External_Validation.R
    """
}

process R_LiverLungOverlapCheck {
    container 'progression_crc_r:latest'
    publishDir "${params.tables_dir}", mode: 'copy', pattern: "*.txt"

    input:
    path q8_csv

    output:
    path "Liver_Lung_Overlap_Check.txt"

    script:
    """
    Rscript "/project/scripts/r/9.LIVER&LUNG_OVERLAP_CHECK.R"
    """
}

process R_MetastasisOrganMarkerCheck {
    container 'progression_crc_r:latest'
    publishDir "${params.tables_dir}", mode: 'copy', pattern: "*.csv"
    publishDir "${params.tables_dir}", mode: 'copy', pattern: "*.txt"

    input:
    path gz_file
    path q8_csv

    output:
    path "Metastasis_Genes_Organ_Baseline_Check.csv"
    path "Metastasis_Genes_Organ_Baseline_Summary.txt"

    script:
    """
    Rscript /project/scripts/r/10.Check_If_METASTASIS_Organ_Markers.R
    """
}

process R_PairedPrimaryMetastasisCheck {
    container 'progression_crc_r:latest'
    publishDir "${params.tables_dir}", mode: 'copy', pattern: "*.csv"
    publishDir "${params.tables_dir}", mode: 'copy', pattern: "*.txt"

    input:
    path gz_file
    path q8_csv

    output:
    path "Metastasis_Genes_Paired_Patient_Check.csv"
    path "Metastasis_Genes_Paired_Agreement_Summary.csv"
    path "Paired_Patient_Check_Summary.txt"

    script:
    """
    Rscript "/project/scripts/r/11.PRIMARY-VS-METASTASIS_Matched_Patients.R"
    """
}

process R_FeatureSizeSensitivity {
    container 'progression_crc_r:latest'
    publishDir "${params.tables_dir}",  mode: 'copy', pattern: "*.csv"
    publishDir "${params.figures_dir}", mode: 'copy', pattern: "*.png"

    input:
    path pheno_clean_rds
    path expr_clean_rds
    path q8_csv

    output:
    path "Feature_Size_Sensitivity.csv"
    path "FigureS2_Feature_Size_Sensitivity.png"

    script:
    """
    Rscript "/project/scripts/r/12.FEATURE-SET_SIZE_SENSITIVITY.R"
    """
}

process R_BootstrapFeatureStability {
    container 'progression_crc_r:latest'
    publishDir "${params.tables_dir}",  mode: 'copy', pattern: "*.csv"
    publishDir "${params.figures_dir}", mode: 'copy', pattern: "*.png"

    input:
    path expr_train_rds
    path pheno_clean_rds
    path common_genes_rds

    output:
    path "Bootstrap_Gene_Stability_All.csv"
    path "Bootstrap_Gene_Stability_FinalPanel.csv"
    path "FigureS3_Bootstrap_Gene_Stability.png"

    script:
    """
    Rscript "/project/scripts/r/13.BOOTSTRAP_FEATURE-SELECTION_STABILITY.R"
    """
}

process R_BootstrapCI {
    container 'progression_crc_r:latest'
    publishDir "${params.tables_dir}", mode: 'copy', pattern: "*.csv"

    input:
    path rf_model
    path expr_test_rds
    path pheno_test_rds
    path common_genes_rds
    path x_test_rds
    path y_test_rds
    path predictions_gse_rds
    path accuracy_rds
    path accuracy_primary_rds

    output:
    path "Test_Set_Bootstrap_CI.csv"
    path "External_Validation_Bootstrap_CI.csv"

    script:
    """
    Rscript /project/scripts/r/14.Bootstrap_95_CI.R
    """
}

process R_BalancedMetricsClassWeight {
    container 'progression_crc_r:latest'
    publishDir "${params.tables_dir}", mode: 'copy', pattern: "*.csv"
    publishDir "${params.tables_dir}", mode: 'copy', pattern: "*.txt"

    input:
    path x_train_rds
    path y_train_rds
    path x_test_rds
    path y_test_rds
    path rf_model

    output:
    path "Class_Weighting_Comparison.csv"
    path "Balanced_Metrics_Class_Weight_Summary.txt"

    script:
    """
    Rscript "/project/scripts/r/15.BALANCED_METRICS_and_CLASS-WEIGHT_SENSITIVITY_CHECK.R"
    """
}

process Python_EnsemblMapping {
    container 'progression_crc_python:latest'
    publishDir "${params.tables_dir}", mode: 'copy', pattern: "*.json"
    publishDir "${params.tables_dir}", mode: 'copy', pattern: "*.txt"

    input:
    path importance_csv

    output:
    path "target_gene_ensg_mapping.json", emit: ensg_json
    path "unresolved_genes.txt"

    script:
    """
    python3 /project/scripts/python/8.Ensembl_Mapping.py \
        \$PWD/${importance_csv} \$PWD/
    """
}

process Python_MultiCohortValidation {
    container 'progression_crc_python:latest'
    publishDir "${params.figures_dir}", mode: 'copy', pattern: "*.png"

    input:
    path importance_csv
    path train_matrix
    path matrix_110223
    path matrix_164191
    path matrix_142987
    path gpl_platform
    path ensg_json

    output:
    path "ROC.png"
    path "SHAP_top10.png"
    path "Figure7_Tissue_Blood_Plasma_ROC.png"

    script:
    """
    python3 /project/scripts/python/9.Multi_Cohort_Validation.py \
        \$PWD/${importance_csv} \$PWD/${train_matrix} \
        \$PWD/${matrix_110223} \$PWD/${matrix_164191} \
        \$PWD/${matrix_142987} \$PWD/${gpl_platform} \
        \$PWD/${ensg_json} \$PWD/
    """
}

process Python_BloodCrossValidation {
    container 'progression_crc_python:latest'
    publishDir "${params.tables_dir}", mode: 'copy', pattern: "*.csv"
    publishDir "${params.figures_dir}", mode: 'copy', pattern: "*.png"

    input:
    path importance_csv
    path matrix_164191
    path gpl_platform

    output:
    path "Validated_Genes_Blood.csv", emit: blood_csv
    path "Final_Blood_ROC.png"
    path "Blood_Panel_Bootstrap_Stability.csv"
    path "Figure8_Blood_Refinement_ROC.png"

    script:
    """
    python3 /project/scripts/python/10.Blood_Cross_Validation.py \
        \$PWD/${importance_csv} \$PWD/${matrix_164191} \
        \$PWD/${gpl_platform} \$PWD/
    """
}

process Python_KEGGAnalysis {
    container 'progression_crc_python:latest'
    publishDir "${params.tables_dir}", mode: 'copy', pattern: "*.csv"
    publishDir "${params.figures_dir}", mode: 'copy', pattern: "*.png"

    input:
    path blood_csv

    output:
    path "Full_KEGG_Analysis_62_Genes.csv"
    path "KEGG_Pathways_Top10.png"

    script:
    """
    python3 /project/scripts/python/11.KEGG_Analysis.py \$PWD/${blood_csv} \$PWD/
    """
}

process Python_GOAnalysis {
    container 'progression_crc_python:latest'
    publishDir "${params.tables_dir}", mode: 'copy', pattern: "*.csv"
    publishDir "${params.figures_dir}", mode: 'copy', pattern: "*.png"

    input:
    path blood_csv

    output:
    path "Full_GO_Analysis_62_Genes.csv"
    path "GO_DotPlot_Top10.png"

    script:
    """
    python3 /project/scripts/python/12.GO_Analysis.py \$PWD/${blood_csv} \$PWD/
    """
}

workflow {
    input_file_ch = Channel.fromPath(params.input_csv)

    R_LoadData(input_file_ch)

    R_LimmaAnalysis(
        input_file_ch,
        R_LoadData.out.metadata
    )

    R_Annotation(R_LimmaAnalysis.out.deg_csvs.flatten())

     R_StepMiner(
        R_LimmaAnalysis.out.expr_clean_rds,
        R_LimmaAnalysis.out.pheno_clean_rds,
        R_LimmaAnalysis.out.limma_results_dir 
    )

        R_TransitionMapping(
        R_StepMiner.out.ordered_samples_rds,
        R_StepMiner.out.pheno_clean_rds_out,
        R_StepMiner.out.stepminer_csv
    )

       R_StageClassification(
        R_LimmaAnalysis.out.expr_train_rds,
        R_LimmaAnalysis.out.expr_test_rds,
        R_LimmaAnalysis.out.pheno_test_rds,
        R_LimmaAnalysis.out.expr_clean_rds,
        R_TransitionMapping.out.pheno_clean_rds_out,
        R_TransitionMapping.out.q8_csv
    )
        R_PermutationStabilityCheck(
        R_StageClassification.out.pheno_clean_rds_final,
        R_TransitionMapping.out.ordered_samples_rds_out,
        R_StepMiner.out.expr_final_rds,
        R_TransitionMapping.out.stage_boundaries_rds,
        R_StepMiner.out.step_results_rds,
        R_StageClassification.out.selected_genes_rds
    )

       R_ExternalValidation(
        R_StageClassification.out.trained_rf_model,
        R_StageClassification.out.x_train_rds,
        Channel.fromPath(params.valid_matrix),
        Channel.fromPath(params.valid_clinical)
    )
    
        R_LiverLungOverlapCheck(R_TransitionMapping.out.q8_csv)
        
        R_MetastasisOrganMarkerCheck(
        input_file_ch,
        R_TransitionMapping.out.q8_csv
    )
    
        R_PairedPrimaryMetastasisCheck(
        input_file_ch,
        R_TransitionMapping.out.q8_csv
    )
    
        R_FeatureSizeSensitivity(
        R_StageClassification.out.pheno_clean_rds_final,
        R_StageClassification.out.expr_clean_rds_final,
        R_TransitionMapping.out.q8_csv
    )
    
    R_BootstrapFeatureStability(
        R_LimmaAnalysis.out.expr_train_rds,
        R_StageClassification.out.pheno_clean_rds_final,
        R_StageClassification.out.common_genes_rds
    )
    
        R_BootstrapCI(
        R_StageClassification.out.trained_rf_model,
        R_StageClassification.out.expr_test_rds_final,
        R_StageClassification.out.pheno_test_rds_final,
        R_StageClassification.out.common_genes_rds,
        R_StageClassification.out.x_test_rds,
        R_StageClassification.out.y_test_rds,
        R_ExternalValidation.out.predictions_gse_rds,
        R_StageClassification.out.accuracy_rds,
        R_ExternalValidation.out.accuracy_primary_rds
    )
    
        R_BalancedMetricsClassWeight(
        R_StageClassification.out.x_train_rds,
        R_StageClassification.out.y_train_rds,
        R_StageClassification.out.x_test_rds,
        R_StageClassification.out.y_test_rds,
        R_StageClassification.out.trained_rf_model
    )
    

    Python_EnsemblMapping(R_StageClassification.out.importance_csv)

    Python_MultiCohortValidation(
        R_StageClassification.out.importance_csv,
        input_file_ch,
        Channel.fromPath(params.matrix_110223),
        Channel.fromPath(params.matrix_164191),
        Channel.fromPath(params.matrix_142987),
        Channel.fromPath(params.gpl_platform),
        Python_EnsemblMapping.out.ensg_json
    )

    Python_BloodCrossValidation(
        R_StageClassification.out.importance_csv,
        Channel.fromPath(params.matrix_164191),
        Channel.fromPath(params.gpl_platform)
    )

    Python_KEGGAnalysis(Python_BloodCrossValidation.out.blood_csv)
    Python_GOAnalysis(Python_BloodCrossValidation.out.blood_csv)
}
