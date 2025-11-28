nextflow.enable.dsl=2

params.samplesheet = '01_metadata/samplesheet.csv' // Provide a default
params.outdir = './results'

// --- Channel for sample information ---
Channel
    .fromPath(params.samplesheet)
    .splitCsv(header:true)
    .map { row -> tuple(row.samples, file(row.path)) }
    .set { samples_ch }

workflow {
    // Step 1: Run QC, filtering, and normalization on each sample in parallel
    QC_AND_FILTER_NORMALISE(samples_ch)

    // Step 2: Collect all the normalized Rds files and integrate them
    // Note: We need to specify which output from the first process we want.
    // .out[0] will be the Rds files, .out[1] will be the metrics.
    INTEGRATE(QC_AND_FILTER_NORMALISE.out.rds_files.collect())
}

process QC_AND_FILTER_NORMALISE {
    // NEW: Publish different outputs to different directories
    publishDir "${params.outdir}/00_qc_reports", mode: 'copy', pattern: "*.qc_metrics.csv"
    publishDir "${params.outdir}/01_normalized_rds", mode: 'copy', pattern: "*.normalised.Rds"


    input:
    tuple val(sample_id), path(raw_data_file)

    // NEW: Define named output channels for clarity
    output:
    path "${sample_id}.normalised.Rds", emit: rds_files
    path "${sample_id}.qc_metrics.csv", emit: qc_metrics
    cpus 4
    memory '16.GB'

    script:
    """
    01_qc_and_filtering.R \\
        --sample_id ${sample_id} \\
        --input_file ${raw_data_file} \\
        --output_rds ${sample_id}.normalised.Rds \\
        --output_qc_metrics ${sample_id}.qc_metrics.csv
    """
}

process INTEGRATE {
    publishDir "${params.outdir}/02_integration", mode: 'copy'

    input:
    path(normalized_rds_list) // This will be a list of files from the previous step

    output:
    path "integrated_seurat_object.Rds", emit: integrated_rds
    //path "integrated_umap.png"
    path "elbow_plot.png"

    // Request more memory and CPUs for this step as it's the most intensive
    cpus 8
    memory '64.GB'

    script:
    """
    02_integrate_samples.R \\
        --input_files ${normalized_rds_list.join(' ')} \\
        --output_rds "integrated_seurat_object.Rds" \\
        --output_plot "integrated_umap.png"
    """
}

process CLUSTER {
    publishDir "${params.outdir}/03_cluster", mode: 'copy'

    input:
    path(integrated_rds)

    output:
    path "clustered_integrated_seurat_dietseurat.rds", emit: rds_clustered_seurat_file
    path "clustered_integrated_sce_dietseurat.rds", emit: rds_clustered_sce_file
    path "*.csv"
    path "cellular-composition-by-sample-integrated_snn_res.0.5*"

    cpus 8
    memory '20.GB'

    script:
    """
    03_cluster.R \\
      --input_rds ${integrated_rds}
    """
}
/*
process ANNOTATION {
    input:
    path(rds_clustered_sce_file)


    output:
    path "annotated_clustered_sce.rds"

    cpus 8
    memory '128.GB'

    script:
    """
    04_annotation.R \\
      --input_rds ${rds_clustered_sce_file} \\
      --output_rds "annotated_clustered_sce.rds"
    """ 

}
*/
