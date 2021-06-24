#!/usr/bin/env nextflow

params.nant = 20
antennas = Channel
    .from(0..params.nant-1)

params.pols = ['x', 'y']
polarisations = Channel
    .fromList(params.pols)

process beamform {
    executor 'slurm'
    cpus 4
    time '30m'
    memory '128 GB'

    input:
    tuple val(pol), val(antnum) from polarisations.combine(antennas)

    output:
    tuple val(pol), path("${params.label}_${antnum}_${pol}_f.npy") into spectra

    """
    module load python/2.7.14 
    module load numpy/1.16.3-python-2.7.14 
    module load scipy/1.0.0-python-2.7.14 
    module load astropy/2.0.3-python-2.7.14 
    module load matplotlib/2.2.2-python-2.7.14 
    module load joblib/0.11

    if [ "$pol" = "x" ]; then
        f_vcraft="${params.x_vcraft}"
    else
        f_vcraft="${params.y_vcraft}"
    fi

    python $baseDir/craftcor_tab.py -i ${params.numints} \
                           -n ${params.intlen} \
                           --offset ${params.offset} \
                           --calcfile ${params.calcfile} \
                           --par_set ${params.fcm} \
                           --hwfile ${params.hwfile} \
                           --aips_c ${params.bandpass} \
                           --an $antnum \
                           -o ${params.label}_${antnum}_${pol}_f.npy \
                           --tab \$f_vcraft
    """
}

process sum {
    executor 'slurm'
    cpus 1
    time '15m'
    memory '16 GB'

    input:
    tuple val(pol), path(spectra) from spectra.groupTuple()

    output:
    tuple val(pol), path("${params.label}_sum_${pol}_f.npy") into summed_spectrum

    """
    module load gcc/7.3.0 
    module load openmpi/3.0.0
    module load python/3.7.4
    module load numpy/1.18.2-python-3.7.4

    python3 $baseDir/sum.py --f_dir $workDir \
                   -f ${params.label} \
                   -p $pol \
                   -o ${params.label}_sum_${pol}_f.npy
    """
}

process deripple {
    executor 'slurm'
    cpus 1
    time '1h'
    memory '64 GB'

    input:
    tuple val(pol), path(spectrum) from summed_spectrum

    output:
    tuple val(pol), path("${params.label}_sum_${pol}_f_derippled.npy") into derippled_spectrum

    """
    module load python/2.7.14
    module load numpy/1.16.3-python-2.7.14
    module load scipy/1.0.0-python-2.7.14

    fftlen=\$(( ${params.intlen} * 64 ))

    python $baseDir/deripple.py -f $spectrum \
                       -l $fftlen \
                       -o ${params.label}_sum_${pol}_f_derippled.npy
    """
}

process dedisperse {
    executor 'slurm'
    cpus 1
    time '10m'
    memory '64 GB'

    input:
    tuple val(pol), path(spectrum) from derippled_spectrum

    output:
    tuple val(pol), path("${params.label}_sum_${pol}_f_dedispersed_${params.DM}.npy") into dedispersed_spectrum

    """
    module load gcc/7.3.0
    module load openmpi/3.0.0
    module load python/3.7.4
    module load numpy/1.18.2-python-3.7.4

    python3 $baseDir/dedisperse.py -f $spectrum \
                          --DM ${params.DM} \
                          --f0 ${params.f0} \
                          --bw 336 \
                          -o ${params.label}_sum_${pol}_f_dedispersed_${params.DM}.npy
    """
}

process ifft {
    executor 'slurm'
    cpus 1
    time '10m'
    memory '32 GB'

    input:
    tuple val(pol), path(spectrum) from dedispersed_spectrum

    output:
    path "${params.label}_sum_${pol}_t_${params.DM}.npy" into pol_time_series

    """
    module load gcc/7.3.0
    module load openmpi/3.0.0
    module load python/3.7.4
    module load numpy/1.18.2-python-3.7.4
    module load scipy/1.4.1-python-3.7.4

    python3 $baseDir/ifft.py -f $spectrum -o ${params.label}_sum_${pol}_t_${params.DM}.npy
    """
}

process generate_dynspecs {
    executor 'slurm'
    cpus 1
    time '1h'
    memory '64 GB'

    input:
    path pol_time_series from pol_time_series.collect()

    """
    module load gcc/7.3.0
    module load openmpi/3.0.0
    module load python/3.7.4
    module load numpy/1.18.2-python-3.7.4

    python3 $baseDir/dynspecs.py -x ${params.label}_sum_x_t_${params.DM}.npy \
                        -y ${params.label}_sum_y_t_${params.DM}.npy \
                        -o ${params.label}_sum_!_@_${params.DM}.npy
    
    tar czvf ${params.label}_n${params.intlen}_fulltimeres.tar.gz ${params.label}_sum*${params.DM}.npy

    if [ ! -d $baseDir/output ]; then
        mkdir $baseDir/output
    fi

    if [ ! -d $baseDir/output/${params.label} ]; then
        mkdir $baseDir/output/${params.label}
    fi

    cp ${params.label}_n${params.intlen}_fulltimeres.tar.gz $baseDir/output/${params.label}
    """
}