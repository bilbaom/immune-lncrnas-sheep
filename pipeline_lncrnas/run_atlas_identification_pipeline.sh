#FEELNC paths
export FEELNCPATH=/home/labo/anaconda3/
export PERL5LIB=${FEELNCPATH}/lib/:$PERL5LIB
export PATH=$PATH:${FEELNCPATH}/scripts/
export PATH=$PATH:${FEELNCPATH}/utils/
export PATH=$PATH:${FEELNCPATH}/bin/LINUX/

snakemake -s atlas_identification_pipeline.smk \
          --use-conda \
          --conda-prefix /home/labo/anaconda3/envs \
          --configfile atlas_identification_config.yaml $1
