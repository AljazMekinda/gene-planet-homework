# Runner for GENE PLANET HOMEWORK
import os
import pysam

from functions import get_coverages
from utils.general import setup_logging, load_config

config = load_config(os.path.join('config', "config.yaml"))

project_name = config['project_name']

logger = setup_logging(os.path.join('config', 'logging_config.yaml'),
                       logger_name="runner",
                       filename_prefix=f"{project_name}")

logger.info("Opening bam file")
samfile = pysam.AlignmentFile(os.path.join(config["bam_file"]), "rb")

logger.info("Calculating coverage across the genome")
base_order = config["coverage"]["base_order"]
total_coverage_dict, total_coverage, total_length = get_coverages(samfile,
                                                                  base_order)

average_coverage = total_coverage / total_length
# for read in samfile:
#     print(read.query_name, read.reference_name, read.cigarstring)
#
#
# for read in samfile.fetch('1', 10000, 10010):
#     print(read)

logger.info("Closing bam file")
samfile.close()
