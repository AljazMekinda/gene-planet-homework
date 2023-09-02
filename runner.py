# Runner for GENE PLANET HOMEWORK
import os
import pysam
import matplotlib.pyplot as plt
import numpy as np

from functions import GenomeExplorer
from utils.general import setup_logging, load_config

config = load_config(os.path.join('config', "config.yaml"))

project_name = config['project_name']

logger = setup_logging(os.path.join('config', 'logging_config.yaml'),
                       logger_name="runner",
                       filename_prefix=f"{project_name}")

logger.info("Opening bam file")
samfile = pysam.AlignmentFile(os.path.join(config["bam_file"]), "rb")

ge = GenomeExplorer(config)

statistics = ge.get_statistics(samfile)

coverage_dict = ge.get_average_coverage(samfile)

ge.generate_pdf_report(statistics, coverage_dict)
logger.info("Closing bam file")
samfile.close()
