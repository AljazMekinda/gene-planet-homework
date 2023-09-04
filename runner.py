# Runner for GENE PLANET HOMEWORK
import os
import argparse
import pysam
import matplotlib.pyplot as plt
import numpy as np

from functions import GenomeExplorer
from utils.general import setup_logging, load_config

parser = argparse.ArgumentParser(
    description="Runenr for Gene-Planet Assignment project")

parser.add_argument("--config", help="Choose main config file", type=str,
                    default="config.yaml")
args = parser.parse_args()



config = load_config(os.path.join('config', "runner_config", "config.yaml"))

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
