import logging
import matplotlib.pyplot as plt
import numpy as np
import pysam
# from pyspark.sql import SparkSession
from reportlab.lib import colors
from reportlab.lib.pagesizes import landscape, letter
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Image, \
    Table, TableStyle
from reportlab.lib.styles import getSampleStyleSheet
from io import BytesIO

from utils.general import create_if_not_exists


class GenomeExplorer:
    def __init__(self, config):
        """
        Initialize project runner
        """

        self.logger = logging.getLogger("genome_explorer")
        self.config = config
        self.report_path = config["report_path"]

    def get_statistics(self, samfile):
        """
        Get the total number of mapped reads in samfile
        :param samfile:
        :return:
        """
        # Count the number of reads for each chromosome.
        statistics = dict()
        chromosome_read_counts = {}
        chromosome_lengths = {}
        gc_bases = 0
        total_bases = 0
        for read in samfile:
            # Get the reference name (chromosome) of the read
            # Sum only mapped reads
            if read.is_unmapped:
                continue
            chromosome = samfile.get_reference_name(read.reference_id)
            # if chromosome is None:
            #     continue
            if chromosome not in chromosome_lengths:
                chromosome_lengths[chromosome] = samfile.get_reference_length(
                    chromosome)

            sequence = read.query_sequence
            gc_bases += sequence.count("G") + sequence.count("C")
            total_bases += len(sequence)

            # Increment the read count for the corresponding chromosome
            if chromosome in chromosome_read_counts:
                chromosome_read_counts[chromosome] += 1
            else:
                self.logger.info(f"New chromosome found: {chromosome}")
                chromosome_read_counts[chromosome] = 1

        # Calculate the GC percentage (only mapped reads)
        gc_percentage = (gc_bases / total_bases) * 100
        # Calculate the total number of mapped reads
        total_read_count = sum(chromosome_read_counts.values())
        # Calculate the total length of the genome
        total_length = sum(chromosome_lengths.values())
        statistics["total_read_count"] = total_read_count
        statistics["total_length"] = total_length
        statistics["gc_percentage"] = gc_percentage
        chromosome_read_counts_1 = {f"total_read_count - {key}": value for
                                    key, value in
                                    chromosome_read_counts.items()}
        statistics.update(chromosome_read_counts_1)
        self.logger.info(f"Total number of reads: {total_read_count}")
        self.logger.info(f"Total length of the genome: {total_length}")
        return statistics

    def mean_coverage(self, samfile, contig=None):
        """
        Calculate the mean coverage for a given contig
        :param samfile:
        :param contig:
        :return:
        """
        self.logger.info(f"Calculating coverage for contig {contig}")
        # length = samfile.get_reference_length(contig)
        # coverage = samfile.count_coverage(contig, start=0, stop=length)
        # base_sums = map(sum, coverage)
        # total_coverage = sum(base_sums)
        # average_coverage = total_coverage / length
        coverage_values = []
        for pileupcolumn in samfile.pileup(contig):
            coverage_values.append(pileupcolumn.n)

        average_coverage = sum(coverage_values) / len(coverage_values)

        return contig, average_coverage

    def get_average_coverage(self, samfile):
        self.logger.info("Calculating average coverage for each chromosome")
        chromosomes = []
        av_coverage_values = []

        # Iterate through the references (chromosomes)
        for contig in samfile.references:
            # Check if the chromosome is one of the autosomes (1-22), X, or Y
            if contig.isdigit() and 1 <= int(contig) <= 22 or contig in ["X",
                                                                         "Y"]:
                contig, contig_mean_coverage = self.mean_coverage(samfile,
                                                                  contig)

                chromosomes.append(contig)
                av_coverage_values.append(contig_mean_coverage)

        coverage_dict = dict(zip(chromosomes, av_coverage_values))
        return coverage_dict

    # def get_average_coverage_spark(self, samfile):
    #     self.logger.info("Calculating average coverage for each chromosome")
    #
    #     contigs = samfile.references
    #     # Initialize Spark
    #     spark = SparkSession.builder.appName(
    #         "CoverageCalculation").getOrCreate()
    #     # Iterate through the references (chromosomes)
    #     results = spark.sparkContext.parallelize(contigs).map(
    #         lambda contig: self.mean_coverage(samfile, contig)).collect()
    #
    #     spark.stop()
    #
    #     chromosomes = []
    #     coverage_values = []
    #
    #     for contig, average_coverage in results:
    #         chromosomes.append(contig)
    #         coverage_values.append(average_coverage)
    #     coverage_dict = dict(zip(chromosomes, coverage_values))
    #     return coverage_dict

    def generate_pdf_report(self, statistics, coverage_dict):
        # Create a PDF document
        self.logger.info("Generating PDF report")
        create_if_not_exists("reports")
        # Create a PDF document
        doc = SimpleDocTemplate(self.report_path, pagesize=landscape(letter))

        # Create a story (content) for the PDF
        story = []
        # Add a title
        title_style = getSampleStyleSheet()['Title']
        title = Paragraph("Gene Planet Data Scientist Job Assignment",
                          title_style)
        story.append(title)

        # Add text content
        text_content = """
        This is an auto generated report created by GenomeExplorer. It contains
        statistics about the bam file and a graph showing the average coverage
        across the genome for selected chromosomes. 
        """

        text_style = getSampleStyleSheet()['Normal']
        text = Paragraph(text_content, text_style)
        story.append(text)

        # Add a spacer
        story.append(Spacer(1, 12))

        for statistic, value in statistics.items():
            text_content = f"{statistic}: {value}"
            text = Paragraph(text_content, text_style)
            story.append(text)

        # Add a spacer
        story.append(Spacer(1, 12))

        # Create the graph
        chromosomes = list(coverage_dict.keys())
        av_coverage_values = list(coverage_dict.values())
        plt.figure(figsize=(10, 6))
        plt.plot(chromosomes, av_coverage_values, marker='o', linestyle='-')
        plt.title(
            "Average Coverage Across The Genome For Selected Chromosomes")
        plt.xlabel("Chromosome")
        plt.ylabel("Average Coverage")
        plt.xticks(rotation=45)
        plt.grid(True)
        plt.tight_layout()

        # Save the graph to a BytesIO buffer
        buf = BytesIO()
        plt.savefig(buf, format='png',
                    bbox_inches='tight')  # Use 'bbox_inches' to avoid cutting off labels
        buf.seek(0)

        # Create an Image object from the graph and adjust the width and height
        graph = Image(buf, width=doc.width, height=doc.height / 2)

        # Add the graph to the PDF
        story.append(Spacer(1, 12))
        story.append(
            Paragraph("Average coverage across the genome by chromosome",
                      text_style))
        story.append(Spacer(1, 6))
        story.append(graph)  # Add the graph as an Image

        # Build the PDF document
        doc.build(story)
