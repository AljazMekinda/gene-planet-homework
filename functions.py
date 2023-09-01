import logging
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import pysam


class GenomeExplorer:
    def __init__(self, config):
        """
        Initialize project runner
        """

        self.logger = logging.getLogger("genome_explorer")
        self.config = config
        self.base_order = config["coverage"]["base_order"]

    # def coverages_by_bases(self, samfile, contig):
    #     """
    #     Calculate average coverage for each contig in samfile
    #     and each base ACGT in this contig
    #     :param samfile:
    #     :return: {contig: {A: sum_coverage, C: sum_coverage, G:
    #                             sum_coverage, T: sum_coverage}
    #                             ...}
    #     """
    #     total_coverage = 0
    #     total_length = 0
    #
    #     coverage_dict = {}
    #
    #     self.logger.info(f"Calculating coverage for contig {contig}")
    #     coverage_dict[contig] = dict()
    #     reference_start = samfile.get_reference_start(contig)
    #     reference_end = samfile.get_reference_end(contig)
    #     length = samfile.get_reference_length(contig)
    #     total_length += length
    #     # coverage for the bases A, C, G, T (in this order)
    #     coverage = samfile.count_coverage(contig, start=reference_start,
    #                                       stop=reference_end)
    #
    #     base_sums = map(sum, coverage)
    #
    #     for i, value in enumerate(base_sums):
    #         self.logger.info(f"Coverage for base {self.base_order[i]}: "
    #                          f"{value}")
    #         coverage_dict[contig][self.base_order[i]] = value
    #
    #     total_coverage += sum(base_sums)
    #
    #     return coverage_dict, total_coverage, total_length
    #
    # def coverage_contig(self, samfile, contig, genome_length,
    #                     chunk_size=1000000):
    #     """
    #     Calculate average coverage for a selected contig in samfile.
    #     average coverage in a region is calculated as:
    #     n_mapped_reads * average_length_of_reads_in_region/total_length_of_region
    #
    #     :param samfile: samfile to calculate average coverage
    #     :param contig: contig to calculate average coverage
    #     :param genome_length: total length of the genome
    #     :param chunk_size: size of a region to calculate average coverage
    #     :return:
    #     """
    #     self.logger.info(f"Calculating average coverage for contig {contig}")
    #     # list of n_reads * average_length_of_reads_in a region (chunk)
    #     coverage_sums = []
    #     # global position
    #     pos_index = list()
    #     for i in range(0, genome_length, chunk_size):
    #         print(i)
    #         try:
    #             # try if contig is present in the region [i:i+chunk_size]
    #             contig_region = samfile.fetch(contig, i, i + chunk_size)
    #         except ValueError:
    #             pos_index.append(i)
    #             coverage_sums.append(0)
    #             continue
    #
    #         # nmb of reads in that chunk
    #         n_reads = 0
    #         lengths = np.array([])
    #         for read in contig_region:
    #             if read.is_unmapped:
    #                 continue
    #             else:
    #                 n_reads += 1
    #                 lengths = np.append(lengths, read.query_alignment_length)
    #         if n_reads > 0:
    #             coverage_sums.append(
    #                 n_reads * np.mean(lengths) / genome_length)
    #             pos_index.append(i)
    #     return coverage_sums, pos_index
    #
    # def calculate_coverages(self, samfile, total_length, chunk_size=1000000):
    #     """
    #     Visualize coverage for each contig in samfile
    #     :param samfile:
    #     :return:
    #     """
    #     contigs = samfile.header.references
    #
    #     coverage_dict = {}
    #     for contig in contigs:
    #         self.logger.info(f"Calculating coverage for contig {contig}")
    #         coverage_sums, pos_index = self.coverage_contig(samfile, contig,
    #                                                         total_length,
    #                                                         chunk_size=1000000)
    #         coverage_dict[contig] = (coverage_sums, pos_index)
    #     self.logger.info(f"Coverage for all contigs calculated")
    #     return coverage_dict

    def get_statistics(self, samfile):
        """
        Get the total number of mapped reads in samfile
        :param samfile:
        :return:
        """
        # Count the number of reads for each chromosome.
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
        self.logger.info(f"Total number of reads: {total_read_count}")
        self.logger.info(f"Total length of the genome: {total_length}")
        return (total_read_count, total_length, chromosome_read_counts,
                gc_percentage)

    def get_average_coverage(self, samfile, visualize=True):
        self.logger.info("Calculating average coverage for each chromosome")
        chromosomes = []
        coverage_values = []

        # Iterate through the references (chromosomes)
        for contig in samfile.references:
            # Check if the chromosome is one of the autosomes (1-22), X, or Y
            if contig.isdigit() and 1 <= int(contig) <= 22 or contig in ["X",
                                                                         "Y"]:
                self.logger.info(f"Calculating coverage for contig {contig}")
                length = samfile.get_reference_length(contig)
                coverage = samfile.count_coverage(contig, start=0, stop=length)
                base_sums = map(sum, coverage)
                total_coverage = sum(base_sums)
                average_coverage = total_coverage / length

                chromosomes.append(contig)
                coverage_values.append(average_coverage)
        if visualize:
            self.logger.info("Visualizing average coverage")
            self.visualize_coverage(chromosomes, coverage_values)
        coverage_dict = dict(zip(chromosomes, coverage_values))
        return coverage_dict

    def visualize_coverage(self, chromosomes, coverage_values):
        # Create a line plot
        plt.figure(figsize=(10, 6))
        plt.plot(chromosomes, coverage_values, marker='o', linestyle='-')
        plt.title("Average Coverage Across Selected Chromosomes")
        plt.xlabel("Chromosome")
        plt.ylabel("Average Coverage")
        plt.xticks(rotation=45)
        plt.grid(True)

        # Show the plot or save it to a file
        plt.tight_layout()
        plt.show()
