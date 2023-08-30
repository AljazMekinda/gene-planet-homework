import logging

class GenomeExplorer:
    def __init__(self, config):
        """
        Initialize project runner
        """

        self.logger = logging.getLogger("genome_explorer")
        self.config = config
        self.base_order = config["coverage"]["base_order"]

    def get_coverages(self, samfile):
        """
        Calculate average coverage for each contig in samfile
        and each base ACGT in this contig
        :param samfile:
        :return: {contig: {A: sum_coverage, C: sum_coverage, G:
                                sum_coverage, T: sum_coverage}
                                ...}
        """
        total_coverage = 0
        total_length = 0

        total_coverage_dict = {}
        contigs = samfile.header.references

        for contig in contigs:
            self.logger.info(f"Calculating coverage for contig {contig}")
            total_coverage_dict[contig] = dict()
            length = samfile.get_reference_length(contig)
            total_length += length
            # coverage for the bases A, C, G, T (in this order)
            coverage = samfile.count_coverage(contig, start=0, stop=length)
            base_sums = map(sum, coverage)

            for i, value in enumerate(base_sums):
                self.logger.info(f"Coverage for base {self.base_order[i]}: "
                                 f"{value}")
                total_coverage_dict[contig][self.base_order[i]] = value

            total_coverage += sum(base_sums)

        return total_coverage_dict, total_coverage, total_length