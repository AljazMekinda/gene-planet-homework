

def get_coverages(samfile, base_order):
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
        total_coverage_dict[contig] = dict()
        length = samfile.get_reference_length(contig)
        total_length += length
        # coverage for the bases A, C, G, T (in this order)
        coverage = samfile.count_coverage(contig, start=0, stop=length)
        base_sums = map(sum, coverage)

        for i, value in enumerate(base_sums):
            total_coverage_dict[contig][base_order[i]] = value

        total_coverage += sum(base_sums)


    return total_coverage_dict, total_coverage, total_length