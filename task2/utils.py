import pandas as pd


def count_unique_icgc_mutations(data: pd.DataFrame):
    """Count unique mutation ids for all possible combinations of mutated_from_allele and mutated_to_allele.
    """
    # Create new DF with unique mutation counts
    data = data.groupby(["mutated_from_allele", "mutated_to_allele"])["icgc_mutation_id"].nunique().reset_index()
    # Make the name of the column with unique mutation counts more descriptive
    return data.rename({"icgc_mutation_id": "count_unique_icgc_mutation_ids"}, axis=1)


def get_sorted_mutation_count(data: pd.DataFrame):
    """Count unique icgc_mutation_id-s per each icgc_sample_id and sort records by count.
    For sample ids with equal number of mutations, sort alphanumerically
    """
    # Create new DF with unique mutation counts
    data = data.groupby(["icgc_sample_id"])["icgc_mutation_id"].nunique().reset_index()
    # Sort values by unique mutation counts and by icgc_sample_id (alphanumerically)
    data = data.sort_values(["icgc_mutation_id", "icgc_sample_id"], ascending=[False, True])
    # Reset index after sorting and make the name of the column with unique mutation counts more descriptive
    return data.reset_index(drop=True).rename({"icgc_mutation_id": "uniq_icgc_mutation_count"}, axis=1)


def get_max_and_min_icgc_mutation_count(data: pd.DataFrame):
    """Get icgc_sample_id with maximum and minimum unique icgc_mutation_id count.
    If more than one sample ids have highest/lowest number of mutations,
    a sample id with higher alphanumeric precedence is returned
    """
    icgc_sorted_count = get_sorted_mutation_count(data)

    max_count = icgc_sorted_count.head(1).uniq_icgc_mutation_count.values[0]
    max_count_sample_id = icgc_sorted_count[icgc_sorted_count.uniq_icgc_mutation_count == max_count].head(1).icgc_sample_id.values[0]

    min_count= icgc_sorted_count.tail(1).uniq_icgc_mutation_count.values[0]
    min_count_sample_id = icgc_sorted_count[icgc_sorted_count.uniq_icgc_mutation_count == min_count].head(1).icgc_sample_id.values[0]

    return (
        {"icgc_sample_id": max_count_sample_id, "mutation_count": max_count},
        {"icgc_sample_id": min_count_sample_id, "mutation_count": min_count},
    )
