import pandas as pd
from task2.utils import count_unique_icgc_mutations, get_max_and_min_icgc_mutation_count, get_sorted_mutation_count


def test_count_unique_icgc_mutations():
    """count_unique_icgc_mutations should accept a dataframe that contains columns: "mutated_from_allele", "mutated_to_allele", "icgc_mutation_id"
    and return the number of unique icgc_mutation_id records per each combination of mutated_from_allele and mutated_to_allele
    """
    data = pd.DataFrame(
        [
            ["A", "B", "MU1"],
            ["A", "B", "MU2"],
            ["A", "B", "MU2"],
            ["C", "D", "MU3"],
            ["C", "D", "MU4"],
            ["C", "D", "MU5"],
            ["F", "G", "MU6"]
        ],
        columns=["mutated_from_allele", "mutated_to_allele", "icgc_mutation_id"]
    )
    expected_result = [("A", "B", 2), ("C", "D", 3), ("F", "G", 1)]

    result = count_unique_icgc_mutations(data)
    for mutated_from, mutated_to, expected_count in expected_result:
        uniq_mutation_count = result[(result.mutated_from_allele == mutated_from) & (result.mutated_to_allele == mutated_to)]
        # Records are grouped by mutation pairs, so each mutation pair must have only one record in the aggregated DF
        assert len(uniq_mutation_count.uniq_count_icgc_mutation.values) == 1
        # Confirm that the unique count matches static expected value
        assert uniq_mutation_count.uniq_count_icgc_mutation.values[0] == expected_count


def test_get_sorted_mutation_count():
    """get_sorted_mutation_count accepts the source dataframe and returns unique mutation counts per each sample id,
    sorted by counts and alphanumerically by sample id
    """
    data = pd.DataFrame(
        [
            ["SA1", "MU1"],
            ["SA1", "MU2"],
            ["SA2", "MU3"],
            ["SA2", "MU4"],
            ["SA2", "MU5"],
            ["SA3", "MU6"],
            ["SA4", "MU7"]
        ],
        columns=["icgc_sample_id", "icgc_mutation_id"]
    )
    expected_result = [("SA2", 3), ("SA1", 2), ("SA3", 1), ("SA4", 1)]

    sorted_mutation_count = get_sorted_mutation_count(data)
    for sample_id, mutation_count in expected_result:
        selection = sorted_mutation_count[sorted_mutation_count.icgc_sample_id == sample_id]
        # Each sample id should have one record in the aggregated DF
        assert len(selection) == 1
        # Confirm that the mutation count matches static expected value
        assert selection.uniq_count_icgc_mutation.values[0] == mutation_count


def test_get_max_and_min_icgc_mutation_count():
    """get_max_and_min_icgc_mutation_count accepts a dataframe that contains columns: "icgc_sample_id", "icgc_mutation_id"
    and returns the number of unique icgc_mutation_id records per each icgc_sample_id
    """
    data = pd.DataFrame(
        [
            ["SA1", "MU1"],
            ["SA1", "MU2"],
            ["SA2", "MU3"],
            ["SA2", "MU4"],
            ["SA2", "MU5"],
            ["SA3", "MU6"]
        ],
        columns=["icgc_sample_id", "icgc_mutation_id"]
    )
    expected_max = ("SA2", 3)
    expected_min = ("SA3", 1)

    max_record, min_record = get_max_and_min_icgc_mutation_count(data)
    # Check max record
    assert max_record["icgc_sample_id"] == expected_max[0]
    assert max_record["mutation_count"] == expected_max[1]
    # Check min record
    assert min_record["icgc_sample_id"] == expected_min[0]
    assert min_record["mutation_count"] == expected_min[1]


def test_get_max_and_min_icgc_mutation_count_alphanumeric_ordering():
    """If more than one icgc_sample_id-s have highest/lowest number of mutations,
    then return the one with higher alphanumeric precedence
    """
    data = pd.DataFrame(
        [
            ["SA1", "MU1"],
            ["SA2", "MU2"],
            ["SA3", "MU3"],
            ["SA3", "MU4"],
            ["SA4", "MU5"],
            ["SA4", "MU6"],
        ],
        columns=["icgc_sample_id", "icgc_mutation_id"]
    )
    expected_max_sample_id = "SA3"
    expected_min_sample_id = "SA1"

    max_record, min_record = get_max_and_min_icgc_mutation_count(data)
    # Check max record
    assert max_record["icgc_sample_id"] == expected_max_sample_id
    # Check min record
    assert min_record["icgc_sample_id"] == expected_min_sample_id
