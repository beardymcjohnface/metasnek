import os
import warnings
import pytest
import shutil
from metasnek.fastq_finder import (
    parse_directory,
    parse_tsv_file,
    parse_samples,
    convert_to_dictionary,
    parse_samples_to_dictionary,
    write_samples_tsv,
)


@pytest.fixture(scope="session")
def temp_directory(request, tmpdir_factory):
    temp_dir = tmpdir_factory.mktemp("temp_directory")
    temp_dir_path = str(temp_dir)

    def fin():
        shutil.rmtree(temp_dir_path, ignore_errors=True)

    request.addfinalizer(fin)

    return temp_dir_path


@pytest.fixture(scope="session")
def file_list(temp_directory):
    _file_list = [
        os.path.join(str(temp_directory), "sample1_R1.fastq"),
        os.path.join(str(temp_directory), "sample1_R2.fastq"),
        os.path.join(str(temp_directory), "sample2_R1_001.fastq.gz"),
        os.path.join(str(temp_directory), "sample2_R2_001.fastq.gz"),
        os.path.join(str(temp_directory), "sample2_RS_001.fastq.gz"),
        os.path.join(str(temp_directory), "sample3.fastq"),
        os.path.join(str(temp_directory), "sample4_R2.fastq"),
        os.path.join(str(temp_directory), "sample5.fasta.gz"),
        os.path.join(str(temp_directory), "sample6.fastq.gz"),
        os.path.join(str(temp_directory), "sample7.R1.fastq"),
        os.path.join(str(temp_directory), "sample7.R2.fastq"),
        os.path.join(str(temp_directory), "sample8_1.FASTQ.GZ"),
        os.path.join(str(temp_directory), "sample8_2.FASTQ.GZ"),
        os.path.join(str(temp_directory), "sample9.1.001.fq.gz"),
        os.path.join(str(temp_directory), "sample9.2.001.fq.gz"),
    ]

    for file_path in _file_list:
        open(file_path, "w").close()

    return _file_list


def assert_parsed_files(paired_files, unpaired_files, temp_directory):
    # Assert paired files
    assert (
        "sample1",
        os.path.join(str(temp_directory), "sample1_R1.fastq"),
        os.path.join(str(temp_directory), "sample1_R2.fastq"),
        None,
    ) in paired_files
    assert (
        "sample2",
        os.path.join(str(temp_directory), "sample2_R1_001.fastq.gz"),
        os.path.join(str(temp_directory), "sample2_R2_001.fastq.gz"),
        os.path.join(str(temp_directory), "sample2_RS_001.fastq.gz"),
    ) in paired_files
    assert (
        "sample7",
        os.path.join(str(temp_directory), "sample7.R1.fastq"),
        os.path.join(str(temp_directory), "sample7.R2.fastq"),
        None,
    ) in paired_files
    assert (
        "sample8",
        os.path.join(str(temp_directory), "sample8_1.FASTQ.GZ"),
        os.path.join(str(temp_directory), "sample8_2.FASTQ.GZ"),
        None,
    ) in paired_files
    assert (
        "sample9",
        os.path.join(str(temp_directory), "sample9.1.001.fq.gz"),
        os.path.join(str(temp_directory), "sample9.2.001.fq.gz"),
        None,
    ) in paired_files
    # Assert unpaired files
    assert (
        "sample3",
        os.path.join(str(temp_directory), "sample3.fastq"),
    ) in unpaired_files
    assert (
        "sample4_R2",
        os.path.join(str(temp_directory), "sample4_R2.fastq"),
    ) in unpaired_files
    assert (
        "sample5",
        os.path.join(str(temp_directory), "sample5.fasta.gz"),
    ) in unpaired_files
    assert (
        "sample6",
        os.path.join(str(temp_directory), "sample6.fastq.gz"),
    ) in unpaired_files


@pytest.fixture(scope="session")
def create_sample_tsv(temp_directory, file_list):
    tsv_file_path = os.path.join(temp_directory, "sample.tsv")

    with open(tsv_file_path, "w") as tsv_file:
        tsv_file.write(
            "sample1\t"
            + os.path.join(str(temp_directory), "sample1_R1.fastq")
            + "\t"
            + os.path.join(str(temp_directory), "sample1_R2.fastq")
            + "\n"
            "sample2\t"
            + os.path.join(str(temp_directory), "sample2_R1_001.fastq.gz")
            + "\t"
            + os.path.join(str(temp_directory), "sample2_R2_001.fastq.gz")
            + "\t"
            + os.path.join(str(temp_directory), "sample2_RS_001.fastq.gz")
            + "\n"
            "sample3\t" + os.path.join(str(temp_directory), "sample3.fastq") + "\n"
            "sample4_R2\t"
            + os.path.join(str(temp_directory), "sample4_R2.fastq")
            + "\n"
            "sample5\t" + os.path.join(str(temp_directory), "sample5.fasta.gz") + "\n"
            "sample6\t" + os.path.join(str(temp_directory), "sample6.fastq.gz") + "\n"
            "sample7\t"
            + os.path.join(str(temp_directory), "sample7.R1.fastq")
            + "\t"
            + os.path.join(str(temp_directory), "sample7.R2.fastq")
            + "\n"
            "sample8\t"
            + os.path.join(str(temp_directory), "sample8_1.FASTQ.GZ")
            + "\t"
            + os.path.join(str(temp_directory), "sample8_2.FASTQ.GZ")
            + "\n"
            "sample9\t"
            + os.path.join(str(temp_directory), "sample9.1.001.fq.gz")
            + "\t"
            + os.path.join(str(temp_directory), "sample9.2.001.fq.gz")
        )

    return tsv_file_path


@pytest.fixture(scope="session")
def expected_dictionary(temp_directory, file_list):
    directory_dictionary = {
        "sample1": {
            "R1": os.path.join(temp_directory, "sample1_R1.fastq"),
            "R2": os.path.join(temp_directory, "sample1_R2.fastq"),
            "S": None,
        },
        "sample2": {
            "R1": os.path.join(temp_directory, "sample2_R1_001.fastq.gz"),
            "R2": os.path.join(temp_directory, "sample2_R2_001.fastq.gz"),
            "S": os.path.join(temp_directory, "sample2_RS_001.fastq.gz"),
        },
        "sample3": {
            "R1": os.path.join(temp_directory, "sample3.fastq"),
            "R2": None,
            "S": None,
        },
        "sample4_R2": {
            "R1": os.path.join(temp_directory, "sample4_R2.fastq"),
            "R2": None,
            "S": None,
        },
        "sample5": {
            "R1": os.path.join(temp_directory, "sample5.fasta.gz"),
            "R2": None,
            "S": None,
        },
        "sample6": {
            "R1": os.path.join(temp_directory, "sample6.fastq.gz"),
            "R2": None,
            "S": None,
        },
        "sample7": {
            "R1": os.path.join(str(temp_directory), "sample7.R1.fastq"),
            "R2": os.path.join(str(temp_directory), "sample7.R2.fastq"),
            "S": None,
        },
        "sample8": {
            "R1": os.path.join(str(temp_directory), "sample8_1.FASTQ.GZ"),
            "R2": os.path.join(str(temp_directory), "sample8_2.FASTQ.GZ"),
            "S": None,
        },
        "sample9": {
            "R1": os.path.join(str(temp_directory), "sample9.1.001.fq.gz"),
            "R2": os.path.join(str(temp_directory), "sample9.2.001.fq.gz"),
            "S": None,
        },
    }
    return directory_dictionary


@pytest.mark.filterwarnings("ignore:Orphaned paired")
def test_parse_directory(file_list, temp_directory):
    paired_files, unpaired_files = parse_directory(file_list)
    assert_parsed_files(paired_files, unpaired_files, temp_directory)


def test_parse_directory_orphand_r():
    # Assert warning for orphaned paired read
    with warnings.catch_warnings(record=True) as warning_list:
        parse_directory(["sample7_R2.fastq"])
        assert len(warning_list) == 1
        assert "Possible orphaned paired read detected for" in str(
            warning_list[0].message
        )


def test_parse_tsv_file(create_sample_tsv, temp_directory):
    paired_reads, unpaired_reads = parse_tsv_file(create_sample_tsv)
    assert_parsed_files(paired_reads, unpaired_reads, temp_directory)

    # Assert file not found error for missing R1 file
    invalid_tsv_file = os.path.join(temp_directory, "invalid_reads.tsv")
    with open(invalid_tsv_file, "w") as f:
        f.write("sample1\tnonexistent_R1.fastq\t")

    with pytest.raises(FileNotFoundError):
        parse_tsv_file(str(invalid_tsv_file))

    # Assert file not found error for missing R2 file
    invalid_tsv_file = os.path.join(temp_directory, "invalid_reads.tsv")
    with open(invalid_tsv_file, "w") as f:
        f.write(
            "sample1\t"
            + os.path.join(str(temp_directory), "sample1_R1.fastq")
            + "\tnonexistent_R2.fastq"
        )

    with pytest.raises(FileNotFoundError):
        parse_tsv_file(str(invalid_tsv_file))


@pytest.mark.filterwarnings("ignore:Orphaned paired")
def test_parse_samples(temp_directory, create_sample_tsv):
    paired_reads_dir, unpaired_reads_dir = parse_samples(temp_directory)
    assert_parsed_files(paired_reads_dir, unpaired_reads_dir, temp_directory)

    # Test parsing a TSV file
    paired_reads_tsv, unpaired_reads_tsv = parse_samples(create_sample_tsv)
    assert_parsed_files(paired_reads_tsv, unpaired_reads_tsv, temp_directory)

    # Test parsing an invalid file or directory
    invalid_path = os.path.join(temp_directory, "nonexistent_file.tsv")
    with pytest.raises(ValueError):
        parse_samples(invalid_path)

    invalid_path = os.path.join(temp_directory, "nonexistent_directory")
    with pytest.raises(ValueError):
        parse_samples(invalid_path)

    # Test parsing a TSV file with FileNotFoundError
    invalid_tsv_file = os.path.join(temp_directory, "invalid_reads.tsv")
    invalid_tsv_content = "sample1\tnonexistent_R1.fastq\t"

    with open(invalid_tsv_file, "w") as tsv_file:
        tsv_file.write(invalid_tsv_content)

    with pytest.raises(ValueError):
        parse_samples(invalid_tsv_file)


def test_convert_to_dictionary():
    paired_reads = {
        ("sample1", "sample1_R1.fastq", "sample1_R2.fastq", None),
        ("sample2", "sample2_R1.fastq", "sample2_R2.fastq", "sample2_RS.fastq"),
    }
    unpaired_reads = {("sample3", "sample3_R1.fastq"), ("sample4", "sample4_R1.fastq")}

    expected_dictionary = {
        "sample1": {"R1": "sample1_R1.fastq", "R2": "sample1_R2.fastq", "S": None},
        "sample2": {
            "R1": "sample2_R1.fastq",
            "R2": "sample2_R2.fastq",
            "S": "sample2_RS.fastq",
        },
        "sample3": {"R1": "sample3_R1.fastq", "R2": None, "S": None},
        "sample4": {"R1": "sample4_R1.fastq", "R2": None, "S": None},
    }

    assert convert_to_dictionary(paired_reads, unpaired_reads) == expected_dictionary


@pytest.mark.filterwarnings("ignore:Orphaned paired")
def test_parse_samples_to_dictionary(
    temp_directory, create_sample_tsv, expected_dictionary
):
    assert parse_samples_to_dictionary(temp_directory) == expected_dictionary

    assert parse_samples_to_dictionary(create_sample_tsv) == expected_dictionary

    # Test parsing an invalid file or directory
    invalid_path = os.path.join(temp_directory, "nonexistent_file.tsv")
    with pytest.raises(ValueError):
        parse_samples_to_dictionary(invalid_path)

    invalid_path = os.path.join(temp_directory, "nonexistent_directory")
    with pytest.raises(ValueError):
        parse_samples_to_dictionary(invalid_path)

    # Test parsing a TSV file with FileNotFoundError
    invalid_tsv_file = os.path.join(temp_directory, "invalid_reads.tsv")
    invalid_tsv_content = "sample1\tnonexistent_R1.fastq\t"

    with open(invalid_tsv_file, "w") as tsv_file:
        tsv_file.write(invalid_tsv_content)

    with pytest.raises(ValueError):
        parse_samples_to_dictionary(invalid_tsv_file)


def test_write_samples_tsv(temp_directory):
    samples_dictionary = {
        "sample1": {
            "R1": os.path.join(temp_directory, "sample1_R1.fastq"),
            "R2": os.path.join(temp_directory, "sample1_R2.fastq"),
        },
        "sample2": {
            "R1": os.path.join(temp_directory, "sample2_R1_001.fastq.gz"),
            "R2": os.path.join(temp_directory, "sample2_R2_001.fastq.gz"),
            "S": os.path.join(temp_directory, "sample2_RS_001.fastq.gz"),
        },
        "sample3": {"R1": os.path.join(temp_directory, "sample3.fastq"), "R2": None},
        "sample4_R2": {
            "R1": os.path.join(temp_directory, "sample4_R2.fastq"),
            "R2": None,
        },
        "sample5": {"R1": os.path.join(temp_directory, "sample5.fasta.gz"), "R2": None},
        "sample6": {"R1": os.path.join(temp_directory, "sample6.fastq.gz"), "R2": None},
    }

    output_file = os.path.join(temp_directory, "output.tsv")
    write_samples_tsv(samples_dictionary, output_file)

    with open(output_file, "r") as file:
        content = file.read()

    expected_content = (
        "sample1\t"
        + os.path.join(temp_directory, "sample1_R1.fastq")
        + "\t"
        + os.path.join(temp_directory, "sample1_R2.fastq")
        + "\n"
        "sample2\t"
        + os.path.join(temp_directory, "sample2_R1_001.fastq.gz")
        + "\t"
        + os.path.join(temp_directory, "sample2_R2_001.fastq.gz")
        + "\t"
        + os.path.join(temp_directory, "sample2_RS_001.fastq.gz")
        + "\n"
        "sample3\t" + os.path.join(temp_directory, "sample3.fastq") + "\n"
        "sample4_R2\t" + os.path.join(temp_directory, "sample4_R2.fastq") + "\n"
        "sample5\t" + os.path.join(temp_directory, "sample5.fasta.gz") + "\n"
        "sample6\t" + os.path.join(temp_directory, "sample6.fastq.gz") + "\n"
    )

    assert content == expected_content
