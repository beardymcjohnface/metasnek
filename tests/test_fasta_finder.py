import os
import warnings
import pytest
import shutil

from metasnek.fastq_finder import (
    fastas_from_directory,
    parse_tsv_file,
    parse_fastas,
    write_fastas_tsv,
    combine_fastas,
)


@pytest.fixture
def dir_test_files():
    temp_dir = tempfile.mkdtemp()
    files = [
        "sequence.fasta",
        "protein.fa",
        "data.fna",
        "output.ffn",
        "results.faa",
        "report.frn",
        "data.txt",
        "code.py",
        "image.png",
    ]
    for file_name in files:
        open(os.path.join(temp_dir, file_name), 'w').close()
    yield temp_dir
    shutil.rmtree(temp_dir)


@pytest.fixture
def expected_dictionary(dir_test_files):
    expected_result = {
        'sequence.fasta': os.path.join(dir_test_files, 'sequence.fasta'),
        'protein.fa': os.path.join(dir_test_files, 'protein.fa'),
        'data.fna': os.path.join(dir_test_files, 'data.fna'),
        'output.ffn': os.path.join(dir_test_files, 'output.ffn'),
        'results.faa': os.path.join(temp_dir, 'results.faa'),
        'report.frn': os.path.join(temp_dir, 'report.frn'),
    }
    return expected_result


def test_fastas_from_directory_empty():
    temp_dir = tempfile.mkdtemp()
    assert fastas_from_directory(temp_dir) == {}


def test_fastas_from_directory_mixed_formats(dir_test_files, expected_dictionary):
    fasta_files = fastas_from_directory(dir_test_files)
    assert fasta_files == expected_dictionary


@pytest.fixture
def tsv_file_path():
    content = "ref1\tfile1.fasta\nref2\tfile2.fasta\nref3\tfile3.fasta"
    with tempfile.NamedTemporaryFile(mode='w', delete=False) as temp_file:
        temp_file.write(content)
        temp_file_path = temp_file.name
    yield temp_file_path
    os.remove(temp_file_path)


def test_parse_tsv_valid_file(tsv_file_path):
    expected_output = {
        "ref1": "file1.fasta",
        "ref2": "file2.fasta",
        "ref3": "file3.fasta"
    }
    assert parse_tsv_file(tsv_file_path) == expected_output


def test_parse_tsv_empty_file():
    with tempfile.NamedTemporaryFile(mode='w', delete=False) as temp_file:
        temp_file_path = temp_file.name
    assert parse_tsv_file(temp_file_path) == {}
    os.remove(temp_file_path)


def test_parse_tsv_file_not_found():
    with pytest.raises(FileNotFoundError):
        parse_tsv_file("non_existent_file.tsv")


# def test_parse_tsv_malformed_lines(tsv_file_path):
#     with open(tsv_file_path, 'a') as tsv_file:
#         tsv_file.write("ref4")
#     with pytest.raises(ValueError):
#         parse_tsv_file(tsv_file_path)


@pytest.fixture
def tsv_file_path():
    content = "ref1\tfile1.fasta\nref2\tfile2.fasta\nref3\tfile3.fasta"
    with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix=".tsv") as temp_file:
        temp_file.write(content)
        temp_file_path = temp_file.name
    yield temp_file_path
    os.remove(temp_file_path)


@pytest.fixture
def fasta_directory(tmpdir):
    fasta_content = ">seq1\nACGT\n>seq2\nTGCA"
    tsv_content = "ref1\tfile1.fasta\nref2\tfile2.fasta\nref3\tfile3.fasta"

    tmpdir.mkdir("fasta_files")
    fasta_file = tmpdir.join("fasta_files", "test.fasta")
    tsv_file = tmpdir.join("test.tsv")

    fasta_file.write(fasta_content)
    tsv_file.write(tsv_content)

    return tmpdir


def test_parse_fastas_valid_fasta(fasta_file_path):
    expected_output = {"test": fasta_file_path}
    assert parse_fastas(fasta_file_path) == expected_output


def test_parse_fastas_valid_tsv(tsv_file_path):
    expected_output = {
        "ref1": "file1.fasta",
        "ref2": "file2.fasta",
        "ref3": "file3.fasta"
    }
    assert parse_fastas(tsv_file_path) == expected_output


def test_parse_fastas_valid_directory(fasta_directory):
    expected_output = {
        "test": str(fasta_directory.join("fasta_files", "test.fasta"))
    }
    assert parse_fastas(str(fasta_directory)) == expected_output


def test_parse_fastas_invalid_input():
    assert parse_fastas("non_existent_path") == {}


def test_write_fastas_tsv():
    fasta_dict = {
        "ref1": "file1.fasta",
        "ref2": "file2.fasta",
        "ref3": "file3.fasta"
    }
    expected_content = "ref1\tfile1.fasta\nref2\tfile2.fasta\nref3\tfile3.fasta\n"

    with tempfile.NamedTemporaryFile(mode='w', delete=False) as temp_file:
        temp_file_path = temp_file.name
        write_fastas_tsv(fasta_dict, temp_file_path)

        with open(temp_file_path, 'r') as tsv_file:
            content = tsv_file.read()

        assert content == expected_content

    os.remove(temp_file_path)


def test_write_fastas_tsv_empty():
    fasta_dict = {}
    expected_content = ""

    with tempfile.NamedTemporaryFile(mode='w', delete=False) as temp_file:
        temp_file_path = temp_file.name
        write_fastas_tsv(fasta_dict, temp_file_path)

        with open(temp_file_path, 'r') as tsv_file:
            content = tsv_file.read()

        assert content == expected_content

    os.remove(temp_file_path)


@pytest.fixture
def fasta_files(tmpdir):
    fasta_content_1 = ">seq1\nACGT\n>seq2\nTGCA"
    fasta_content_2 = ">seq3\nGGGG\n>seq4\nCCCC"

    fasta_file_1 = tmpdir.join("file1.fasta")
    fasta_file_2 = tmpdir.join("file2.fasta")

    fasta_file_1.write(fasta_content_1)
    fasta_file_2.write(fasta_content_2)

    return {
        "ref1": str(fasta_file_1),
        "ref2": str(fasta_file_2)
    }

def test_combine_fastas(fasta_files):
    expected_content = (
        ">ref1:seq1\nACGT\n>ref1:seq2\nTGCA\n"
        ">ref2:seq3\nGGGG\n>ref2:seq4\nCCCC\n"
    )

    with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix=".fasta") as temp_file:
        temp_file_path = temp_file.name
        combine_fastas(fasta_files, temp_file_path)

        with open(temp_file_path, 'r') as combined_fasta:
            content = combined_fasta.read()

        assert content == expected_content

    os.remove(temp_file_path)