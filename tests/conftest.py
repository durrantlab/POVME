import os

import pytest

TEST_DIR = os.path.dirname(__file__)


@pytest.fixture
def path_4nss_input():
    return os.path.join(TEST_DIR, "files/4nss/povme-input.ini")


@pytest.fixture
def path_4nss_output():
    return os.path.join(TEST_DIR, "tmp/", "4nss/POVME_")
