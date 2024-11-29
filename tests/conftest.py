import os

import pytest

TEST_DIR = os.path.dirname(__file__)


@pytest.fixture
def path_4nss_config():
    return os.path.join(TEST_DIR, "files/4nss/povme.yml")


@pytest.fixture
def path_4nss_output():
    return os.path.join(TEST_DIR, "tmp/", "4nss/POVME_")
