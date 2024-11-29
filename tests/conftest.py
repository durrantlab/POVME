import os

import pytest

from povme import enable_logging

TEST_DIR = os.path.dirname(__file__)

@pytest.fixture(scope="session", autouse=True)
def turn_on_logging():
    enable_logging(10)


@pytest.fixture
def path_4nss_config():
    return os.path.join(TEST_DIR, "files/4nss/povme.yml")


@pytest.fixture
def path_4nss_output():
    return os.path.join(TEST_DIR, "tmp/", "4nss/POVME_")
