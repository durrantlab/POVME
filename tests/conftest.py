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


@pytest.fixture
def path_rel1_config():
    return os.path.join(TEST_DIR, "files/rel1/pocket-id.yml")


@pytest.fixture
def path_rel1_output():
    return os.path.join(TEST_DIR, "tmp/", "rel1/")


@pytest.fixture
def path_rogfp2_pdb():
    return os.path.join(TEST_DIR, "files/rogfp2/rogfp2.pdb")


@pytest.fixture
def path_rogfp2_config():
    return os.path.join(TEST_DIR, "files/rogfp2/pocket-id.yml")


@pytest.fixture
def path_rogfp2_output():
    return os.path.join(TEST_DIR, "tmp/", "rogfp2/")
