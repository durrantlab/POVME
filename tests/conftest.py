import os

import pytest

from povme import enable_logging

TEST_DIR = os.path.dirname(__file__)

try:
    import ray

    HAS_RAY = True
except ImportError:
    HAS_RAY = False


def pytest_sessionstart(session):  # pytest_configure(config)
    """Called after the Session object has been created and
    before performing collection and entering the run test loop.
    """

    # Creates a tmp directory for writing files.
    os.makedirs("./tests/tmp", exist_ok=True)


def pytest_sessionfinish(session, exitstatus):
    # Shutdown ray if possible.
    if HAS_RAY:
        ray.shutdown()


@pytest.fixture(scope="session", autouse=True)
def turn_on_logging():
    enable_logging(10)


@pytest.fixture
def path_4nss_config():
    return os.path.join(TEST_DIR, "files/4nss/povme.yml")


@pytest.fixture
def path_4nss_output():
    return os.path.join(TEST_DIR, "tmp/", "4nss/")


@pytest.fixture
def path_rel1_config():
    return os.path.join(TEST_DIR, "files/rel1/pocket-id.yml")


@pytest.fixture
def path_rel1_output():
    return os.path.join(TEST_DIR, "tmp/", "rel1/")


@pytest.fixture
def path_5csn_config():
    return os.path.join(TEST_DIR, "files/5csn/pocket-id.yml")


@pytest.fixture
def path_5csn_output():
    return os.path.join(TEST_DIR, "tmp/", "5csn/")


@pytest.fixture
def path_rogfp2_pdb():
    return os.path.join(TEST_DIR, "files/rogfp2/rogfp2.pdb")


@pytest.fixture
def path_rogfp2_config():
    return os.path.join(TEST_DIR, "files/rogfp2/pocket-id.yml")


@pytest.fixture
def path_rogfp2_output():
    return os.path.join(TEST_DIR, "tmp/", "rogfp2/")


@pytest.fixture
def path_rogfp2_traj_config():
    return os.path.join(TEST_DIR, "files/rogfp2/povme.yml")


@pytest.fixture
def path_rogfp2_traj_output():
    return os.path.join(TEST_DIR, "tmp/", "rogfp2-traj/")
