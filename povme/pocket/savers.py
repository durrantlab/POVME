import csv
import os

from loguru import logger


def init_vol_csv(output_prefix: str) -> None:
    csv_filename = f"{output_prefix}volumes.csv"
    logger.info(f"Initializing CSV volume file: {csv_filename}")
    if not os.path.exists(csv_filename):
        with open(csv_filename, mode="w", newline="") as csv_file:
            writer = csv.writer(csv_file)
            writer.writerow(["frame_idx", "volume"])
        logger.info(f"Initialized CSV file with headers: {csv_filename}")


def write_vol_csv(results_vol: list[tuple[int, float]], output_prefix: str) -> None:
    """
    Append frame indices and volumes to a CSV file.

    Args:
        results_vol: A list of tuples, each containing (frame_idx, volume).
        output_prefix: Path to the directory where the CSV file will be saved.
    """
    csv_filename = os.path.join(output_prefix, "volumes.csv")
    file_exists = os.path.isfile(csv_filename)

    with open(csv_filename, mode="a", newline="") as csv_file:
        writer = csv.writer(csv_file)
        if not file_exists:
            writer.writerow(["frame_idx", "volume"])
        for result in results_vol:
            writer.writerow([result[0], result[1]])
