import argparse

from povme.config import PocketDetectConfig, PocketVolumeConfig
from povme.pocket.detect import PocketDetector
from povme.pocket.volume import PocketVolume


def detect_pockets(args):
    """Handler for the 'detect' subcommand."""
    config = PocketDetectConfig()
    config.from_yaml(args.config)
    detector = PocketDetector(config)
    detector.run(path_pdb=args.input_pdb, output_prefix=args.output_prefix)


def calculate_volume(args):
    """Handler for the 'volume' subcommand."""
    config = PocketVolumeConfig()
    config.from_yaml(args.config)
    volume_calculator = PocketVolume(config)
    volume_calculator.run(
        path_pdb=args.input_pdb,
        output_prefix=args.output_prefix,
        chunk_size=args.chunk_size,
    )


def povme_cli():
    parser = argparse.ArgumentParser(
        description="POVME: Pocket detection and volume estimation tool."
    )
    subparsers = parser.add_subparsers(
        title="Commands", description="Available subcommands", dest="command"
    )

    # Subcommand: detect
    detect_parser = subparsers.add_parser(
        "detect", help="Detect pockets in a given PDB file."
    )
    detect_parser.add_argument(
        "-c", "--config", required=True, help="Path to the configuration YAML file."
    )
    detect_parser.add_argument(
        "-i", "--input-pdb", required=True, help="Path to the input PDB file."
    )
    detect_parser.add_argument(
        "-o",
        "--output-prefix",
        required=False,
        help="Prefix for the output files.",
        default="./",
    )
    detect_parser.set_defaults(func=detect_pockets)

    # Subcommand: volume
    volume_parser = subparsers.add_parser(
        "volume", help="Calculate pocket volumes from a PDB file."
    )
    volume_parser.add_argument(
        "-c", "--config", required=True, help="Path to the configuration YAML file."
    )
    volume_parser.add_argument(
        "-i", "--input-pdb", required=True, help="Path to the input PDB file."
    )
    volume_parser.add_argument(
        "-o",
        "--output-prefix",
        required=False,
        help="Prefix for the output files.",
        default="./",
    )
    volume_parser.add_argument(
        "--chunk-size",
        type=int,
        default=10,
        help="Number of frames to process at once.",
    )
    volume_parser.set_defaults(func=calculate_volume)

    # Parse the arguments and call the appropriate handler
    args = parser.parse_args()
    if args.command:
        args.func(args)
    else:
        parser.print_help()
