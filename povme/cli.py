#!/usr/bin/env python3

import argparse
import sys

from . import __version__
from .pocket_id import run_pocket_id
from .povme import RunPOVME


def parse_inclusion_regions(args):
    inclusion_regions = []
    if args.points_inclusion_sphere:
        for sphere in args.points_inclusion_sphere:
            inclusion_regions.append(
                f"PointsInclusionSphere {' '.join(map(str, sphere))}"
            )
    if args.points_inclusion_box:
        for box in args.points_inclusion_box:
            inclusion_regions.append(f"PointsInclusionBox {' '.join(map(str, box))}")
    return inclusion_regions


def parse_exclusion_regions(args):
    exclusion_regions = []
    if args.points_exclusion_sphere:
        for sphere in args.points_exclusion_sphere:
            exclusion_regions.append(
                f"PointsExclusionSphere {' '.join(map(str, sphere))}"
            )
    if args.points_exclusion_box:
        for box in args.points_exclusion_box:
            exclusion_regions.append(f"PointsExclusionBox {' '.join(map(str, box))}")
    return exclusion_regions


def povme_cli():
    parser = argparse.ArgumentParser(
        description=f"POVME {__version__}",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    subparsers = parser.add_subparsers(
        title="subcommands",
        description="Valid subcommands",
        help="Additional help",
        dest="command",
    )
    subparsers.required = True

    # Subparser for run_povme
    povme_parser = subparsers.add_parser(
        "run_povme", help="Run the primary POVME analysis"
    )

    # Required arguments
    povme_parser.add_argument(
        "--pdb",
        required=True,
        help="Path to the receptor PDB file (can be a trajectory with multiple frames)",
    )

    # Inclusion regions
    povme_parser.add_argument(
        "--inclusion-sphere",
        type=float,
        nargs=4,
        action="append",
        metavar=("X", "Y", "Z", "RADIUS"),
        help="Define an inclusion sphere with X Y Z coordinates and RADIUS. Can be used multiple times.",
    )
    povme_parser.add_argument(
        "--inclusion-box",
        type=float,
        nargs=6,
        action="append",
        metavar=("X", "Y", "Z", "DX", "DY", "DZ"),
        help="Define an inclusion box with X Y Z coordinates and dimensions DX DY DZ. Can be used multiple times.",
    )

    # Exclusion regions
    povme_parser.add_argument(
        "--exclusion-sphere",
        type=float,
        nargs=4,
        action="append",
        metavar=("X", "Y", "Z", "RADIUS"),
        help="Define an exclusion sphere with X Y Z coordinates and RADIUS. Can be used multiple times.",
    )
    povme_parser.add_argument(
        "--exclusion-box",
        type=float,
        nargs=6,
        action="append",
        metavar=("X", "Y", "Z", "DX", "DY", "DZ"),
        help="Define an exclusion box with X Y Z coordinates and dimensions DX DY DZ. Can be used multiple times.",
    )

    # Grid spacing
    povme_parser.add_argument(
        "--grid-spacing",
        type=float,
        default=1.0,
        help="Distance separating each equidistant point (default: 1.0 Å)",
    )

    # Save points
    povme_parser.add_argument(
        "--save-points",
        type=str,
        choices=["true", "false"],
        default="false",
        help="Whether to save the generated point field (true/false)",
    )

    # Distance cutoff
    povme_parser.add_argument(
        "--distance-cutoff",
        type=float,
        default=1.09,
        help="Distance cutoff for removing points near receptor atoms (default: 1.09 Å)",
    )

    # Convex hull exclusion
    povme_parser.add_argument(
        "--convex-hull-exclusion",
        type=str,
        choices=["true", "false"],
        default="false",
        help="Activate convex hull exclusion (true/false)",
    )

    # Contiguous points criteria
    povme_parser.add_argument(
        "--contiguous-points-criteria",
        type=int,
        default=3,
        help="Minimum number of neighboring points to retain a point (default: 3)",
    )

    # Contiguous pocket seed
    povme_parser.add_argument(
        "--pocket-seed-sphere",
        type=float,
        nargs=4,
        action="append",
        metavar=("X", "Y", "Z", "RADIUS"),
        help="Define a pocket seed sphere with X Y Z coordinates and RADIUS. Can be used multiple times.",
    )
    povme_parser.add_argument(
        "--pocket-seed-box",
        type=float,
        nargs=6,
        action="append",
        metavar=("X", "Y", "Z", "DX", "DY", "DZ"),
        help="Define a pocket seed box with X Y Z coordinates and dimensions DX DY DZ. Can be used multiple times.",
    )

    # Additional POVME parameters
    povme_parser.add_argument(
        "--num-processors",
        type=int,
        default=1,
        help="Number of processors to use (default: 1)",
    )
    povme_parser.add_argument(
        "--use-disk-not-memory",
        type=str,
        choices=["true", "false"],
        default="false",
        help="Use disk space instead of memory (true/false)",
    )
    povme_parser.add_argument(
        "--output-prefix",
        type=str,
        default="./POVME_output/",
        help="Prefix for output files",
    )
    povme_parser.add_argument(
        "--save-individual-volumes",
        type=str,
        choices=["true", "false"],
        default="false",
        help="Save individual pocket volumes (true/false)",
    )
    povme_parser.add_argument(
        "--save-trajectory",
        type=str,
        choices=["true", "false"],
        default="false",
        help="Save pocket volumes as a trajectory (true/false)",
    )
    povme_parser.add_argument(
        "--equal-num-points",
        type=str,
        choices=["true", "false"],
        default="false",
        help="Ensure equal number of points per frame (true/false)",
    )
    povme_parser.add_argument(
        "--save-tabbed",
        type=str,
        choices=["true", "false"],
        default="false",
        help="Save volumes in tabbed format (true/false)",
    )
    povme_parser.add_argument(
        "--save-density-map",
        type=str,
        choices=["true", "false"],
        default="false",
        help="Save volumetric density map (true/false)",
    )
    povme_parser.add_argument(
        "--compress-output",
        type=str,
        choices=["true", "false"],
        default="false",
        help="Compress output files (true/false)",
    )

    # Subparser for pocket_id
    pocket_id_parser = subparsers.add_parser(
        "pocket_id",
        help="Run Pocket ID to identify protein pockets and generate inclusion spheres",
    )

    # Pocket ID specific arguments
    pocket_id_parser.add_argument(
        "--input-pdb", required=True, help="Path to the input PDB file for Pocket ID"
    )
    pocket_id_parser.add_argument(
        "--output-prefix",
        type=str,
        default="./PocketID_output/",
        help="Prefix for Pocket ID output files",
    )
    # Add more Pocket ID specific arguments here as needed

    args = parser.parse_args()

    if args.command == "run_povme":
        # Build the argument list for POVME
        povme_args = ["RunPOVME"]

        povme_args.append(f"PDBFileName {args.pdb}")

        # Inclusion regions
        inclusion_regions = parse_inclusion_regions(args)
        povme_args.extend(inclusion_regions)

        # Exclusion regions
        exclusion_regions = parse_exclusion_regions(args)
        povme_args.extend(exclusion_regions)

        # Other parameters
        povme_args.append(f"GridSpacing {args.grid_spacing}")
        povme_args.append(f"SavePoints {args.save_points}")
        povme_args.append(f"DistanceCutoff {args.distance_cutoff}")
        povme_args.append(f"ConvexHullExclusion {args.convex_hull_exclusion}")
        povme_args.append(f"ContiguousPointsCriteria {args.contiguous_points_criteria}")

        # Pocket seed regions
        if args.pocket_seed_sphere:
            for sphere in args.pocket_seed_sphere:
                povme_args.append(
                    f"ContiguousPocketSeedSphere {' '.join(map(str, sphere))}"
                )
        if args.pocket_seed_box:
            for box in args.pocket_seed_box:
                povme_args.append(f"ContiguousPocketSeedBox {' '.join(map(str, box))}")

        # Additional parameters
        povme_args.append(f"NumProcessors {args.num_processors}")
        povme_args.append(f"UseDiskNotMemory {args.use_disk_not_memory}")
        povme_args.append(f"OutputFilenamePrefix {args.output_prefix}")
        povme_args.append(f"SaveIndividualPocketVolumes {args.save_individual_volumes}")
        povme_args.append(f"SavePocketVolumesTrajectory {args.save_trajectory}")
        povme_args.append(f"OutputEqualNumPointsPerFrame {args.equal_num_points}")
        povme_args.append(f"SaveTabbedVolumeFile {args.save_tabbed}")
        povme_args.append(f"SaveVolumetricDensityMap {args.save_density_map}")
        povme_args.append(f"CompressOutput {args.compress_output}")

        # Run POVME
        RunPOVME(povme_args)

    elif args.command == "pocket_id":
        # Build the argument list for Pocket ID
        pocket_id_args = ["pocket_id"]

        pocket_id_args.append(f"PDBFileName {args.input_pdb}")
        pocket_id_args.append(f"OutputFilenamePrefix {args.output_prefix}")

        # Add more Pocket ID specific arguments here

        # Run Pocket ID
        run_pocket_id(pocket_id_args)

    else:
        parser.print_help()
        sys.exit(1)
