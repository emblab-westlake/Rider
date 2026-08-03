"""Command-line interface for versioned Rider resource downloads."""

import argparse

from rider.utils.download import RESOURCES, default_rider_root, download_resource, resource_table


def configure_parser(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "resources",
        nargs="*",
        metavar="RESOURCE",
        help="one or more resources to download; use --list to inspect sizes",
    )
    parser.add_argument(
        "--rider-root",
        default=None,
        help=(
            "Rider repository root. By default, use the current/source checkout; "
            "otherwise use ~/.rider_data"
        ),
    )
    parser.add_argument(
        "--download-only",
        action="store_true",
        help="verify and retain archives without installing them",
    )
    parser.add_argument(
        "--list",
        action="store_true",
        help="list available Zenodo resources without downloading",
    )


def run(args: argparse.Namespace) -> int:
    if args.list:
        print(resource_table())
        return 0
    if not args.resources:
        raise ValueError("select at least one RESOURCE or use --list")
    unknown = sorted(set(args.resources) - set(RESOURCES))
    if unknown:
        raise ValueError(
            f"unknown resource(s): {', '.join(unknown)}; "
            f"choose from {', '.join(sorted(RESOURCES))}"
        )
    rider_root = args.rider_root or default_rider_root()
    print(f"Rider root: {rider_root}")
    for resource in args.resources:
        download_resource(
            resource,
            rider_root=rider_root,
            install=not args.download_only,
        )
    return 0


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(
        prog="rider-download-db",
        description="Download versioned Rider resources from Zenodo.",
    )
    configure_parser(parser)
    args = parser.parse_args(argv)
    try:
        return run(args)
    except (OSError, RuntimeError, ValueError) as error:
        parser.error(str(error))


if __name__ == "__main__":
    main()
