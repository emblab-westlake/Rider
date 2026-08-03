import argparse

from rider.tools.download_cli import configure_parser as configure_download_parser
from rider.tools.download_cli import run as run_download

def main():
    parser = argparse.ArgumentParser(prog="rider", description="Rider CLI tool")
    subparsers = parser.add_subparsers(dest="command", required=True, help="Sub-commands")

    parser_download = subparsers.add_parser(
        "download-db",
        help="Download versioned databases and dependencies from Zenodo",
    )
    configure_download_parser(parser_download)
    parser_download.set_defaults(func=run_download)

    args = parser.parse_args()
    try:
        return args.func(args)
    except (OSError, RuntimeError, ValueError) as error:
        parser.error(str(error))
