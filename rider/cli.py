import argparse
from rider.predict_pipeline import main as predict_main
from rider.tools.download_cli import main as download_main

def main():
    parser = argparse.ArgumentParser(prog="rider", description="Rider CLI tool")
    subparsers = parser.add_subparsers(dest="command", help="Sub-commands")

    # predict 子命令
    parser_predict = subparsers.add_parser("predict", help="Run prediction pipeline")
    parser_predict.set_defaults(func=predict_main)

    # download-db 子命令
    parser_download = subparsers.add_parser("download-db", help="Download the DB")
    parser_download.set_defaults(func=download_main)

    args = parser.parse_args()

    if hasattr(args, 'func'):
        args.func()
    else:
        parser.print_help()