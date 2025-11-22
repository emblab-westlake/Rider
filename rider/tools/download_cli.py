# rider/tools/download_cli.py

from rider.utils.download import ensure_data_exists

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--db", action="store_true")
    parser.add_argument("--submodule", action="store_true")
    args = parser.parse_args()

    ensure_data_exists(download_db=args.db, download_submodule=args.submodule)