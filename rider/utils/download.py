import os
import urllib.request
import tarfile

USER_DATA_DIR = os.path.expanduser("~/.rider_data")

def ensure_data_exists(download_db=True, download_submodule=True):
    os.makedirs(USER_DATA_DIR, exist_ok=True)

    # download Rider PDB database
    if download_db:
        db_final_dir = os.path.join(USER_DATA_DIR, "Rider_pdb_database", "database")
        if os.path.exists(db_final_dir):
            print(f"✅ Rider database already exists: {db_final_dir}")
        else:
            print("📥 Downloading Rider database...")

            url = "https://zenodo.org/record/XXXXXX/files/database.tar.gz"  
            tar_path = os.path.join(USER_DATA_DIR, "database.tar.gz")

            urllib.request.urlretrieve(url, tar_path)

            print("📦 Extracting database...")
            with tarfile.open(tar_path, "r:gz") as tar:
                tar.extractall(path=USER_DATA_DIR)

            os.remove(tar_path)
            print(f"✅ Database extracted to {db_final_dir}")

    # download submodule
    if download_submodule:
        submodule_dir = os.path.join(USER_DATA_DIR, "submodule")
        if os.path.exists(submodule_dir):
            print(f"✅ Submodule already exists: {submodule_dir}")
        else:
            print("📥 Downloading submodule...")

            url = "https://zenodo.org/record/YYYYYY/files/submodule.tar.gz"  
            tar_path = os.path.join(USER_DATA_DIR, "submodule.tar.gz")

            urllib.request.urlretrieve(url, tar_path)

            print("📦 Extracting submodule...")
            with tarfile.open(tar_path, "r:gz") as tar:
                tar.extractall(path=USER_DATA_DIR)

            os.remove(tar_path)

            # 设置 Foldseek 可执行权限（如果存在）
            foldseek_exec = os.path.join(submodule_dir, "foldseek")
            if os.path.exists(foldseek_exec):
                os.chmod(foldseek_exec, 0o755)

            print(f"✅ Submodule extracted to {submodule_dir}")