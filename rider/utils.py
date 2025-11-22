import urllib.request
import tarfile

def ensure_data_exists(download_db=True, download_submodule=True):
    os.makedirs(USER_DATA_DIR, exist_ok=True)

    if download_db:
        db_path = os.path.join(USER_DATA_DIR, "Rider_pdb_database")
        if not os.path.exists(db_path):
            print("📥 Downloading Foldseek database...")
            url = "https://zenodo.org/record/XXXXXX/files/database.tar.gz"  
            local_file = "/tmp/foldseek_db.tar.gz"
            urllib.request.urlretrieve(url, local_file)
            with tarfile.open(local_file, "r:gz") as tar:
                tar.extractall(path=db_path)
            print(f"✅ Database extracted to {db_path}")

    if download_submodule:
        submodule_path = os.path.join(USER_DATA_DIR, "submodule")
        if not os.path.exists(submodule_path):
            print("📥 Downloading submodule files...")
            os.makedirs(submodule_path, exist_ok=True)
            # 示例：下载 foldseek 可执行文件
            urllib.request.urlretrieve("https://example.com/foldseek", os.path.join(submodule_path, "foldseek"))
            os.chmod(os.path.join(submodule_path, "foldseek"), 0o755)
            # 你也可以下载模型权重压缩包并解压
            print(f"✅ Submodule downloaded to {submodule_path}")