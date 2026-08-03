"""Download versioned Rider resources from Zenodo."""

from __future__ import annotations

import hashlib
import os
from pathlib import Path
import posixpath
import tarfile
import tempfile
import urllib.error
import urllib.request


ZENODO_RECORD = "19247869"
ZENODO_DOI = "10.5281/zenodo.19247869"
ZENODO_CONTENT_ROOT = f"https://zenodo.org/api/records/{ZENODO_RECORD}/files"
CHUNK_SIZE = 8 * 1024 * 1024

RESOURCES = {
    "rdsdb30": {
        "filename": "Rider_RDSDB30.tar.gz",
        "size": 587_257_550,
        "md5": "ce75285dee2a2933fc45b800b3d95c07",
        "description": "30%-identity non-redundant RdRp structure database (recommended)",
        "root": "Rider_RDSDB30",
        "archive_prefix": "root/gaoyang/westlake_emblab/Rider/Rider_RDSDB_30",
        "required": ("Rider_RDSDB30/pdbs",),
    },
    "rdsdb": {
        "filename": "Rider_RDSDB.tar.gz",
        "size": 6_749_359_564,
        "md5": "97a74266f002edd8c080a44744f70415",
        "description": "full RdRp-domain-specific structure database",
        "root": "Rider_RDSDB",
        "archive_prefix": "usr/commondata/public/gaoyang/Rider_pdb_database/Rider_RDSDB",
        "required": ("Rider_RDSDB/pdbs",),
    },
    "rsdb": {
        "filename": "Rider_RSDB.tar.gz",
        "size": 7_975_030_936,
        "md5": "29b010487910c91722832a40fee46c51",
        "description": "broad Rider structure library",
        "root": "Rider_RSDB",
        "archive_prefix": "usr/commondata/public/gaoyang/Rider_pdb_database/Rider_RSDB",
        "required": ("Rider_RSDB/pdbs",),
    },
    "submodule": {
        "filename": "submodule.tar.gz",
        "size": 8_166_034_733,
        "md5": "856f8c607bc2944f51ca9b7190b80a54",
        "description": "ESM2, ESMFold and Foldseek dependencies",
        "root": "submodule",
        "archive_prefix": "submodule",
        "required": (
            "submodule/esm2_t12_35M_UR50D/model.safetensors",
            "submodule/esmfold_v1/pytorch_model.bin",
            "submodule/foldseek/bin/foldseek",
        ),
    },
}


def default_rider_root() -> Path:
    """Locate a source checkout without installing data inside site-packages."""
    current = Path.cwd().resolve()
    package_root = Path(__file__).resolve().parents[2]
    for candidate in (current, package_root):
        if (candidate / "setup.py").is_file() and (candidate / "rider").is_dir():
            return candidate
    return Path.home() / ".rider_data"


def format_size(size: int) -> str:
    value = float(size)
    for unit in ("B", "KiB", "MiB", "GiB"):
        if value < 1024 or unit == "GiB":
            return f"{value:.1f} {unit}"
        value /= 1024
    raise AssertionError("unreachable")


def resource_table() -> str:
    lines = [f"Zenodo DOI: https://doi.org/{ZENODO_DOI}"]
    for name, metadata in RESOURCES.items():
        lines.append(
            f"  {name:9s} {format_size(metadata['size']):>10s}  "
            f"{metadata['description']}"
        )
    return "\n".join(lines)


def _md5(path: Path) -> str:
    digest = hashlib.md5()  # nosec B324: used only for Zenodo file integrity
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(CHUNK_SIZE), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _download(url: str, destination: Path, expected_size: int) -> None:
    partial = destination.with_name(destination.name + ".part")
    offset = partial.stat().st_size if partial.exists() else 0
    headers = {"User-Agent": "Rider/1.0.0"}
    if offset:
        headers["Range"] = f"bytes={offset}-"

    request = urllib.request.Request(url, headers=headers)
    try:
        response = urllib.request.urlopen(request)
    except urllib.error.HTTPError as error:
        if error.code == 416 and offset == expected_size:
            partial.replace(destination)
            return
        raise

    status = getattr(response, "status", response.getcode())
    if offset and status != 206:
        print("Server did not accept the resume request; restarting the download.")
        offset = 0
    mode = "ab" if offset and status == 206 else "wb"

    downloaded = offset
    with response, partial.open(mode) as output:
        while True:
            chunk = response.read(CHUNK_SIZE)
            if not chunk:
                break
            output.write(chunk)
            downloaded += len(chunk)
            percent = 100 * downloaded / expected_size
            print(
                f"\rDownloading {destination.name}: "
                f"{format_size(downloaded)} / {format_size(expected_size)} "
                f"({percent:.1f}%)",
                end="",
                flush=True,
            )
    print()

    if downloaded != expected_size:
        raise RuntimeError(
            f"Incomplete download for {destination.name}: "
            f"expected {expected_size} bytes, received {downloaded}. "
            "Run the command again to resume."
        )
    partial.replace(destination)


def _validate_member(member: tarfile.TarInfo, output_dir: Path) -> None:
    root = output_dir.resolve()
    target = (root / member.name).resolve()
    if os.path.commonpath((root, target)) != str(root):
        raise RuntimeError(f"Unsafe archive path: {member.name}")
    if member.issym() or member.islnk():
        link = (target.parent / member.linkname).resolve()
        if os.path.commonpath((root, link)) != str(root):
            raise RuntimeError(f"Unsafe archive link: {member.name}")


def _normalize_member(member: tarfile.TarInfo, resource: str) -> tarfile.TarInfo:
    metadata = RESOURCES[resource]
    prefix = metadata["archive_prefix"].rstrip("/")
    original = member.name.rstrip("/")
    if original == prefix:
        relative = ""
    elif original.startswith(prefix + "/"):
        relative = original[len(prefix) + 1 :]
    else:
        raise RuntimeError(
            f"Unexpected path in {metadata['filename']}: {member.name}"
        )

    member.name = posixpath.join(metadata["root"], relative)
    if member.islnk():
        link = member.linkname.rstrip("/")
        if link == prefix:
            member.linkname = metadata["root"]
        elif link.startswith(prefix + "/"):
            member.linkname = posixpath.join(
                metadata["root"], link[len(prefix) + 1 :]
            )
    return member


def _resource_is_installed(resource: str, rider_root: Path) -> bool:
    metadata = RESOURCES[resource]
    for relative in metadata["required"]:
        required = rider_root / relative
        if not required.exists():
            return False
        if required.is_dir() and not any(required.iterdir()):
            return False
    return True


def extract_archive(resource: str, archive_path: Path, rider_root: Path) -> None:
    metadata = RESOURCES[resource]
    target_root = rider_root / metadata["root"]
    if target_root.exists():
        if _resource_is_installed(resource, rider_root):
            print(f"Verified installed resource; extraction skipped: {target_root}")
            return
        raise RuntimeError(
            f"Incomplete installation already exists at {target_root}. "
            "Move it aside before extracting the verified archive."
        )

    staging_root = Path(tempfile.mkdtemp(prefix=".rider_extracting_", dir=rider_root))
    print(f"Extracting {archive_path.name} into a temporary directory ...")
    seen_expected_root = False
    with tarfile.open(archive_path, "r:gz") as archive:
        for member in archive:
            member = _normalize_member(member, resource)
            _validate_member(member, staging_root)
            seen_expected_root = True
            archive.extract(member, staging_root)
    if not seen_expected_root or not _resource_is_installed(resource, staging_root):
        raise RuntimeError(f"Extracted {archive_path.name}, but required Rider files are missing")

    foldseek = staging_root / "submodule" / "foldseek" / "bin" / "foldseek"
    if resource == "submodule" and foldseek.is_file():
        foldseek.chmod(foldseek.stat().st_mode | 0o111)
    (staging_root / metadata["root"]).replace(target_root)
    staging_root.rmdir()
    print(f"Installed {resource} at {target_root}")


def download_resource(
    resource: str,
    rider_root: Path | str | None = None,
    install: bool = True,
) -> Path:
    metadata = RESOURCES[resource]
    rider_root = Path(rider_root or default_rider_root()).expanduser().resolve()
    rider_root.mkdir(parents=True, exist_ok=True)
    archive_dir = rider_root / ".rider_downloads"
    archive_dir.mkdir(parents=True, exist_ok=True)
    destination = archive_dir / metadata["filename"]

    if destination.exists():
        if destination.stat().st_size == metadata["size"] and _md5(destination) == metadata["md5"]:
            print(f"Using verified existing file: {destination}")
        else:
            raise RuntimeError(
                f"Existing file does not match the Zenodo record: {destination}. "
                "Move it aside and rerun the command."
            )
    else:
        url = f"{ZENODO_CONTENT_ROOT}/{metadata['filename']}/content"
        _download(url, destination, metadata["size"])
        checksum = _md5(destination)
        if checksum != metadata["md5"]:
            raise RuntimeError(
                f"MD5 mismatch for {destination}: expected {metadata['md5']}, got {checksum}"
            )
        print(f"Verified MD5: {checksum}")

    if install:
        extract_archive(resource, destination, rider_root)
    return destination
