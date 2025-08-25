# setup.py
from setuptools import setup, find_packages
from pathlib import Path

# readme
HERE = Path(__file__).resolve().parent
README = (HERE / "README.md").read_text(encoding="utf-8") if (HERE / "README.md").exists() else ""

# NOTE:
# - We intentionally avoid pinning GPU / Cuda-specific packages like `torch==2.2.0` here
#   because users should install the appropriate PyTorch build for their GPU/OS
#   via conda/pytorch channel before installing this package.
# - External command-line tools (mmseqs2, foldseek, diamond, etc.) are not bundled;
#   users should install them separately (conda-forge / bioconda recommended).

setup(
    name="rider",
    version="1.0.0",
    description="Rider: RNA virus identification pipeline",
    long_description=README,
    long_description_content_type="text/markdown",
    author="Gaoyang Luo",
    author_email="lgyjsnjhit@gmail.com",
    url="https://github.com/your-org/Rider",  # 修改为你的仓库 URL
    license="MIT",
    packages=find_packages(exclude=("tests", "docs")),
    include_package_data=True,
    python_requires=">=3.10",
    install_requires=[
        # Core Python packages (selected from your requirements.txt)
        "absl-py>=2.1.0",
        "accelerate>=0.26.1",
        "aiohttp>=3.9.5",
        "datasets>=2.19.2",
        "einops>=0.8.0",
        "fair-esm>=2.0.0",
        "fsspec>=2024.2.0",
        "gitpython>=3.1.41",
        "h5py>=3.11.0",
        "huggingface-hub>=0.23.3",
        "matplotlib>=3.8.4",
        "numpy>=1.21",
        "pandas>=2.2.0",
        "pytorch-lightning>=2.4.0",
        "pyyaml>=6.0.1",
        "scikit-learn>=1.4.0",
        "scipy>=1.12.0",
        "tokenizers>=0.15.1",
        "transformers>=4.37.2",
        "umap-learn>=0.5.6",
        "wandb>=0.16.4",
        "biopython>=1.80",
        "safetensors>=0.3.0",
        "psutil>=5.9.0",
    ],
    extras_require={
        # Optional extras: developer / heavy-GPU dependencies can be documented here.
        "dev": ["pytest", "black", "flake8"],
    },
    entry_points={
        "console_scripts": [
            # these two assume modules will be at rider/predict_pipline.py and rider/predict_pipline_light.py
            "rider-predict = rider.predict_pipline:main",
            "rider-predict-light = rider.predict_pipline_light:main",
        ]
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
    ],
)