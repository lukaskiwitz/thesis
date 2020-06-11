import setuptools

with open("./README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="thesisproject",  # Replace with your own username
    version="0.0.1",
    author="Lukas Kiwitz",
    author_email="kiwitz@hu-berlin.de",
    description="",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/lukaskiwitz/thesis",
    packages=["thesis/main", "thesis/scenarios", "thesis/example_scripts"],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: Posix",
        "Development Status :: 3 - Alpha",
        "Environment :: Console"
    ],
    python_requires='>=3.6',
)
