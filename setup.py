from setuptools import setup, find_packages

setup(
    name="closeread",
    version="0.1.0",
    author="Yixin Zhu",
    author_email="yz3398@cornell.edu",
    description="CloseRead: A pipeline for evaluating IG assemblies.",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    url="https://github.com/phylo-lab-usc/CloseRead_ig-assembly-eval",
    packages=find_packages(),
    entry_points={
        "console_scripts": [
            "closeread-pipeline=closeread.pipeline:run_pipeline_cli",
            "closeread-plot=closeread.finalize_plot:run_plot_cli",
        ],
    },
    include_package_data=True,
    package_data={
        "closeread": ["scripts/*.sh"],
    },   
    python_requires=">=3.10",  # Enforce Python 3.10
    classifiers=[
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.10",
        "Operating System :: OS Independent",
    ],
)
