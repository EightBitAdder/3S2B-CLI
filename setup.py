from setuptools import setup, find_packages
import os
from pathlib import Path


def install_dependencies():

	req_path = Path(__file__).parent / "requirements.txt"

	if (not req_path.exists()):

		return[]

	with req_path.open("r") as file:

		return [
			line.strip()
			for line in file
			if line.strip() and not line.startswith("#")
		]


setup(
	name="3s2b",
	version="1.0.0",
	packages=find_packages(where="src"),
	package_dir={"": "src"},
	include_package_data=True,
	install_requires=install_dependencies(),
	entry_points={
		"console_scripts": [
			"3s2b=main.cli:cli",
			"make-db=db.make:makeDB"
		]
	}
)