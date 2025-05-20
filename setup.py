from setuptools import setup, find_packages
import os


def install_dependencies():

	req_path = os.path.join(os.path.dirname(__file__), "requirements.txt")

	if (not os.path.exists(req_path)):

		return[]

	with open(req_path, "r") as file:

		return [line.strip() for line in file if line.strip() and not line.startswith("#")]


setup(
	name="3s2b",
	version="1.0.0",
	packages=find_packages(where="src"),
	package_dir={"": "src"},
	include_package_data=True,
	install_requires=install_dependencies(),
	entry_points={
		"console_scripts": [
			"3s2b=cli:cli",
			"make-db=db.make:makeDB"
		]
	}
)