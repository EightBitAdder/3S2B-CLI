import subprocess
import sys
import os
from pathlib import Path
import platform


VENV_DIR = Path(".venv")


def make_venv():

    if (not VENV_DIR.exists()):

        print("[INFO] Creating Virtual Environment . . .")

        subprocess.check_call([sys.executable, "-m", "venv", str(VENV_DIR)])

    else:

        print("[INFO] Virtual Environment already exists.")


def get_venv_py():

    if (os.name == "nt"):

        return VENV_DIR / "Scripts" / "python.exe"

    else:

        return VENV_DIR / "bin" / "python"


def venv_run(args):

    py = get_venv_py()

    if (not py.exists()):

        print("[ERROR] Virtual Environment broken.")

        sys.exit(1)

    return subprocess.check_call([str(py)] + args)


def e_install():

    print("[INFO] Upgrading pip . . .")

    venv_run(["-m", "pip", "install", "--upgrade", "pip"])

    print("[INFO] Installing project in editable mode . . .")

    venv_run(["-m", "pip", "install", "-e", "."])


def make_db():

    print("[INFO] Executing `make-db`. . .")

    if (os.name == "nt"):

        exe = VENV_DIR / "Scripts" / "make-db.exe"

    else:

        exe = VENV_DIR / "bin" / "make-db"


    if (not exe.exists()):

        print("[ERROR] `make-db` not found.")

        sys.exit(1)

    subprocess.check_call([str(exe)])


def main():

    print("[SETUP] Bootstrapping Project Environment . . .\n")

    make_venv()
    e_install()
    make_db()

    print("\n[SUCCESS] Setup complete!")


if __name__ == "__main__":

    main()