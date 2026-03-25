import argparse
import platform
import sys
from importlib.metadata import version as get_pip_version

os_name = platform.system()
arch = platform.machine()
python_version = platform.python_version()

def get_version():
    try:
        return get_pip_version("braid")
    except Exception:
        return "Unknown (Dev)"

__version__ = get_version()
BUILD_DATE = "2026-03-25"
GIT_COMMIT = "Latest"

class FancyVersionAction(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        BOLD = "\033[1m"
        DIM = "\033[2m"
        BLUE = "\033[38;5;39m"
        CYAN = "\033[38;5;51m"
        GREEN = "\033[38;5;47m"
        RESET = "\033[0m"

        logo = f"""{BLUE}{BOLD}
  ____            _     _
 | __ ) _ __ __ _(_) __| |
 |  _ \\| '__/ _` | |/ _` |
 | |_) | | | (_| | | (_| |
 |____/|_|  \\__,_|_|\\__,_|{RESET}
"""
        print(logo)
        print(f"{BOLD}Braid {GREEN}v{__version__}{RESET}")
        print(f"{DIM}Block Resolution and Annotation of Integrated DNA{RESET}\n")

        print(f"{BOLD}Build Information:{RESET}")
        print(f"  {CYAN}Revision:{RESET}    {GIT_COMMIT}")
        print(f"  {CYAN}Built:{RESET}       {BUILD_DATE}")
        print(f"  {CYAN}Platform:{RESET}    {os_name} ({arch})")
        print(f"  {CYAN}Python Env:{RESET}  v{python_version}\n")

        sys.exit(0)