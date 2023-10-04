"""Provides global variables."""

import os

from kakapo import __script_name__, __version__
from kakapo.utils.misc import python_version
from kakapo.utils.os_diffs import check_os, cpu_count, sys_ram

DEBUG_MODE = False
DEBUG_PROCESSES = False

# System information ---------------------------------------------------------
_ = check_os()
MACHINE_TYPE = _['machine_type']
OS_ID = _['os_id']
OS_STR = _['os_str']
RELEASE_ID = _['release_id']
RELEASE_NAME = _['release_name']
DIST_ID = _['dist_id']
DEBIAN_DISTS = _['debian_dists']
REDHAT_DISTS = _['redhat_dists']
SUPPORTED_DISTS = _['supported_dists']
OS_INFO = OS_STR
ENV: dict[str, str] = dict(os.environ)

if OS_ID == 'mac':
    OS_INFO += ' ' + RELEASE_ID
    if RELEASE_NAME != '':
        OS_INFO += ' (' + RELEASE_NAME + ')'

_, PY_V_STR = python_version()

_ = cpu_count()
NCPU = 2 if _[0] is None else _[0]
NCPUL = 2 if _[1] is None else _[1]
RAM = sys_ram()

os.environ['KKP_OS_ID'] = OS_ID
os.environ['KKP_MACHINE_TYPE'] = MACHINE_TYPE

# Other ----------------------------------------------------------------------
PICKLE_PROTOCOL = 2

# Color support --------------------------------------------------------------
CONSRED = '\033[0;91m'
CONGREE = '\033[0;92m'
CONYELL = '\033[0;93m'
CONBLUE = '\033[0;94m'
CONSDFL = '\033[0m'

# Filesystem paths for kakapo data directories -------------------------------
DIR_USR = os.path.expanduser('~')
DIR_DAT = os.path.join(DIR_USR, '.local', 'share', __script_name__)  # XDG_DATA_HOME
DIR_DEP = os.path.join(DIR_DAT, 'dependencies')
DIR_TAX = os.path.join(DIR_DAT, 'ncbi-taxonomy')
DIR_KRK = os.path.join(DIR_DAT, 'kraken2_dbs')

# Logo -----------------------------------------------------------------------

D = '\033[2;92m'  # Dark Green
G = '\033[0;92m'  # Green
T = '\033[0;91m'  # Red
X = '\033[2;91m'  # Dark Red
W = '\033[0;91m'  # Light Red
Y = '\033[0;93m'  # Yellow
B = '\033[2;90m'  # Black
E = '\033[0m'     # Reset

if ((('COLORTERM' in ENV) and (ENV['COLORTERM'] in ('truecolor', 'xterm-truecolor')))
        or ('TERM' in ENV and ENV['TERM'] in ('alacritty', ))):

    D = f'\033[38;2;{69};{121};{40}m'    # Dark Green
    G = f'\033[38;2;{127};{183};{70}m'   # Green
    T = f'\033[38;2;{254};{116};{85}m'   # Red
    X = f'\033[38;2;{240};{110};{85}m'   # Dark Red
    W = f'\033[38;2;{254};{129};{111}m'  # Light Red
    Y = f'\033[38;2;{255};{183};{82}m'   # Yellow
    B = f'\033[38;2;{0};{0};{0}m'        # Black
    E = '\033[m'                         # Reset

LOGO = f'''              {D}██████████
            {D}██{G}██████████{D}████
          {D}██{G}████████████████{D}██
        {D}██{G}████████████████████{D}██
      {D}██{G}██████{B}██{G}████████████{B}██{D}██
    {D}██{G}████████{B}██{G}██{W}████████{G}██{B}██{G}██{D}██
    {D}██{G}████████████{W}█{B}█{T}██████{B}█{X}█{G}████{D}██
  {D}████{G}████████████{W}██{T}██████{X}██{G}██████{D}██
  {D}██{G}██████████████{W}██{T}██████{X}██{G}████████{D}██
{D}██{G}████████████████{W}██{T}████{X}████{G}████████{D}██
{D}██{G}██████████████████{W}██{T}██{X}██{G}██████████{D}██
{D}██{G}██████████████████{W}██{T}██{X}██{G}██████████{D}██
  {D}██{G}██████████████████{X}██{G}████████████{D}██
    {D}██{G}████████████████████████████████{D}██
    {D}██{G}████████████{Y}████{G}████████████{Y}████{G}█{D}██
    {D}████{G}███████████{Y}██{G}██████████████{Y}██{G}████{D}████
    {D}██{G}██{D}██████{G}█████{Y}██{G}█{Y}████{G}█{Y}████{G}████{Y}██{G}█{Y}████{G}█{Y}████{D}█ {Y}██████   ████
    {D}██{G}████████{D}██{G}███{Y}██{G}██{Y}██{G}█████{Y}██{G}███{Y}██{G}██{Y}██{G}█████{Y}██{D}██{Y}██  ██ ██  ██
    {D}███████████████{Y}█████{D}███{Y}█████{D}███{Y}█████{D}███{Y}█████{D}██{Y}██{D}█ {Y}██ ██  ██
                   {Y}██  ██ ███ ██   ██  ██ ███ ██  ██  ██ ██  ██
                  {Y}████ ███ ██████ ████ ███ ██████ █████   ████
                                                  ██{E}'''

print()
print(LOGO)

# Script Info ----------------------------------------------------------------
SCRIPT_INFO = (''
               + f'{__script_name__.title()} version: {__version__}\n'
               + f'Python version: {PY_V_STR}\n'
               + f'Operating system: {OS_INFO}\n'
               + f'System info: {NCPU} physical ({NCPUL} logical) cores, '
                 f'{RAM:.2f} GB RAM ({MACHINE_TYPE})\n')
