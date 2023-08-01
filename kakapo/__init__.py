"""Kakapo."""

from datetime import datetime

date_time = datetime.now()
y = str(date_time.year)

__script_name__ = 'kakapo'
__version__ = '0.9.5'
__description__ = ('Extract and annotate protein family members from '
                   'transcriptomes.')
__author__ = 'Karolis Ramanauskas'
__author_email__ = 'kraman2@uic.edu'
__copyright__ = 'Copyright \u00A9 ' + __author__ + ', ' + y
__license__ = 'GNU General Public License Version 3'
__url__ = 'https://github.com/karolisr/kakapo'

__all__ = ['__author__', '__author_email__', '__description__', '__script_name__', '__url__', '__version__']
