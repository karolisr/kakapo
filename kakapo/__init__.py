# -*- coding: utf-8 -*-

from datetime import datetime

date_time = datetime.now()
y = str(date_time.year)

__script_name__ = 'kakapo'
__version__ = '0.2.2'
__description__ = ('Extract and annotate protein family members from '
                   'transcriptomes.')
__author__ = 'Karolis Ramanauskas'
__author_email__ = 'kraman2@uic.edu'
__copyright__ = 'Copyright \u00A9 ' + __author__ + ', ' + y
__license__ = 'cc-by-sa-4.0'
__url__ = 'https://github.com/karolisr/kakapo'

__all__ = []
