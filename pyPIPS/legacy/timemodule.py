from datetime import timedelta
import pytz

# Time zone stuff
Central = pytz.timezone('US/Central')

timefmt = '%Y-%m-%d %H:%M:%S %Z%z'
timefmt2 = '%Y-%m-%d %H:%M:%S'
timefmt3 = '%Y%m%d%H%M%S'
timefmt4 = '%m%d'
timefmt5 = '%Y,%m,%d'

timedelta5s = timedelta(seconds=5)
timedelta10s = timedelta(seconds=10)
timedelta30s = timedelta(seconds=30)
timedelta1min = timedelta(minutes=1)
