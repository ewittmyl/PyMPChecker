#!/bin/csh
#
# Run this from the PyMPChecker directory, as
# ./scripts/MP_Database_Update.csh
#
# If you normally startup using pipenv shell, running this way seems to
# bypass the environment. Hence the need to call python3 specifically
# (rather than just python) in the call below. Not sure how to fix this...
# - dkg 2019 Jun

echo I am going to update the MP database ....

# wget http://www.minorplanetcenter.org/iau/MPCORB/MPCORB.DAT
wget --no-check-certificate https://www.minorplanetcenter.org/iau/MPCORB/MPCORB.DAT

python3 ./MPCORB2MonthlyCatalog.py ./MPCORB.DAT

rm MPCORB.DAT

exit
