import pandas as pd
import os
import sys

# usage: python reupload_fastqs.py <encode file submission filename>

ifile = sys.argv[1]

df = pd.read_csv(ifile, sep='\t')

import encode_utils as eu
from encode_utils.connection import Connection
conn = Connection("www.encodeproject.org")
for ind, entry in df.iterrows():
    file_id = entry.aliases
    file_path=entry.submitted_file_name
    conn.upload_file(file_id=file_id,
                     file_path=file_path)
