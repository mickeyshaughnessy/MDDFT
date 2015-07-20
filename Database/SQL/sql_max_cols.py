#!/usr/bin/env python
import itertools
import sqlite3

db = sqlite3.connect(':memory:')
try:
    for num_columns in itertools.count(1):
        db.execute('CREATE TABLE T%d (%s)' % (num_columns, ','.join('C%d' % i for i in range(num_columns))))
except sqlite3.DatabaseError as ex:
    if 'too many columns' in str(ex):
        print('Max columns = %d' % (num_columns - 1))
