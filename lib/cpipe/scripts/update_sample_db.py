#!/usr/bin/env python
'''
###########################################################################
#
# This file is part of Cpipe.
#
# Cpipe is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, under version 3 of the License, subject
# to additional terms compatible with the GNU General Public License version 3,
# specified in the LICENSE file that is part of the Cpipe distribution.
#
# Cpipe is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Cpipe.  If not, see <http:#www.gnu.org/licenses/>.
#
###########################################################################
#
# writes run information to a database for each sample
# see main() for usage info
#
###########################################################################
'''

import datetime
import os
import sqlite3
import sys

# use this to manage schema versions
SCHEMA_VERSION=0

def update_sample_db(db, sample, run_id, analysis, capture, pipeline_version, log):
    '''
        create table if required 
    '''
    conn = sqlite3.connect(db)

    # ensure database schema is up to date
    conn.execute('create table if not exists version (version, created)')
    # check version
    version = conn.execute('select max(version) from version').fetchone()[0]
    if version is None: # create new
        log.write('update_sample_db: creating new schema\n')
        conn.execute('create table if not exists analysis (sample, run_id, location, created, analysis_type, capture, pipeline_version)')
        conn.execute('insert into version values (?, ?)', (SCHEMA_VERSION, datetime.datetime.now()))
    elif version < SCHEMA_VERSION: # upgrade
        log.write('update_sample_db: upgrading schema from {0} to {1}\n'.format(version, SCHEMA_VERSION))
        conn.execute('insert into version values (?, ?)', (SCHEMA_VERSION, datetime.datetime.now()))
    else:
        log.write('update_sample_db: schema is up to date at version {0}\n'.format(version))


    # now write record
    log.write('update_sample_db: writing record to {0}...\n'.format(db))
    location = os.getcwd()
    conn.execute('insert into analysis values (?, ?, ?, ?, ?, ?, ?)', (sample, run_id, location, datetime.datetime.now(), analysis, capture, pipeline_version))
    conn.commit()
    log.write('update_sample_db: writing record to {0}: done\n'.format(db))

def main():
    '''
        run from command line
    '''
    import argparse
    parser = argparse.ArgumentParser(description='write details of analysis for this sample to a sqlite database')
    parser.add_argument('--db', required=True, help='database file')
    parser.add_argument('--sample', required=True, help='sample name')
    parser.add_argument('--run_id', required=True, help='pipeline run id')
    parser.add_argument('--analysis', required=True, help='type of analysis')
    parser.add_argument('--capture', required=True, help='exome capture file')
    parser.add_argument('--pipeline_version', required=True, help='version of Cpipe')
    args = parser.parse_args()
    update_sample_db(args.db, args.sample, args.run_id, args.analysis, args.capture, args.pipeline_version, sys.stderr)

if __name__ == '__main__':
    main()
