#!/usr/bin/env python
####################################################################################
#
# Melbourne Genomics Pipeline Annotation Script
#
# Copyright Melbourne Genomics Health Alliance members. All rights reserved.
#
# DISTRIBUTION:
#
# This source code should not be distributed to a third party without prior
# approval of the Melbourne Genomics Health Alliance steering committee (via
# Natalie Thorne - natalie.thorne@melbournegenomics.org.au).
#
####################################################################################
#
# Purpose:
# * add a pipeline run ID to the sample metadata file, read from a pipeline run ID file
# * takes incoming meta from stdin and writes to stdout
####################################################################################

import argparse
import os.path
import random
import sys

def generate_new_id( f ):
  '''
    given a file, reads the current ID, appends to it, and writes it back to the same file.
    if the file doesn't exist, a random ID is generated.
    format of the ID is site_000000000
  '''
  current_id = get_current_id( f )
  site, run = current_id.rsplit("_", 1)
  # increment run ID
  run = int(run) + 1

  new_id = '%s_%09i' % (site, run)

  fh = open( f, 'w' )
  fh.write( new_id )
  return new_id

def get_current_id( f ):
  '''
    given a file, reads the current ID. if the file doesn't exist, a random ID is generated.
    format of the ID is site_000000000
    @param: filename containing pipeline ID
    @returns: current ID
  '''
  # NOTE! we don't do any file locking. 
  # parallel pipelines could potentially attempt to update the ID simulatenously, resulting in a non-unique ID
  if os.path.isfile( f ):
    fh = open( f, 'r' )
    current = fh.readline().strip()
    fh.close()
    if "_" in current:
      site, run = current.rsplit("_", 1)
      try:
        run = int(run)
      except ValueError: # rightmost section isn't an id after all
        site = current_id
        run = 1
      current_id = '%s_%09i' % (site, run)
    elif len(current) == 0: # empty file
      run = 0
      current_id = 'site%i_%09i' % (random.randint(0, 1e6), run )
    else: # no _
      current_id = '%s_%09i' % (current, 0)
  else: # no file
    run = 0
    current_id = 'site%i_%09i' % (random.randint(0, 1e6), run )

  return current_id

def write( src, target, new_id ):
  '''
    reads lines from src and writes to target, appending the pipeline ID at the start as a new column of a tab separated file
  '''
  first = True
  for line in src:
    if first:
      target.write( 'Pipeline_Run_ID\t%s' % line )
      first = False
    else:
      target.write( '%s\t%s' % ( new_id, line ) )

if __name__ == "__main__":
  parser = argparse.ArgumentParser(description='Generate sample metadata file with pipeline ID')
  parser.add_argument('--id', required=True, help='ID file to read/write')
  parser.add_argument('--increment', type=bool, required=False, default=False, help='Increment the pipeline ID')
  parser.add_argument('--parse', type=bool, required=False, default=False, help='Parse metadata file')
  args = parser.parse_args() 
  if args.increment:
    new_id = generate_new_id( args.id )
  else:
    new_id = get_current_id( args.id )
  if args.parse:
    write( sys.stdin, sys.stdout, new_id )
  else:
    print new_id

