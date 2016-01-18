#!/usr/bin/env python
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

import argparse
import datetime
import glob
import os
import os.path
import sys

DEFAULT_PRIORITY = '1'

def write_log( log, msg ):
  log.write( '%s: %s\n' % ( datetime.datetime.now().strftime( '%y%m%d-%H%M%S' ), msg ) )

def add_profile( profile_name, fh_out ):
  '''
    add ./designs/profile_name/profile_name.genes.txt
  '''
  directory = './designs/{0}'.format( profile_name )
  target = './designs/{0}/{0}.genes.txt'.format( profile_name )
  if os.path.isfile(target):
    fh_out.write( 'ERROR: profile "{0}" already exists\n'.format( profile_name ) )
  else:
    if not os.path.exists( directory ):
      os.makedirs( directory )
      with open( target, 'w' ) as fh:
        fh.write( '# profile created {0}'.format( datetime.datetime.now().strftime( '%y%m%d' ) ) )

def list_profiles( fh_out ):
  '''
    look for files of the form designs/X/X.genes.txt
  '''
  for profile in glob.glob( './designs/*/*.genes.txt' ):
    if 'genelists' not in profile:
      fh_out.write( '{0}\n'.format( os.path.basename(profile).split('.')[0] ) )

def list_genes( profile, fh_out ):
  for line in open( './designs/{0}/{0}.genes.txt'.format( profile ), 'r' ):
    if line.startswith( '#' ):
      continue
    fh_out.write( '{0}\n'.format( line.strip().split()[0] ) )

def add_genes( profile, fh_in, fh_out, force=False ):
  # load existing genes from the list
  existing = {}
  for line in open( './designs/{0}/{0}.genes.txt'.format( profile ), 'r' ):
    if line.startswith( '#' ):
      continue
    gene, priority = [ x.strip() for x in line.strip().split() ]
    existing[ gene.upper() ] = priority

  # build incidentalome and exon list
  validation = build_validation_sets()

  new_genes = {}
  requested = 0
  proceed = True
  for line in fh_in:
    gene = line.strip().upper()
    requested += 1
    if gene not in existing:
      new_genes[ gene ] = DEFAULT_PRIORITY
    # validate
    if gene in validation['incidentalome']:
      fh_out.write( 'WARNING: {0} is on the incidentalome\n'.format( gene ) )
      proceed = False
    if gene not in validation['exons']:
      fh_out.write( 'WARNING: {0} not found in the all genes list\n'.format( gene ) )
      proceed = False
      
  if proceed or force:
    if len( new_genes ) > 0:
      fh_out.write( 'Adding {0} new gene(s) out of {1} requested...\n'.format( len( new_genes ), requested ) )
      existing.update( new_genes )
      with open( './designs/{0}/{0}.genes.txt'.format( profile ), 'w' ) as target:
        target.write( '# version {0}\n'.format( datetime.datetime.now().strftime( '%y%m%d' ) ) )
        target.write( '# added {0}\n'.format( ' '.join( sorted( list( new_genes ) )  ) ) )
        for gene in sorted( existing.keys() ):
          target.write( '{0}\t{1}\n'.format( gene, existing[gene] ) )
      fh_out.write( 'Done\n' )
    else:
      fh_out.write( 'No new genes to add\n' )
        
  else: # !proceed and !force
    fh_out.write( 'Not adding genes due to previous warnings. Use --force if you are sure you want to add these genes\n' )

def remove_genes( profile, fh_in, fh_out, force=False ):
  # load existing genes from the list
  existing = {}
  for line in open( './designs/{0}/{0}.genes.txt'.format( profile ), 'r' ):
    if line.startswith( '#' ):
      continue
    gene, priority = [ x.strip() for x in line.strip().split() ]
    existing[ gene.upper() ] = priority

  deleted = set()
  proceed = True
  for line in fh_in:
    gene = line.strip().upper()
    if gene in existing:
      deleted.add( gene )
      del existing[gene]
    else:
      proceed = False
      fh_out.write( 'WARNING: {0} not found in existing profile\n'.format( gene ) )

  if proceed or force:
    fh_out.write( 'Removing {0} genes...\n'.format( len( deleted ) ) )
    with open( './designs/{0}/{0}.genes.txt'.format( profile ), 'w' ) as target:
      target.write( '# version {0}\n'.format( datetime.datetime.now().strftime( '%y%m%d' ) ) )
      target.write( '# removed {0}\n'.format( ' '.join( sorted( list( deleted ) ) ) ) )
      for gene in sorted( existing.keys() ):
        target.write( '{0}\t{1}\n'.format( gene, existing[gene] ) )
    fh_out.write( 'Done\n' )
  else:
    fh_out.write( 'Not updating gene list due to previous warnings. Use --force if you are sure you want to add these genes\n' )
 
def build_validation_sets():
  incidentalome = set()
  for line in open( './designs/genelists/incidentalome.genes.txt', 'r' ):
    if line.startswith( '#' ):
      continue
    incidentalome.add( line.strip().split()[0].upper() )

  exons = set()
  for line in open( './designs/genelists/exons.bed', 'r' ):
    if line.startswith( '#' ):
      continue
    exons.add( line.strip().split()[3].upper() )

  return { 'incidentalome': incidentalome, 'exons': exons }

def validate( profile, fh_out ):
  '''
    check against incidentalome and exons
  '''
  validation = build_validation_sets()

  found_incidentalome = set()
  not_found_exons = set()
  all = set()
  for line in open( './designs/{0}/{0}.genes.txt'.format( profile ), 'r' ):
    gene = line.strip().split()[0].upper()
    all.add( gene )
    if gene in validation['incidentalome']:
      found_incidentalome.add( gene )
    if gene not in validation['exons']:
      not_found_exons.add( gene )
  fh_out.write( 'examined:              \t{0}\n'.format( len( all ) ) )
  fh_out.write( 'found on incidentalome:\t{0}\t{1}\n'.format( len( found_incidentalome ), ' '.join( sorted( list( found_incidentalome ) ) ) ) )
  fh_out.write( 'not found:             \t{0}\t{1}\n'.format( len( not_found_exons ), ' '.join( sorted( list( not_found_exons ) ) ) ) )

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='Manage gene lists')
  parser.add_argument('command', help='command to execute', choices=['add_profile', 'list_profiles', 'list_genes', 'add_genes', 'remove_genes', 'validate'] )
  parser.add_argument('--profile', required=False, help='profile to update')
  parser.add_argument('--force', action='store_true', help='force addition of genes')
  args = parser.parse_args()
  if args.command == 'list_profiles':
    list_profiles( sys.stdout )
  else:
    if not args.profile:
      parser.print_help()
      sys.exit(1)

    if args.command == 'add_profile':
      add_profile( args.profile, sys.stdout )
    elif args.command == 'list_genes':
      list_genes( args.profile, sys.stdout )
    elif args.command == 'add_genes':
      add_genes( args.profile, sys.stdin, sys.stdout, args.force )
    elif args.command == 'remove_genes':
      remove_genes( args.profile, sys.stdin, sys.stdout, args.force )
    elif args.command == 'validate':
      validate( args.profile, sys.stdout )
