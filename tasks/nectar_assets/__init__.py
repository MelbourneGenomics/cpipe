#!/usr/bin/python
"""
Defines tasks that are involved in downloading cpipe assets from NECTAR
"""

# Imports
import json
import os
import hashlib
from os import path
from subprocess import check_call
from swiftclient.service import SwiftService
from tasks.nectar_assets.dependencies import *

current_dir = path.dirname(__file__)
root = path.realpath(path.join(current_dir, '..', '..'))
temp = path.join(current_dir, 'temp')
current_manifest = path.join(current_dir, 'current.manifest.json')


def task_setup_manifests():
    return {
        'actions': [create_current_manifest],
        'uptodate': [True]     
    }


def create_current_manifest():
    # Create the current manifest if it doesn't exist
    if not path.exists(current_manifest):
        with open(current_manifest, 'w') as current:
            json.dump({}, current)

def assets_needing_update():
    create_current_manifest()

    with open(path.join(current_dir, 'target.manifest.json'), 'r') as target, \
            open(current_manifest, 'r') as current:

        # Work out which files we need to download
        to_download = []

        # Open the input and output json files
        target_json = json.load(target)
        current_json = json.load(current)

        # Check each key of the target manifest
        for key in target_json:

            # We need to update if the file doesn't exist or is out of date (wrong hash)
            if key not in current_json or target_json[key]['hash'] != current_json[key]['hash']:
                to_download.append('{path}/{version}.tar.gz'.format(**target_json[key]))

        return to_download


def download_nectar_assets():
    print('Updating Cpipe assets...')

    to_download = assets_needing_update()

    with SwiftService() as swift, \
            open(path.join(current_dir, 'target.manifest.json'), 'r') as target, \
            open(current_manifest, 'r+') as current:

        # Open the input and output json files
        current.seek(0)
        target_json = json.load(target)
        try:
            current_json = json.load(current)
        except:
            current_json = {}

        # Maps paths back to keys
        reverse_lookup = {target_json[name]['path']: name for name in target_json}

        # Do the download and update the list of downloaded assets
        for result in swift.download(
                container='cpipe-2.4-assets',
                objects=to_download,
                options={'out_directory': str(root)}
        ):
            zip_file = result['path']
            target_dir = path.dirname(zip_file)
            asset_key = reverse_lookup[path.dirname(result['object'])]
            target_hash = target_json[asset_key]['hash']

            if not result['success']:
                print('\t' + asset_key + '... FAILED! ' + str(result['error']))
                continue

            # sha1hash the zip file to ensure its integrity
            with open(zip_file, 'r') as zip_handle:
                current_hash = hashlib.sha1(zip_handle.read()).hexdigest()
                if current_hash != target_hash:
                    raise "{0} failed hashsum check! Check its integrity or update and commit your target.manifest.json".format(
                        asset_key)

            # Unzip, removing the outer directory
            check_call(['tar', '-xzf', zip_file, '-C', str(target_dir), '--strip-components=1'])

            # And delete the zip file
            os.remove(zip_file)

            # Update the list of currently installed assets
            current_json[asset_key] = target_json[asset_key]

            print('\t' + asset_key + '... done.')

        # Write out the updated list
        current.seek(0)
        current.write(json.dumps(current_json, indent=4))


def task_nectar_assets():
    return {
        'setup': ['setup_manifests'],
        'uptodate': [lambda: len(assets_needing_update()) == 0],
        'actions': [download_nectar_assets]
    }
