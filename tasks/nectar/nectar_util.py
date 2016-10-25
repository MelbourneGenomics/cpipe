import os
import hashlib
from os import path
import json

import shutil
from swiftclient.service import SwiftService
import tempfile

from tasks.common import unzip_todir, has_swift_auth, ROOT

current_dir = path.dirname(__file__)
current_manifest = path.join(current_dir, 'current.manifest.json')
target_manifest = path.join(current_dir, 'target.manifest.json')

def add_to_manifest(key):
    with open(target_manifest, 'r') as target, \
            open(current_manifest, 'r+') as current:
        target_json = json.load(target)
        current_json = json.load(current)
        current_json[key] = target_json[key]
        current.seek(0)
        current.truncate()
        json.dump(current_json, current, indent=4)


def create_current_manifest():
    # Create the current manifest if it doesn't exist
    if not path.exists(current_manifest):
        with open(current_manifest, 'w') as current:
            json.dump({}, current)

def nectar_asset_needs_update(asset_key):
    create_current_manifest()
    with open(target_manifest, 'r') as target, \
            open(current_manifest, 'r') as current:

        # Open the input and output json files
        target_json = json.load(target)
        current_json = json.load(current)

        # We need to update if the file doesn't exist or is out of date (wrong hash)
        if asset_key not in current_json or target_json[asset_key]['hash'] != current_json[asset_key]['hash']:
            return True
    return False

def download_nectar_asset(asset_key, to_temp=True):
    create_current_manifest()
    with SwiftService() as swift, \
            open(target_manifest, 'r') as target, \
            open(current_manifest, 'r+') as current:

        # Open the input and output json files
        current.seek(0)
        target_json = json.load(target)

        # If we're downloading to a temp directory, do so. Otherwise use the path listed in the manifest,
        # clearing the directory if it exists
        if to_temp:
            download_dir = tempfile.mkdtemp()
        else:
            download_dir = os.path.join(ROOT, target_json[asset_key]['path'])
            if os.path.isdir(download_dir):
                shutil.rmtree(download_dir)
            os.makedirs(download_dir)

        for result in swift.download(
                container='cpipe-2.4-assets',
                objects=[os.path.join(target_json[asset_key]['path'], target_json[asset_key]['version'] + '.tar.gz')],
                options={'out_directory': download_dir}
        ):
            zip_file = result['path']
            target_hash = target_json[asset_key]['hash']

            if not result['success']:
                print('\t' + asset_key + '... FAILED! ' + str(result['error']))
                raise IOError(result['error'])

            # sha1hash the zip file to ensure its integrity
            with open(zip_file, 'r') as zip_handle:
                current_hash = hashlib.sha1(zip_handle.read()).hexdigest()
                if current_hash != target_hash:
                    raise "{0} failed hashsum check! Check its integrity or update and commit your target.manifest.json".format(
                        asset_key)

                # Unzip into the temporary folder removing the outer directory
                zip_handle.seek(0)
                unzip_todir(zip_handle, download_dir, 'tgz')

                # And delete the zip file
                os.remove(zip_file)

                # Log success
                print('\t' + asset_key + '... done.')

                # Return the download location so we can pass it to the install tasks
                return download_dir

def nectar_task(asset_key, to_temp=True):
    def action():
        # Download the asset from nectar
        dir = download_nectar_asset(asset_key, to_temp)
        # If it's not being downloaded to a temporary directory, its already installed, so update the manifest
        if not to_temp:
            add_to_manifest(asset_key)

        # If it is being downloaded to a temporary directory, it's going to be installed, so save the download directory for the installer
        else:
            return {'dir': dir}
    return {
        # Download the asset and return the directory as a doit arg
        'actions': [action],
        'uptodate': [not nectar_asset_needs_update(asset_key)]
    }
