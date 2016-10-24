import os
import hashlib
from os import path
import json
from swiftclient.service import SwiftService
import tempfile

from tasks.common import unzip_todir, has_swift_auth

current_dir = path.dirname(__file__)
current_manifest = path.join(current_dir, 'current.manifest.json')
target_manifest = path.join(current_dir, 'target.manifest.json')


def nectar_asset_needs_update(asset_key):
    with open(target_manifest, 'r') as target, \
            open(current_manifest, 'r') as current:

        # Open the input and output json files
        target_json = json.load(target)
        current_json = json.load(current)

        # We need to update if the file doesn't exist or is out of date (wrong hash)
        if asset_key not in current_json or target_json[asset_key]['hash'] != current_json[asset_key]['hash']:
            return True
    return False

def download_nectar_asset(asset_key):
    with SwiftService() as swift, \
            open(target_manifest, 'r') as target, \
            open(current_manifest, 'r+') as current:

        # Open the input and output json files
        current.seek(0)
        target_json = json.load(target)
        try:
            current_json = json.load(current)
        except:
            current_json = {}

        # Do the download and update the list of downloaded assets
        download_dir = tempfile.mkdtemp()
        for result in swift.download(
                container='cpipe-2.4-assets',
                objects=target_manifest[asset_key],
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

                # Update the list of currently installed assets
                current_json[asset_key] = target_json[asset_key]

                # Log sucess
                print('\t' + asset_key + '... done.')

                # Delete the temp dir
                # shutil.rmtree(download_dir)

                # Write out the updated list
                current.seek(0)
                current.write(json.dumps(current_json, indent=4))

                # Return the download location so we can pass it to the install tasks
                return download_dir

def nectar_task(asset_key):
    return {
        # Download the asset and return the directory as a doit arg
        'actions': [lambda: {asset_key + '_dir': download_nectar_asset(asset_key)}],
        'uptodate': [nectar_asset_needs_update(asset_key)]
    }
