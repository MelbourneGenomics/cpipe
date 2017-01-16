import os
import hashlib
from os import path
import json

import shutil
from swiftclient.service import SwiftService
import tempfile
from urllib.parse import urlparse, urljoin
from urllib.request import urlretrieve

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

        # Calculate the download directory and file names
        # Directory to download into
        download_dir = tempfile.mkdtemp()
        # The download file name without any directories
        file_name = target_json[asset_key]['version'] + '.tar.gz'
        # The download file with path relative to the cpipe directory, e.g. tools/gatk/3.6.tar.gz
        file_path = os.path.join(target_json[asset_key]['path'], file_name)
        container = target_json[asset_key]['container']

        # If we're downloading to a temp directory, do so. Otherwise use the path listed in the manifest,
        # clearing the directory if it exists
        if to_temp:
            unzip_dir = tempfile.mkdtemp()
        else:
            unzip_dir = os.path.join(ROOT, target_json[asset_key]['path'])
            if os.path.isdir(unzip_dir):
                shutil.rmtree(unzip_dir)
            os.makedirs(unzip_dir)

        # Do the download using either the swift client for private files, or urlretrieve for public ones
        parsed = urlparse(container)
        if parsed.scheme.startswith('http'):
            full_url = urljoin(container, file_path)
            zip_file = os.path.join(download_dir, file_name)
            urlretrieve(full_url, zip_file)
        else:
            for result in swift.download(
                    container=container,
                    objects=[file_path],
                    options={'out_directory': download_dir}
            ):
                zip_file = result['path']

                if not result['success']:
                    print(('\t' + asset_key + '... FAILED! ' + str(result['error'])))
                    raise IOError(result['error'])

        # sha1hash the zip file to ensure its integrity
        target_hash = target_json[asset_key]['hash']
        with open(zip_file, 'rb') as zip_handle:
            current_hash = hashlib.sha1(zip_handle.read()).hexdigest()
            if current_hash != target_hash:
                raise "{0} failed hashsum check! Check its integrity or update and commit your target.manifest.json".format(
                    asset_key)

            # Rewind the zip handle
            zip_handle.seek(0)

            # Do the actual unzipping
            unzip_todir(zip_handle, unzip_dir, 'tgz')

        # And delete the zip file + parent directories
        shutil.rmtree(download_dir)

        # Return the download location so we can pass it to the install tasks
        return unzip_dir

def nectar_download(asset_key):
    """
    Downloads the asset to a temporary directory for later installation
    """
    return {
        # Download the asset and return the directory as a doit arg
        'actions': [lambda: {'dir': download_nectar_asset(asset_key, to_temp=True)}],
        'uptodate': [False]
    }


def nectar_install(asset_key, extra_keys={}):
    """
    Downloads the asset to its final destination in the cpipe installation, not using a temporary file
    or requiring installation
    """

    def action():
        # Download the asset from nectar
        download_nectar_asset(asset_key, to_temp=False)

        # This is an installation, so update the manifest to reflect that
        add_to_manifest(asset_key)

    # Allow the user to add extra keys to the task dictionary
    task = {
        # Download the asset and return the directory as a doit arg
        'actions': [action],
        'uptodate': [not nectar_asset_needs_update(asset_key)]
    }
    task.update(extra_keys)
    return task
