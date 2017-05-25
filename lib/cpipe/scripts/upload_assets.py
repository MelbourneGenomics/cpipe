#!/usr/bin/python

# Imports
import json
import hashlib
import argparse
from os import path, remove
from subprocess import check_call
from swiftclient.service import SwiftService, SwiftUploadObject

# Argument parsing
def valid_dir(dir):
    if not path.isdir(dir):
        raise argparse.ArgumentTypeError('{0} must be an existing directory to upload'.format(dir))
    return dir

def get_args():
    parser = argparse.ArgumentParser(description='Uploads a given directory to NECTAR and updates the manifest. Note that you must have sourced the Swift Credentials file in order to run this')
    parser.add_argument('directory', type=valid_dir, help='The directory to upload as an asset')
    parser.add_argument('--version', '-v', required=True, help='The version of the asset that is being uploaded ("1.0", "b215eb", "v0.3" etc.)')
    parser.add_argument('--name', '-n', required=True, action='store', help='What to name the asset on nectar. ".tar.gz" will be automatically added')
    parser.add_argument('--manifest', '-m', required=True, action='store', help='Which file to write the new file information to')
    parser.add_argument('--container', '-c', required=False, default='cpipe-2.4-public', action='store', help='The nectar container to upload to')
    return parser.parse_args()

def main():

    # Variables
    current_dir = path.dirname(__file__)
    args = get_args()
    manifest_path = args.manifest
    object_name = args.name # Name of the file in NECTAR, without the file extension
    object_name_ext = '{}/{}.tar.gz'.format(object_name, args.version) # Name of the file in NECTAR, with the file extension
    human_name = path.basename(object_name) # Key for the asset in the manifest
    zip_file = path.basename(object_name_ext) # Zip file that will be created
    target_dir = args.directory # Directory that will be zipped

    # Main script

    # Read in the current manifest
    with open(manifest_path, 'r') as manifest_file, SwiftService() as swift:

        try:
            manifest = json.load(manifest_file)
        except:
            manifest = {}

        # Check this hasn't already been done
        if human_name in manifest and manifest[human_name]['version'] == args.version:
            raise ValueError('An asset named {0} version {1} has already been uploaded and added to the manifest'.format(human_name, args.version))

        # Zip the file
        check_call(["tar", "-zcf", path.basename(zip_file), "-C", path.dirname(target_dir), path.basename(target_dir)])

        # Hash it
        with open(zip_file, 'rb') as zip_handle:
            sha = hashlib.sha1(zip_handle.read()).hexdigest()

        # Upload it
        for result in swift.upload(
                args.container,
                [SwiftUploadObject(zip_file, object_name=object_name_ext)],
                options={'segment_size': 5368709120}
        ):
            if result['action'] != 'upload_object':
                continue

            if result['success']:
                print('{} successfully uploaded.'.format(result['object']))
            else:
                print('{} upload failed! {}'.format(result['object'], result['error']))

        # Update the manifest
        manifest[human_name] = {
            'path': object_name,
            'hash': sha,
            'version': args.version
        }

    # Write out the new manifest
    with open(manifest_path, 'w') as manifest_file:
        json.dump(manifest, manifest_file, indent=4)

        # Delete the zip file
        remove(zip_file)

if __name__ == '__main__':
    main()
